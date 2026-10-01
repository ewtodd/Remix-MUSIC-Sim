#include "Simulator.hpp"

#include <filesystem>
#include <map>
#include <memory>
#include <set>
#include <system_error>

#include <TObjString.h>

Simulator::Simulator(Int_t workerId) {
  Name = "Simulator";
  this->workerId_ = workerId;
  // Master (workerId 0) is chatty; MT workers stay quiet.
  this->verbose_ = (workerId == 0);

  CMEMax = CMEMin = 0;
  Kb_after_window = 0;
  NTraces = 0;
  PrintLevel = 0;
  tracesCreated = false;
  energeticsWritten_ = false;

  Beam = Target = Compound = 0;
  Trace = TraceUB = TraceB = nullptr;
  TraceER = TraceEP = nullptr;
  // Trace-background histograms only exist for the interactive visualizer
  // (ctf.Update != 0). Workers leave these null so they don't collide on
  // gROOT's name registry.
  HCT = HCTB = HPT = 0;

  maxEvaporations = 20;
  numEvaporations = 0;
  EvaP = new Particle *[maxEvaporations]();
  EvaR = new Particle *[maxEvaporations]();
  evap_energy = new Float_t[maxEvaporations];
  residue_energy = new Float_t[maxEvaporations];
  evap_energy_exit = new Float_t[maxEvaporations];
  residue_energy_exit = new Float_t[maxEvaporations];
  theta_cm = new Float_t[maxEvaporations];
  phi_cm = new Float_t[maxEvaporations];
  evap_theta = new Float_t[maxEvaporations];
  evap_phi = new Float_t[maxEvaporations];
  residue_theta = new Float_t[maxEvaporations];
  residue_phi = new Float_t[maxEvaporations];
  evap_stop_x = new Float_t[maxEvaporations];
  evap_stop_y = new Float_t[maxEvaporations];
  evap_stop_z = new Float_t[maxEvaporations];
  evap_stop_strip = new Int_t[maxEvaporations];
  evap_termination = new Int_t[maxEvaporations];
  minEx = new Double_t[maxEvaporations];
  for (Int_t er = 0; er < maxEvaporations; er++)
    minEx[er] = 0.0;

  NuF = new NuclideFinder();

  // TGeoManager is built only on the master; workers do anode-cell lookups
  // directly from AnodeDX/AnodeDZ and never touch ROOT's geometry singletons.
  Geo = nullptr;
  MatVacuum = nullptr;
  Vacuum = nullptr;
  VolTop = nullptr;
  if (workerId == 0) {
    Geo = new TGeoManager("Geo", "MUSIC geometry manager");
    MatVacuum = new TGeoMaterial("Vac", 0, 0, 0);
    Vacuum = new TGeoMedium("Vacuum", 1, MatVacuum);
    VolTop = Geo->MakeBox("VolTop", Vacuum, 300., 300., 300.);
    Geo->SetTopVolume(VolTop);
  }
  VolAnode = 0;

  AnodeRows = AnodeCols = 0;
  AnodeDepth = AnodeLength = AnodeHeight = 0.0;
  AnodeColor = nullptr;
  AnodeDX = AnodeDY = AnodeDZ = nullptr;
  AnodeSegName = nullptr;
  AnodeStpID = nullptr;
  DeltaEB_ave = DeltaEB = nullptr;
  DeltaE_EvaP = DeltaE_EvaR = nullptr;

  TrackBeam = nullptr;
  TrackEvaP = TrackEvaR = nullptr;
  Eve = nullptr;
  TraceCan = nullptr;
  LegCol = LegPart = nullptr;
  LabelKine = nullptr;
  TopNode = nullptr;

  Rdm = new TRandom3(1);

  SimTree = 0;
  MCTree = 0;
  ResetBranches();

  InitCTF();

  gSystem = ::gSystem;
}

Simulator::~Simulator() {
  for (Int_t i = 0; i < maxEvaporations; ++i) {
    delete EvaP[i];
    delete EvaR[i];
  }
  delete[] EvaP;
  delete[] EvaR;
  delete Beam;
  delete Target;
  delete Compound;
  delete NuF;
  delete Rdm;

  delete[] evap_energy;
  delete[] residue_energy;
  delete[] evap_energy_exit;
  delete[] residue_energy_exit;
  delete[] theta_cm;
  delete[] phi_cm;
  delete[] evap_theta;
  delete[] evap_phi;
  delete[] residue_theta;
  delete[] residue_phi;
  delete[] evap_stop_x;
  delete[] evap_stop_y;
  delete[] evap_stop_z;
  delete[] evap_stop_strip;
  delete[] evap_termination;
  delete[] minEx;

  auto freeRows = [&](auto **&rows, Int_t count) {
    if (!rows)
      return;
    for (Int_t row = 0; row < count; ++row)
      delete[] rows[row];
    delete[] rows;
    rows = nullptr;
  };
  if (DeltaE_EvaP) {
    for (Int_t er = 0; er < maxEvaporations; ++er)
      freeRows(DeltaE_EvaP[er], AnodeRows);
    delete[] DeltaE_EvaP;
  }
  if (DeltaE_EvaR) {
    for (Int_t er = 0; er < maxEvaporations; ++er)
      freeRows(DeltaE_EvaR[er], AnodeRows);
    delete[] DeltaE_EvaR;
  }
  freeRows(DeltaEB_ave, AnodeRows);
  freeRows(DeltaEB, AnodeRows);
  freeRows(AnodeColor, AnodeRows);
  freeRows(AnodeDX, AnodeRows);
  freeRows(AnodeDY, AnodeRows);
  freeRows(AnodeDZ, AnodeRows);
  freeRows(AnodeSegName, AnodeRows);
  freeRows(AnodeStpID, AnodeRows);
  freeRows(VolAnode, AnodeRows);

  auto freeGraphs = [&](TGraph **&graphs) {
    if (!graphs)
      return;
    for (Int_t col = 0; col < AnodeCols + 1; ++col)
      delete graphs[col];
    delete[] graphs;
    graphs = nullptr;
  };
  freeGraphs(Trace);
  freeGraphs(TraceUB);
  freeGraphs(TraceB);
  if (TraceER) {
    for (Int_t er = 0; er < numEvaporations; ++er)
      freeGraphs(TraceER[er]);
    delete[] TraceER;
  }
  if (TraceEP) {
    for (Int_t er = 0; er < numEvaporations; ++er)
      freeGraphs(TraceEP[er]);
    delete[] TraceEP;
  }
  // TEveManager owns the arrow objects; these two arrays only store borrowed
  // pointers to them.
  delete[] TrackEvaP;
  delete[] TrackEvaR;
  // TGeoManager registers itself in ROOT's global geometry list. ROOT owns and
  // destroys registered managers during TApplication teardown; deleting it
  // here after that teardown is a double free. Worker instances never create
  // one, and the master manager follows ROOT's process-lifetime ownership.
  Geo = nullptr;
}

Int_t Simulator::CheckMemoryUsage(Int_t Print) {
  if (gSystem == nullptr)
    return 1;
  MemInfo_t mem;
  gSystem->GetMemInfo(&mem);
  const Double_t memoryLimit = 0.95 * static_cast<Double_t>(mem.fMemTotal);
  if (Print) {
    std::cout << "> Total memory used: " << mem.fMemUsed
              << " MB = " << 100.0 * mem.fMemUsed / mem.fMemTotal
              << " % of max memory (" << mem.fMemTotal << " MB)" << std::endl;
  }
  if (static_cast<Double_t>(mem.fMemUsed) > memoryLimit) {
    std::cout << "Memory limit exceeded! Limit at " << memoryLimit << " MB"
              << std::endl;
    return 0;
  }
  return 1;
}

void Simulator::SeedRandom(UInt_t seed) {
  delete Rdm;
  Rdm = new TRandom3(seed);
}

void Simulator::SetPrintLevel(Int_t Level /*0-2*/) {
  if (Level < 0 || Level > 2) {
    std::cout
        << "Warning: Invalid information printing level.\n"
        << "Valid options are:\n"
        << "\t0 - minimum printing\n"
        << "\t1 - info per event (e.g. propagator initial and final conditions)\n"
        << "\t2 - print all" << std::endl;
    PrintLevel = 0;
    return;
  }
  if (verbose_)
    std::cout << "See musicsim.log file for detailed information" << std::endl;
  PrintLevel = Level;
  const TString logName = workerId_ == 0
                              ? "musicsim.log"
                              : TString::Format("musicsim_w%d.log", workerId_);
  Log.open(logName.Data());
  Log << "================================================================================"
      << std::endl;
  Log << "|--- MUSIC simulator log file -------------------------------------------------|"
      << std::endl;
}

void Simulator::SetROOTSystemPointer(TSystem *system) {
  gSystem = system;
  std::cout << "gSystem = " << system << std::endl;
  CheckMemoryUsage(1);
}

// Populate catima's global DataPoint cache for every (projectile, material)
// the workers will see. This is a performance optimization; every catima call
// is still protected by CatimaMutex because its process-wide cache can evict
// entries and write again after pre-warming.
void Simulator::PreWarmCatima() {
  if (!NuF)
    NuF = new NuclideFinder();
  BuildGasMaterial();
  BuildWindows();
  BuildDegrader();
  auto warmIon = [&](Int_t A, Int_t Z) {
    if (A <= 0 || Z <= 0)
      return;
    // Warm both the mean (default, atima14) config and the straggling config:
    // EnergyLoss::BuildTables and EnergyThroughWithStraggling now query catima
    // under each. DataPoint keys on Config, so the two land in distinct cache
    // slots. Warming both avoids avoidable cache misses in worker setup; the
    // mutex remains the correctness guard if the finite cache later evicts one.
    catima::Projectile proj{Double_t(A), Double_t(Z)};
    proj.T = 100.0 / A;
    auto warm = [&](Bool_t enabled, const catima::Material &mat) {
      if (!enabled)
        return;
      std::lock_guard<std::mutex> lock(music::CatimaMutex());
      catima::calculate(proj, mat, physics_.mean);
      catima::calculate(proj, mat, physics_.straggling);
    };
    warm(gasEnabled_, gas_);
    warm(entranceWindowEnabled_, entranceWindow_);
    warm(exitWindowEnabled_, exitWindow_);
    warm(hasDegrader_, degrader_);
  };
  auto warmName = [&](const TString &name) {
    if (name.IsNull() || name == "unassigned beam" ||
        name == "unassigned target" || name == "unassigned compound" ||
        name == "unassigned res" || name == "unassigned evap")
      return;
    Int_t Z = NuF->GetZ(name.Data());
    Double_t m_u = NuF->GetMass(name.Data(), "u");
    Int_t A = Int_t(std::round(m_u));
    warmIon(A, Z);
  };
  warmName(ctf.beamName);
  warmName(ctf.target);
  warmName(ctf.compound);
  for (Int_t i = 0; i < ctf.NumEvapPart; ++i) {
    warmName(ctf.res[i]);
    warmName(ctf.evap[i]);
  }
}

Bool_t Simulator::SetupRun() {
  if (ctf.Update) {
    Eve = new TEveManager(960, 1018, kTRUE, "V");
    Eve->GetDefaultGLViewer()->SetClearColor(kWhite);
    TrackBeam = new TEveArrow();
    TrackEvaP = new TEveArrow *[maxEvaporations];
    TrackEvaR = new TEveArrow *[maxEvaporations];
    for (Int_t er = 0; er < maxEvaporations; er++) {
      TrackEvaP[er] = new TEveArrow();
      TrackEvaR[er] = new TEveArrow();
    }
    TraceCan = new TCanvas("TraceCan", "Traces", 0, 0, 960, 1018);
    TraceCan->Divide(2, 1);
    TraceCan->cd(1)->SetGrid();
    TraceCan->cd(2)->SetGrid();
    LegCol = new TLegend(0.692, 0.616, 0.826, 0.861);
    LegPart = new TLegend(0.692, 0.616, 0.826, 0.861);
    LabelKine = new TPaveText(0.152, 0.679, 0.437, 0.875, "NDC");
    if (Log.is_open())
      Log << "\tVisualization objects created." << std::endl;
  } else {
    LabelKine = 0;
  }

  BuildGasMaterial();
  if (Log.is_open())
    Log << "\tGas material configured (" << ctf.gas << ", " << ctf.pressure
        << " Torr, " << ctf.temperature << " K, density " << gas_.density()
        << " g/cm^3)." << std::endl;
  BuildWindows();
  auto unitOf = [](Bool_t byLen) { return byLen ? "um" : "mg/cm^2"; };
  if (Log.is_open())
    Log << "\tEntrance window: " << ctf.entranceMaterial << " "
        << ctf.entranceThickness << " " << unitOf(ctf.entranceByLength)
        << "; exit: " << ctf.exitMaterial << " " << ctf.exitThickness << " "
        << unitOf(ctf.exitByLength) << "." << std::endl;
  BuildDegrader();
  if (hasDegrader_ && Log.is_open())
    Log << "\tDegrader: " << ctf.degraderMaterial << " " << ctf.degraderLength
        << " " << unitOf(ctf.degraderByLength) << "." << std::endl;
  if (SetAnode(90, ctf.ELossBins, ctf.MaxELoss) == 0) {
    std::cerr << "musicsim ERROR: failed to configure the detector geometry."
              << std::endl;
    return kFALSE;
  }
  if (Log.is_open())
    Log << "\tAnode configured." << std::endl;

  SetBeamParticle(ctf.beamName, kBlack, ctf.dEdxScaleBeam);
  if (Log.is_open())
    Log << "\tBeam particle configured." << std::endl;
  // Mean beam KE at the gas surface (after entrance window). The per-event
  // chain — accelerator FWHM, degrader straggling, window straggling — is
  // sampled inside the event loop; this value is only for log output and the
  // CM-energy-range estimate.
  {
    const Double_t amu_MeV = 931.49410242;
    Int_t A_beam =
        (Beam->Mass > 0) ? Int_t(std::round(Beam->Mass / amu_MeV)) : 0;
    Double_t Eaccel = ctf.BeamEnergy;
    Double_t Eafter_degrader =
        hasDegrader_ ? EnergyOutOfMaterial(A_beam, Beam->Z, Eaccel, degrader_)
                     : Eaccel;
    Kb_at_gas = entranceWindowEnabled_
                    ? EnergyOutOfMaterial(A_beam, Beam->Z, Eafter_degrader,
                                          entranceWindow_)
                    : Eafter_degrader;
    if (verbose_) {
      std::cout << "Beam energy: " << Eaccel << " MeV at accelerator";
      if (hasDegrader_)
        std::cout << " -> " << Eafter_degrader << " MeV after degrader ("
                  << ctf.degraderMaterial << " " << ctf.degraderLength << " "
                  << unitOf(ctf.degraderByLength) << ")";
      if (entranceWindowEnabled_)
        std::cout << " -> " << Kb_at_gas << " MeV at gas surface ("
                  << ctf.entranceMaterial << " " << ctf.entranceThickness << " "
                  << unitOf(ctf.entranceByLength) << " window)";
      else
        std::cout << " -> " << Kb_at_gas
                  << " MeV at gas surface (no entrance window)";
      std::cout << std::endl;
    }
    beam_energy_accel = static_cast<Float_t>(ctf.BeamEnergy);
  }
  SetTargetParticle(ctf.target);
  if (Log.is_open())
    Log << "\tTarget particle configured." << std::endl;
  SetCompoundParticle(ctf.compound);
  if (Log.is_open())
    Log << "\tCompound particle configured." << std::endl;
  for (Int_t i = 0; i < ctf.NumEvapPart; i++)
    SetEvapResAndPart(ctf.res[i], ctf.colorRes[i], ctf.evap[i],
                      ctf.colorEvap[i], ctf.dEdxScaleRes[i],
                      ctf.dEdxScaleEvap[i]);
  if (Log.is_open())
    Log << "\tEvaporated particles and residues configured." << std::endl;

  // Minimum parent excitation needed for each remaining decay step. When a
  // residue is configured as "forced", this recursively includes the
  // excitation required to reach every later step in the chain.
  for (Int_t step = numEvaporations - 1; step >= 0; step--) {
    Double_t mb = Beam->Mass;
    Double_t mt = Target->Mass;
    Double_t ml = EvaP[step]->Mass;
    Double_t mh = EvaR[step]->Mass;
    Double_t Q0 =
        (step == 0) ? (ml + mh - mb - mt) : (ml + mh - EvaR[step - 1]->Mass);
    const Double_t childExcitation =
        (step + 1 < numEvaporations && ctf.residueExc[step] == 0)
            ? minEx[step + 1]
            : 0.0;
    minEx[step] = std::max(0.0, Q0 + childExcitation);
    if (Log.is_open()) {
      if (step == 0)
        Log << "Q0(" << Beam->Name << "+" << Target->Name << "->"
            << EvaP[step]->Name << "+" << EvaR[step]->Name << ") = " << Q0
            << " MeV\tminEx" << step << " = " << minEx[step] << std::endl;
      else
        Log << "Q0(" << EvaR[step - 1]->Name << "->" << EvaP[step]->Name << "+"
            << EvaR[step]->Name << ") = " << Q0 << " MeV\tminEx" << step
            << " = " << minEx[step] << std::endl;
    }
  }
  return kTRUE;
}

UInt_t Simulator::EventSeed(Int_t strip, ULong64_t eventIndex) const {
  // SplitMix64 gives each (master seed, strip, event index) tuple an
  // independent, deterministic TRandom3 seed. Event streams therefore do not
  // depend on worker scheduling or worker count. TRandom3 treats zero as a
  // request for nondeterministic auto-seeding, so map it to one explicitly.
  ULong64_t value = ctf.Seed;
  value ^= (static_cast<ULong64_t>(static_cast<UInt_t>(strip + 2)) << 32);
  value ^= eventIndex + 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
  value ^= value >> 31;
  UInt_t seed = static_cast<UInt_t>(value ^ (value >> 32));
  if (seed == 0)
    seed = 1;
  return seed;
}

Int_t Simulator::run() {
  // MT dispatch: only the master fans out. Workers run as single-threaded
  // (their ctf.Threads is forced to 1 before run()).
  if (ctf.Threads > 1 && workerId_ == 0)
    return runMultiThreaded();

  std::unique_ptr<TFile> rootFile;
  TFile *ROOTfile = nullptr;
  std::filesystem::path finalOutput;
  std::filesystem::path temporaryOutput;
  Bool_t outputIsTemporary = kFALSE;
  Long64_t entriesBeforeRun = 0;

  SetPrintLevel(ctf.PrintOpt);
  Log << "musicsim::run() START *********************************************"
      << std::endl;

  if (!SetupRun())
    return 0;

  const TString energeticsName =
      workerId_ == 0 ? "energetics.log"
                     : TString::Format("energetics_w%d.log", workerId_);
  EnergeticsLog.open(energeticsName.Data());
  CalculateCMEnergyRange();
  EnergeticsLog.close();

  if (!ctf.FileName.IsNull()) {
    finalOutput = std::filesystem::path(ctf.FileName.Data());
    temporaryOutput = finalOutput;
    temporaryOutput += ".tmp";
    std::error_code error;
    if (std::filesystem::exists(temporaryOutput, error)) {
      std::cerr << "musicsim ERROR: refusing to overwrite stale temporary "
                   "output '"
                << temporaryOutput.string() << "'." << std::endl;
      return 0;
    }
    const Bool_t updateExisting =
        ctf.FileOpt == "update" && std::filesystem::exists(finalOutput, error);
    if (error) {
      std::cerr << "musicsim ERROR: cannot inspect output path '"
                << finalOutput.string() << "': " << error.message()
                << std::endl;
      return 0;
    }
    if (updateExisting) {
      std::filesystem::copy_file(finalOutput, temporaryOutput,
                                 std::filesystem::copy_options::none, error);
      if (error) {
        std::cerr << "musicsim ERROR: cannot stage update output: "
                  << error.message() << std::endl;
        return 0;
      }
    }
    outputIsTemporary = kTRUE;
    rootFile.reset(TFile::Open(temporaryOutput.c_str(),
                               updateExisting ? "UPDATE" : "RECREATE"));
    ROOTfile = rootFile.get();
    if (!ROOTfile || ROOTfile->IsZombie() || !ROOTfile->IsWritable()) {
      std::cerr << "musicsim ERROR: could not create writable ROOT output '"
                << temporaryOutput.string() << "'." << std::endl;
      rootFile.reset();
      std::filesystem::remove(temporaryOutput, error);
      return 0;
    }
    if (updateExisting) {
      const auto *storedControl = dynamic_cast<TObjString *>(
          ROOTfile->Get("metadata/control_file_toml"));
      std::ifstream currentControl(ctrlFilePath_.Data(), std::ios::binary);
      std::ostringstream currentText;
      if (currentControl)
        currentText << currentControl.rdbuf();
      if (!storedControl || !currentControl ||
          storedControl->GetString() != currentText.str().c_str()) {
        std::cerr << "musicsim ERROR: update output provenance does not match "
                     "the current control file."
                  << std::endl;
        rootFile->Close();
        rootFile.reset();
        std::filesystem::remove(temporaryOutput, error);
        return 0;
      }
    }
    ROOTfile->cd();
    SimTree = InitTree(ROOTfile, updateExisting ? "update" : "recreate");
    if (!SimTree || !MCTree) {
      rootFile->Close();
      rootFile.reset();
      std::filesystem::remove(temporaryOutput, error);
      return 0;
    }
    entriesBeforeRun = SimTree->GetEntries();
    if (updateExisting && entriesBeforeRun > 0) {
      std::map<Int_t, std::set<ULong64_t>> existingIndices;
      for (Long64_t entry = 0; entry < MCTree->GetEntries(); ++entry) {
        if (MCTree->GetEntry(entry) <= 0 || requested_strip < ctf.stripFirst ||
            requested_strip > ctf.stripLast ||
            !existingIndices[requested_strip].insert(event_index).second) {
          std::cerr << "musicsim ERROR: update output has invalid or duplicate "
                       "event identities."
                    << std::endl;
          rootFile->Close();
          rootFile.reset();
          std::filesystem::remove(temporaryOutput, error);
          return 0;
        }
      }
      ULong64_t priorEventsPerStrip = 0;
      for (Int_t strip = ctf.stripFirst; strip <= ctf.stripLast; ++strip) {
        const auto found = existingIndices.find(strip);
        if (found == existingIndices.end() || found->second.empty() ||
            *found->second.begin() != 0 ||
            *found->second.rbegin() + 1 != found->second.size() ||
            (priorEventsPerStrip != 0 &&
             found->second.size() != priorEventsPerStrip)) {
          std::cerr << "musicsim ERROR: update output event indices are not a "
                       "complete, equal per-strip sequence."
                    << std::endl;
          rootFile->Close();
          rootFile.reset();
          std::filesystem::remove(temporaryOutput, error);
          return 0;
        }
        priorEventsPerStrip = found->second.size();
      }
      eventOffset_ = priorEventsPerStrip;
    }

    TDirectory *traceDirectory = ROOTfile->GetDirectory("traces");
    if (!traceDirectory)
      traceDirectory = ROOTfile->mkdir("traces");
    if (!traceDirectory || !traceDirectory->cd()) {
      std::cerr << "musicsim ERROR: could not create or select the ROOT "
                   "trace directory."
                << std::endl;
      rootFile->Close();
      rootFile.reset();
      std::filesystem::remove(temporaryOutput, error);
      return 0;
    }
    Log << "\tROOT file opened." << std::endl;
  }

  // Trace TGraphs are written per event for the EVE visualization; they're
  // wasteful in bulk MC mode where the user just wants the event tree.
  if (ROOTfile != 0 && ctf.Update) {
    CreateTracesAndTrajectories();
    Log << "\tTraces and trajectories created." << std::endl;
  }
  // stripFirst/stripLast and method=0 were validated in loadCtrlFile.
  Log << "\tStarting simulation loop ..." << std::endl;
  for (Int_t stpID = ctf.stripFirst; stpID <= ctf.stripLast; stpID++)
    Simulate(stpID, ctf.NEvents, ctf.MaxTime, ctf.SimStep, ctf.Update, ctf.Wait,
             ROOTfile);
  Log << "\tSimulation loop ended." << std::endl;

  if (ROOTfile && SimTree) {
    ROOTfile->cd();
    if (ioFailed_ || SimTree->Write("", TObject::kOverwrite) <= 0 || !MCTree ||
        MCTree->Write("", TObject::kOverwrite) <= 0 ||
        (workerId_ == 0 && !WriteRunMetadata(ROOTfile)) ||
        ROOTfile->TestBit(TFile::kWriteError)) {
      std::cerr << "musicsim ERROR: failed while writing ROOT output."
                << std::endl;
      ROOTfile->Close();
      rootFile.reset();
      std::error_code error;
      if (outputIsTemporary)
        std::filesystem::remove(temporaryOutput, error);
      return 0;
    }
    ROOTfile->Close();
    rootFile.reset();
    std::error_code error;
    {
      std::unique_ptr<TFile> check(
          TFile::Open(temporaryOutput.c_str(), "READ"));
      TTree *events = check && !check->IsZombie()
                          ? dynamic_cast<TTree *>(check->Get("events_MeV"))
                          : nullptr;
      TTree *truth = check && !check->IsZombie()
                         ? dynamic_cast<TTree *>(check->Get("MC"))
                         : nullptr;
      const Long64_t expected =
          entriesBeforeRun +
          static_cast<Long64_t>(ctf.stripLast - ctf.stripFirst + 1) *
              ctf.NEvents;
      if (!events || !truth || events->GetEntries() != expected ||
          truth->GetEntries() != expected) {
        std::cerr << "musicsim ERROR: staged ROOT output failed validation."
                  << std::endl;
        std::filesystem::remove(temporaryOutput, error);
        return 0;
      }
    }
    std::filesystem::rename(temporaryOutput, finalOutput, error);
    if (error) {
      std::cerr << "musicsim ERROR: could not publish ROOT output atomically: "
                << error.message() << std::endl;
      return 0;
    }
    Log << "\tROOT file written." << std::endl;
  }

  return 1;
}

Int_t Simulator::runMultiThreaded() {
  ROOT::EnableThreadSafety();
  PreWarmCatima();

  if (!SetupRun())
    return 0;

  EnergeticsLog.open("energetics.log");
  CalculateCMEnergyRange();
  EnergeticsLog.close();

  const Int_t nThreads = std::max(1, ctf.Threads);
  const Int_t totalEvents = ctf.NEvents;
  const TString baseOutput = ctf.FileName;
  const TString ctrlPath = ctrlFilePath_;

  if (ctrlPath.IsNull()) {
    std::cerr << "musicsim: multi-threaded mode requires a control file path."
              << std::endl;
    return 0;
  }

  // Strip ".root" so we can insert per-worker tags.
  TString baseStem = baseOutput;
  if (baseStem.EndsWith(".root"))
    baseStem.Remove(baseStem.Length() - 5);

  std::cout << "Multi-threaded run: " << nThreads << " workers x "
            << (totalEvents / nThreads) << "+ events." << std::endl;

  std::vector<TString> workerOutputs;
  std::vector<Int_t> workerEventCounts;
  std::vector<std::future<Int_t>> futures;
  const Int_t evPerWorker = totalEvents / nThreads;
  const Int_t extra = totalEvents % nThreads;

  ULong64_t eventOffset = 0;
  for (Int_t w = 0; w < nThreads; ++w) {
    Int_t slice = evPerWorker + (w < extra ? 1 : 0);
    const ULong64_t workerOffset = eventOffset;
    eventOffset += static_cast<ULong64_t>(slice);
    TString out = TString::Format("%s_w%d.root", baseStem.Data(), w);
    workerOutputs.push_back(out);
    workerEventCounts.push_back(slice);

    futures.push_back(std::async(
        std::launch::async, [ctrlPath, out, slice, workerOffset, w]() -> Int_t {
          Simulator worker(w + 1); // worker ids start at 1 (0 is master)
          // loadCtrlFile takes char* (non-const); copy into a mutable buffer.
          std::vector<char> path(ctrlPath.Data(),
                                 ctrlPath.Data() + ctrlPath.Length() + 1);
          if (worker.loadCtrlFile(path.data()) == 0)
            return 0;
          worker.OverrideNEvents(slice);
          worker.OverrideOutputFile(out);
          worker.OverrideThreads(1);
          worker.OverrideEventOffset(workerOffset);
          worker.DisableVisualization();
          return worker.run();
        }));
  }

  // Block until every worker completes. std::future::get re-throws any
  // exception that crossed the thread boundary.
  for (size_t worker = 0; worker < futures.size(); ++worker) {
    try {
      if (futures[worker].get() == 0) {
        std::cerr << "musicsim ERROR: worker " << worker << " failed."
                  << std::endl;
        return 0;
      }
    } catch (const std::exception &error) {
      std::cerr << "musicsim ERROR: worker " << worker
                << " threw an exception: " << error.what() << std::endl;
      return 0;
    }
  }

  const Long64_t selectedStrips =
      static_cast<Long64_t>(ctf.stripLast - ctf.stripFirst + 1);
  for (size_t worker = 0; worker < workerOutputs.size(); ++worker) {
    std::unique_ptr<TFile> input(TFile::Open(workerOutputs[worker], "READ"));
    TTree *events = input && !input->IsZombie()
                        ? dynamic_cast<TTree *>(input->Get("events_MeV"))
                        : nullptr;
    TTree *truth = input && !input->IsZombie()
                       ? dynamic_cast<TTree *>(input->Get("MC"))
                       : nullptr;
    const Long64_t expected = selectedStrips * workerEventCounts[worker];
    if (!events || !truth || events->GetEntries() != expected ||
        truth->GetEntries() != expected) {
      std::cerr << "musicsim ERROR: worker " << worker
                << " produced an incomplete or invalid ROOT file (expected "
                << expected << " entries)." << std::endl;
      return 0;
    }
  }

  std::cout << "Merging " << workerOutputs.size() << " worker outputs into "
            << baseOutput << " ..." << std::endl;
  const std::filesystem::path finalOutput(baseOutput.Data());
  std::filesystem::path mergeOutput = finalOutput;
  mergeOutput += ".tmp";
  std::error_code fileError;
  if (std::filesystem::exists(mergeOutput, fileError)) {
    std::cerr << "musicsim ERROR: refusing to overwrite stale temporary "
                 "output '"
              << mergeOutput.string() << "'." << std::endl;
    return 0;
  }
  TFileMerger merger(kFALSE);
  if (!merger.OutputFile(mergeOutput.c_str(), "RECREATE")) {
    std::cerr << "musicsim ERROR: could not create merge destination."
              << std::endl;
    return 0;
  }
  for (const auto &p : workerOutputs) {
    if (!merger.AddFile(p.Data())) {
      std::cerr << "musicsim ERROR: could not add worker output '" << p
                << "' to the merge." << std::endl;
      return 0;
    }
  }
  Bool_t ok = merger.Merge();
  if (!ok) {
    std::cerr << "musicsim: merge failed." << std::endl;
    return 0;
  }
  merger.CloseOutputFile();
  {
    std::unique_ptr<TFile> merged(TFile::Open(mergeOutput.c_str(), "UPDATE"));
    TTree *events = merged && !merged->IsZombie()
                        ? dynamic_cast<TTree *>(merged->Get("events_MeV"))
                        : nullptr;
    TTree *truth = merged && !merged->IsZombie()
                       ? dynamic_cast<TTree *>(merged->Get("MC"))
                       : nullptr;
    const Long64_t expected = selectedStrips * totalEvents;
    if (!events || !truth || events->GetEntries() != expected ||
        truth->GetEntries() != expected || !WriteRunMetadata(merged.get()) ||
        merged->TestBit(TFile::kWriteError)) {
      std::cerr << "musicsim ERROR: merged output validation failed."
                << std::endl;
      if (merged)
        merged->Close();
      return 0;
    }
    merged->Close();
  }
  std::filesystem::rename(mergeOutput, finalOutput, fileError);
  if (fileError) {
    std::cerr << "musicsim ERROR: could not publish merged output atomically: "
              << fileError.message() << std::endl;
    return 0;
  }
  for (const auto &p : workerOutputs)
    std::remove(p.Data());
  std::cout << "Multi-threaded run complete." << std::endl;
  // Report success; main.cpp separately decides whether visualization needs the
  // interactive ROOT event loop.
  return 1;
}

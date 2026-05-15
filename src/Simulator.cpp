#include "Simulator.hpp"

Simulator::Simulator(Int_t workerId) {
  Name = "Simulator";
  this->workerId_ = workerId;
  // Master (workerId 0) is chatty; MT workers stay quiet.
  this->verbose_ = (workerId == 0);

  CMEMax = CMEMin = 0;
  EexcMax = EexcMin = 0;
  Kb_after_window = 0;
  NTraces = 0;
  PrintLevel = 0;
  SegLength = SegCMERange = SegEexcRange = 0;
  tracesCreated = false;

  Beam = Target = Compound = Light = Heavy = DeDau1 = DeDau2 = 0;
  Trace = 0;
  // Trace-background histograms only exist for the interactive visualizer
  // (ctf.Update != 0). Workers leave these null so they don't collide on
  // gROOT's name registry.
  HCT = HCTB = HPT = 0;

  maxEvaporations = 20;
  numEvaporations = 0;
  EvaP = new Particle *[maxEvaporations];
  EvaR = new Particle *[maxEvaporations];
  Kl = new Float_t[maxEvaporations];
  Kh = new Float_t[maxEvaporations];
  Kl_exit = new Float_t[maxEvaporations];
  Kh_exit = new Float_t[maxEvaporations];
  theta_CM = new Float_t[maxEvaporations];
  phi_CM = new Float_t[maxEvaporations];
  theta_l = new Float_t[maxEvaporations];
  phi_l = new Float_t[maxEvaporations];
  theta_h = new Float_t[maxEvaporations];
  phi_h = new Float_t[maxEvaporations];
  xfl = new Float_t[maxEvaporations];
  yfl = new Float_t[maxEvaporations];
  zfl = new Float_t[maxEvaporations];
  minEx = new Double_t[maxEvaporations];
  for (Int_t er = 0; er < maxEvaporations; er++) {
    Kl[er] = Kh[er] = -2.0f;
    phi_CM[er] = theta_CM[er] = -1;
    phi_l[er] = theta_l[er] = -1;
    phi_h[er] = theta_h[er] = -1;
    xfl[er] = yfl[er] = 0;
    zfl[er] = -1000;
    minEx[0] = 0.0;
  }
  xr = yr = 0;
  zr = -1000;
  xfe = yfe = 0;
  zfe = -1000;
  resID = -1;

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

  // TRandom3(0) reseeds itself per process, so runs are not deterministic
  // by default — the MT driver overrides this via SeedRandom.
  Rdm = new TRandom3(0);

  SimTree = 0;
  MCTree = 0;
  for (Int_t s = 0; s < N_STRIPS; ++s) {
    LeftdE[s] = RightdE[s] = TotaldE[s] = 0;
  }
  Cathode = 0;

  InitCTF();

  gSystem = 0;
}

Int_t Simulator::CheckMemoryUsage(Int_t Print) {
  if (gSystem == nullptr)
    return 1;
  MemInfo_t mem;
  gSystem->GetMemInfo(&mem);
  Float_t MemoryLimit = 0.95 * mem.fMemTotal;
  if (Print) {
    std::cout << "> Total memory used: " << mem.fMemUsed
              << " MB = " << 100.0 * mem.fMemUsed / mem.fMemTotal
              << " % of max memory (" << mem.fMemTotal << " MB)" << std::endl;
  }
  if (static_cast<Float_t>(mem.fMemUsed) > MemoryLimit) {
    std::cout << "Memory limit exceeded! Limit at " << MemoryLimit << " MB"
              << std::endl;
    return 0;
  }
  return 1;
}

void Simulator::SeedRandom(ULong_t s) {
  delete Rdm;
  Rdm = new TRandom3(s);
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
  Log.open("musicsim.log");
  Log << "================================================================================"
      << std::endl;
  Log << "|--- MUSIC simulator log file -------------------------------------------------|"
      << std::endl;
}

void Simulator::SetROOTSystemPointer(TSystem *gSystem) {
  this->gSystem = gSystem;
  std::cout << "gSystem = " << gSystem << std::endl;
  CheckMemoryUsage(1);
}

// Populate catima's global DataPoint cache for every (projectile, material)
// the workers will see. catima's cache is a process-wide ring buffer;
// concurrent writes race, but read-only lookups (cspline_special, our default
// eval) are thread-safe. Pre-warming on the master makes every worker call a
// pure read.
void Simulator::PreWarmCatima() {
  if (!NuF)
    NuF = new NuclideFinder();
  BuildGasMaterial();
  BuildWindows();
  BuildDegrader();
  auto warmIon = [&](Int_t A, Int_t Z) {
    if (A <= 0 || Z <= 0)
      return;
    if (gasEnabled_)
      EnergyOutOfMaterial(A, Z, 100.0, gas_);
    if (entranceWindowEnabled_)
      EnergyOutOfMaterial(A, Z, 100.0, entranceWindow_);
    if (exitWindowEnabled_)
      EnergyOutOfMaterial(A, Z, 100.0, exitWindow_);
    if (hasDegrader_)
      EnergyOutOfMaterial(A, Z, 100.0, degrader_);
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

Int_t Simulator::run() {
  // MT dispatch: only the master fans out. Workers run as single-threaded
  // (their ctf.Threads is forced to 1 before run()).
  if (ctf.Threads > 1 && workerId_ == 0)
    return runMultiThreaded();

  TFile *ROOTfile = 0;

  SetPrintLevel(ctf.PrintOpt);
  Log << "musicsim::run() START *********************************************"
      << std::endl;

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
    Log << "\tVisualization objects created." << std::endl;
  } else {
    LabelKine = 0;
  }

  if (ctf.IgnoreShortStrips) {
    Log << "\tIgnoring short (4 cm) half-strip electrodes on the anode readout."
        << " They still contribute to the cathode sum." << std::endl;
    if (verbose_)
      std::cerr << "musicsim WARNING: ignore_short_strips=true — short (4 cm) "
                   "half-strip dE is zeroed on Left/Right anode branches, but "
                   "still summed into Cathode."
                << std::endl;
  }
  BuildGasMaterial();
  Log << "\tGas material configured (" << ctf.gas << ", " << ctf.pressure
      << " Torr, " << ctf.temperature << " K, density " << gas_.density()
      << " g/cm^3)." << std::endl;
  BuildWindows();
  auto unitOf = [](Bool_t byLen) { return byLen ? "um" : "mg/cm^2"; };
  Log << "\tEntrance window: " << ctf.entranceMaterial << " "
      << ctf.entranceThickness << " " << unitOf(ctf.entranceByLength)
      << "; exit: " << ctf.exitMaterial << " " << ctf.exitThickness << " "
      << unitOf(ctf.exitByLength) << "." << std::endl;
  BuildDegrader();
  if (hasDegrader_)
    Log << "\tDegrader: " << ctf.degraderMaterial << " " << ctf.degraderLength
        << " " << unitOf(ctf.degraderByLength) << "." << std::endl;
  if (SetAnode(90, ctf.ELossBins, ctf.MaxELoss) == 0)
    std::exit(EXIT_FAILURE);
  Log << "\tAnode configured." << std::endl;

  SetBeamParticle(ctf.beamName, kBlack, ctf.dEdxScaleBeam);
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
    BeamEnergyAccel = ctf.BeamEnergy;
  }
  SetTargetParticle(ctf.target);
  Log << "\tTarget particle configured." << std::endl;
  SetCompoundParticle(ctf.compound);
  Log << "\tCompound particle configured." << std::endl;
  for (Int_t i = 0; i < ctf.NumEvapPart; i++)
    SetEvapResAndPart(ctf.res[i], ctf.colorRes[i], ctf.evap[i],
                      ctf.colorEvap[i], ctf.dEdxScaleRes[i],
                      ctf.dEdxScaleEvap[i]);
  Log << "\tEvaporated particles and residues configured." << std::endl;

  // Minimum excitation energies needed for each step of the reaction/decay
  // chain to be energetically allowed.
  for (Int_t step = numEvaporations - 1; step >= 0; step--) {
    Double_t mb = Beam->Mass;
    Double_t mt = Target->Mass;
    Double_t ml = EvaP[step]->Mass;
    Double_t mh = EvaR[step]->Mass;
    Double_t Q0 =
        (step == 0) ? (ml + mh - mb - mt) : (ml + mh - EvaR[step - 1]->Mass);
    minEx[step] = (Q0 < 0) ? -Q0 : 0;
    if (step == 0)
      Log << "Q0(" << Beam->Name << "+" << Target->Name << "->"
          << EvaP[step]->Name << "+" << EvaR[step]->Name << ") = " << Q0
          << " MeV\tminEx" << step << " = " << minEx[step] << std::endl;
    else
      Log << "Q0(" << EvaR[step - 1]->Name << "->" << EvaP[step]->Name << "+"
          << EvaR[step]->Name << ") = " << Q0 << " MeV\tminEx" << step << " = "
          << minEx[step] << std::endl;
  }

  if (!ctf.FileName.IsNull()) {
    ROOTfile = new TFile(ctf.FileName.Data(), ctf.FileOpt.Data());
    SimTree = InitTree(ROOTfile, ctf.FileOpt);

    TDirectory *trace_dir = 0;
    if (ROOTfile) {
      if (ctf.FileOpt == "update" || ctf.FileOpt == "UPDATE") {
        trace_dir = (TDirectory *)ROOTfile->Get("traces");
        trace_dir->cd();
      } else {
        trace_dir = ROOTfile->mkdir("traces");
        trace_dir->cd();
      }
    }
    Log << "\tROOT file opened." << std::endl;
  }

  if (ctf.Method == 0) {
    // Trace TGraphs are written per event for the EVE visualization; they're
    // wasteful in bulk MC mode where the user just wants the event tree.
    if (ROOTfile != 0 && ctf.Update) {
      CreateTracesAndTrajectories();
      Log << "\tTraces and trajectories created." << std::endl;
    }
    // stripFirst/stripLast were validated and resolved in loadCtrlFile.
    Log << "\tStarting simulation loop ..." << std::endl;
    for (Int_t stpID = ctf.stripFirst; stpID <= ctf.stripLast; stpID++)
      Simulate(stpID, ctf.NEvents, ctf.MaxTime, ctf.SimStep, ctf.Update,
               ctf.Wait, ROOTfile);
    Log << "\tSimulation loop ended." << std::endl;
  } else if (ctf.Method == 1) {
    std::cout << "musicsim warning: GenerateTraceDataBase method not ready."
              << std::endl;
  }

  if (ROOTfile && SimTree) {
    ROOTfile->cd();
    SimTree->Write("", TObject::kSingleKey);
    if (MCTree)
      MCTree->Write("", TObject::kSingleKey);
    ROOTfile->Close();
    Log << "\tROOT file written." << std::endl;
  }

  return ctf.Update;
}

Int_t Simulator::runMultiThreaded() {
  ROOT::EnableThreadSafety();
  PreWarmCatima();

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
  std::vector<std::future<Int_t>> futures;
  const Int_t evPerWorker = totalEvents / nThreads;
  const Int_t extra = totalEvents % nThreads;

  for (Int_t w = 0; w < nThreads; ++w) {
    Int_t slice = evPerWorker + (w < extra ? 1 : 0);
    TString out = TString::Format("%s_w%d.root", baseStem.Data(), w);
    workerOutputs.push_back(out);

    futures.push_back(
        std::async(std::launch::async, [ctrlPath, out, slice, w]() -> Int_t {
          Simulator worker(w + 1); // worker ids start at 1 (0 is master)
          // loadCtrlFile takes char* (non-const); copy into a mutable buffer.
          std::vector<char> path(ctrlPath.Data(),
                                 ctrlPath.Data() + ctrlPath.Length() + 1);
          if (worker.loadCtrlFile(path.data()) == 0)
            return 0;
          worker.OverrideNEvents(slice);
          worker.OverrideOutputFile(out);
          worker.OverrideThreads(1);
          worker.DisableVisualization();
          worker.SeedRandom(0xC1A55EEDULL + ULong64_t(w) * 0xDEADBEEFULL);
          worker.run();
          return 1;
        }));
  }

  // Block until every worker completes. std::future::get re-throws any
  // exception that crossed the thread boundary.
  for (auto &f : futures)
    f.get();

  std::cout << "Merging " << workerOutputs.size() << " worker outputs into "
            << baseOutput << " ..." << std::endl;
  TFileMerger merger(kFALSE);
  merger.OutputFile(baseOutput.Data(), "RECREATE");
  for (const auto &p : workerOutputs)
    merger.AddFile(p.Data());
  Bool_t ok = merger.Merge();
  if (!ok) {
    std::cerr << "musicsim: merge failed." << std::endl;
    return 0;
  }
  for (const auto &p : workerOutputs)
    std::remove(p.Data());
  std::cout << "Multi-threaded run complete." << std::endl;
  // Return 0 so main.cpp doesn't start the interactive ROOT event loop.
  return 0;
}

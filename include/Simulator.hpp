#ifndef SIMULATOR_HPP
#define SIMULATOR_HPP

#include <algorithm>
#include <array>
#include <cstdlib>
#include <fstream>
#include <future>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TEveArrow.h>
#include <TEveGeoNode.h>
#include <TEveManager.h>
#include <TEveStraightLineSet.h>
#include <TFile.h>
#include <TFileMerger.h>
#include <TGLViewer.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGraph.h>
#include <TH2.h>
#include <TLegend.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TStopwatch.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <toml++/toml.hpp>

#include "catima/catima.h"
#include "catima/config.h"

#include "NuclideFinder.hpp"
#include "Particle.hpp"

class Simulator {
public:
  Simulator(Int_t workerId = 0);
  ~Simulator();

  Simulator(const Simulator &) = delete;
  Simulator &operator=(const Simulator &) = delete;

  Int_t loadCtrlFile(const char *fileName);
  Int_t run();
  Bool_t WantsVisualization() const { return ctf.Update != 0; }

  // Used by the MT driver to repurpose a per-worker instance.
  void OverrideNEvents(Int_t n) { ctf.NEvents = n; }
  void OverrideOutputFile(const TString &f) { ctf.FileName = f; }
  void OverrideThreads(Int_t t) { ctf.Threads = t; }
  void OverrideEventOffset(ULong64_t offset) { eventOffset_ = offset; }
  void DisableVisualization() {
    ctf.Update = 0;
    ctf.Wait = 0;
  }
  void SeedRandom(UInt_t seed);

  void CalculateCMEnergyRange();
  Int_t SetAnode(Short_t Trans /*0..100*/, Int_t ELossBins = 400,
                 Float_t MaxELoss = 5);
  void SetBeamParticle(TString particleName, Int_t color,
                       Float_t dEdxScale = 1.0);
  void SetCompoundParticle(TString particleName);
  void SetEvapResAndPart(TString residueName, Int_t residueColor,
                         TString particleName, Int_t particleColor,
                         Float_t dEdxScaleRes = 1.0,
                         Float_t dEdxScalePar = 1.0);
  void SetPrintLevel(Int_t PrintLevel);
  void SetROOTSystemPointer(TSystem *system);
  void SetTargetParticle(TString particleName);
  void Simulate(Int_t StpID, Int_t eventCount, Double_t MaxTime,
                Double_t UserStep, Int_t UpdateVis = 0, Int_t Wait = 0,
                TFile *ROOTfile = 0);

private:
  // Materials and energy-loss propagation through them.
  void BuildGasMaterial();
  void BuildWindows();
  void BuildDegrader();
  catima::Material LookupMaterial(const TString &name);
  catima::Material BuildSolidMaterial(const TString &name,
                                      Double_t thickness_mg_per_cm2);
  catima::Material BuildBulkMaterial(const TString &name,
                                     Double_t thickness_um);
  Double_t EnergyOutOfMaterial(Int_t A, Int_t Z, Double_t Ein_MeV,
                               const catima::Material &mat);
  Double_t EnergyThroughWithStraggling(Int_t A, Int_t Z, Double_t Ein_MeV,
                                       const catima::Material &mat);
  void PreWarmCatima();

  // Geometry.
  void LoadHardcodedAnodeGeometry();
  void DrawMUSIC(TEveManager *eveManager, Short_t Transparency /*0..100*/);
  // Map a z position (cm along the beam axis) to the readout strip at that
  // depth: 0..17, or -1 for the dead layers and z outside the active volume.
  Int_t StripAtZ(Double_t z);

  // Event lifecycle.
  void FinalizeEvent(Int_t eventIndex);
  void ComputeExitEnergies();
  void ComputeDetectorResponse(Int_t event, Int_t reacStp, Int_t UpdateVis);
  enum class TerminationReason : Int_t {
    NotPropagated = 0,
    ReachedVertex = 1,
    Stopped = 2,
    DownstreamExit = 3,
    UpstreamExit = 4,
    SideExit = 5,
    Timeout = 6,
    InvalidState = 7,
  };
  struct PropagationResult {
    TerminationReason reason = TerminationReason::InvalidState;
    Double_t distance = 0.0;
    Double_t depositedEnergy = 0.0;
    Bool_t reachedVertex() const {
      return reason == TerminationReason::ReachedVertex;
    }
  };
  // Propagate `PO` step-by-step through the gas, accumulating per-strip energy
  // deposits in `DE`. If `endZ > 0`, the loop exits when the particle crosses
  // that z (forward sweep up to a reaction vertex); otherwise the bound is the
  // full anode depth. If `reset_DE == false` the caller-provided `DE` array is
  // *added to* rather than zeroed first — used to continue an interrupted
  // sweep past a rejected reaction vertex.
  PropagationResult PropagateParticle(Particle *PO, Int_t Event,
                                      Double_t MaxTime, Double_t UserStep,
                                      Double_t **DE, Double_t endZ = -1.0,
                                      Bool_t reset_DE = true);

  // I/O and visualization.
  void CreateTracesAndTrajectories();
  TTree *InitTree(TFile *ROOTfile, TString FileOpt);
  void ResetBranches();
  Bool_t WriteRunMetadata(TFile *ROOTfile);
  void UpdateVisuals(Int_t event, Double_t Kbr, Double_t zr, Double_t TOF,
                     Int_t Wait = 0);

  // Kinematics.
  void SetInitialKinematics(Double_t Kbi);
  Int_t SetReactionKinematics(Double_t Kbr, Double_t zr, Double_t tof,
                              Double_t theta_CM = -1, Double_t phi_CM = -1);
  void PrintCompoundEexc(Double_t Kb, Double_t **DeltaEB);

  // Lifecycle / control.
  void InitCTF();
  Bool_t SetupRun();
  Int_t CheckMemoryUsage(Int_t Print = 0);
  Int_t runMultiThreaded();
  UInt_t EventSeed(Int_t strip, ULong64_t eventIndex) const;

  TRandom *Rdm;

  // Particles.
  Particle *Beam;
  Particle *Target;
  Particle *Compound;
  Particle **EvaP;
  Particle **EvaR;
  Int_t maxEvaporations;
  Int_t numEvaporations;
  Double_t Kb_after_window;
  Double_t *minEx;
  music::PhysicsConfig physics_;

  // Map (stpid, col) to flat index in per-electrode configuration arrays.
  // User-facing order: S0, L1, R1, L2, R2, ..., L16, R16, S17 (34 entries).
  // In code: col 0 = beam right, col 1 = beam left.
  static constexpr Int_t kNumElectrodes = 34;
  static Int_t ElectrodeIndex(Int_t stpid, Int_t col) {
    if (stpid == 0)
      return 0;
    if (stpid == 17)
      return 33;
    return 2 * stpid - 1 + (col == 0 ? 1 : 0);
  }

  // Per-strip energy deposits.
  Double_t **DeltaEB_ave;
  Double_t **DeltaEB;
  Double_t ***DeltaE_EvaP;
  Double_t ***DeltaE_EvaR;

  TString Name;
  Int_t PrintLevel;

  Double_t CMEMax;
  Double_t CMEMin;

  TH2F *HCT;
  TH2F *HCTB;
  TH2F *HPT;
  TGraph **Trace;
  TGraph **TraceUB;
  TGraph **TraceB;
  TGraph ***TraceER;
  TGraph ***TraceEP;
  Bool_t tracesCreated;
  Bool_t energeticsWritten_;

  Int_t NTraces;
  Int_t NEvents;

  TEveArrow *TrackBeam;
  TEveArrow **TrackEvaP;
  TEveArrow **TrackEvaR;
  TEveManager *Eve;
  TCanvas *TraceCan;
  TLegend *LegCol;
  TLegend *LegPart;
  TPaveText *LabelKine;

  // Geometry.
  Int_t AnodeRows;
  Int_t AnodeCols;
  Double_t AnodeDepth;
  Double_t AnodeLength;
  Double_t AnodeHeight;
  Short_t **AnodeColor;
  Double_t **AnodeDX;
  Double_t **AnodeDY;
  Double_t **AnodeDZ;
  TString **AnodeSegName;
  Int_t **AnodeStpID;
  TGeoManager *Geo;
  TGeoMaterial *MatVacuum;
  TGeoMedium *Vacuum;
  TGeoVolume ***VolAnode;
  TGeoVolume *VolTop;
  TEveGeoTopNode *TopNode;

  NuclideFinder *NuF;

  // Stopping-power materials (catima).
  catima::Material gas_;
  catima::Material entranceWindow_;
  catima::Material exitWindow_;
  catima::Material degrader_;
  Bool_t hasDegrader_ = false;
  // Per-layer enable flags. Disabled = skipped in straggling chain,
  // ComputeExitEnergies, and PreWarmCatima. Triggered by Pressure <= 0 (gas)
  // or Thickness <= 0 (windows).
  Bool_t gasEnabled_ = true;
  Bool_t entranceWindowEnabled_ = true;
  Bool_t exitWindowEnabled_ = true;
  // Beam KE at the gas surface, after entrance window. Derived from
  // ctf.BeamEnergy.
  Double_t Kb_at_gas;

  // TTree output. Layout mirrors the upstream EventBuilderNearestGrid format
  // but with MeV energies rather than ADC counts (hence the "_MeV" suffix).
  // "events_MeV" carries detector-level branches; "MC" carries truth-only
  // branches and is friended onto events_MeV (rows correspond 1:1).
  TTree *SimTree;
  TTree *MCTree;
  static const Int_t N_SEG_STRIPS = 16;
  // Mirrors the experimental "events" tree: the segmented strips 1..16 are read
  // at a left and a right end, held in arrays of 16 indexed by strip - 1; the
  // unsegmented strips 0 and 17 are each one scalar. The strip total is L+R
  // for 1..16; no TotaldE branch, matching the data tree.
  Float_t LeftdE[N_SEG_STRIPS];
  Float_t RightdE[N_SEG_STRIPS];
  Float_t Strip0dE;
  Float_t Strip17dE;
  Float_t Cathode;

  // MC truth branches (live on MCTree).
  // Energy sentinels: -1.0 = transported but did not reach the downstream
  // exit, -2.0 = N/A for this event (e.g. beam_energy_exit on a reacted event,
  // evap slots of a disallowed step). The termination branches distinguish a
  // stop, upstream/side escape, timeout, and invalid state.
  // *_stop_strip: 0..17 = stopped in that readout strip, -1 = did not stop
  // in a readout strip (exited the gas or stopped in a dead layer —
  // disambiguate via *_stop_z and the exit-energy branch), -2 = N/A.
  // *_stop_x/y/z hold the stop point, or the point where the particle left
  // the active volume if it exited.
  Int_t n_steps;         // configured reaction steps; on-disk array length
  Int_t requested_strip; // configured strip for this event (-1 = beam run)
  Int_t reaction_strip;  // actual reaction strip (-1 = no reaction occurred)
  ULong64_t event_index; // per-strip event index, stable across thread counts
  Float_t beam_energy_accel;    // KE at the accelerator (= ctf.BeamEnergy)
  Float_t beam_energy_gas;      // KE at the gas surface (after entrance window)
  Float_t beam_energy_reaction; // KE at the reaction vertex
  Float_t beam_energy_exit;     // KE after the exit window (unreacted events)
  Float_t beam_stop_x;
  Float_t beam_stop_y;
  Float_t beam_stop_z;
  Int_t beam_stop_strip;
  Float_t vertex_x;
  Float_t vertex_y;
  Float_t vertex_z;
  Float_t DeadUS_dE; // total dE in upstream dead gas layer (not read out)
  Float_t DeadDS_dE; // total dE in downstream dead gas layer
  // Per-step arrays (evap = light ejectile, residue = heavy product).
  Float_t *evap_energy; // KE at creation
  Float_t *residue_energy;
  Float_t *evap_energy_exit;
  Float_t *residue_energy_exit; // only [residue_step] can be a real value:
                                // superseded residues decay, they never exit
  Float_t *theta_cm;            // CM emission angles [deg]
  Float_t *phi_cm;
  Float_t *evap_theta; // lab angles [deg]
  Float_t *evap_phi;
  Float_t *residue_theta;
  Float_t *residue_phi;
  Float_t *evap_stop_x;
  Float_t *evap_stop_y;
  Float_t *evap_stop_z;
  Int_t *evap_stop_strip;
  // Only the chain's final residue is transported, so its stop info is
  // scalar; residue_step says which step it came from (-1 = none).
  Float_t residue_stop_x;
  Float_t residue_stop_y;
  Float_t residue_stop_z;
  Int_t residue_stop_strip;
  Int_t residue_step;
  Int_t beam_termination;
  Int_t *evap_termination;
  Int_t residue_termination;

  TerminationReason beamTermination_ = TerminationReason::NotPropagated;
  std::array<TerminationReason, 20> evapTermination_{};
  std::array<TerminationReason, 20> residueTermination_{};
  Bool_t eventReacted_ = kFALSE;
  Bool_t ioFailed_ = kFALSE;

  std::ofstream Log;
  std::ofstream EnergeticsLog;

  TSystem *gSystem;

  // Worker id (0 = master / single-threaded). Used to name per-worker
  // TGeoManagers and per-worker output files.
  Int_t workerId_ = 0;
  // Master (workerId 0) is chatty; MT workers stay quiet so only the
  // runMultiThreaded summary lines reach stdout.
  Bool_t verbose_ = true;
  TString ctrlFilePath_;
  ULong64_t eventOffset_ = 0;

  static constexpr Double_t c = 29.9792458; // speed of light, cm/ns
  static constexpr Double_t pi = 3.14159265359;

  struct controlFileParams {
    TString gas = "4He";
    Double_t pressure = 760.0;    // Torr
    Double_t temperature = 293.0; // K
    Int_t ELossBins = 300;
    Float_t MaxELoss = 10.0;
    TString beamName = "";
    Float_t dEdxScaleBeam = 1.0;
    TString target = "";
    TString compound = "";
    Int_t NumEvapPart = 0;
    // Stopping-power model for the gas: 0 = catima (built in), 1 = SRIM tables
    // read from disk. SRIM tables are per (ion, gas, pressure, temperature) and
    // are generated once by the srim-cache tool; a missing table is an error
    // rather than a silent fallback, since the two models disagree by ~10% and
    // a quiet substitution would invalidate any dedx_scale calibrated on one.
    Int_t stoppingModel = 0;

    static const Int_t MaxNumEvapPart = 10;
    std::array<TString, MaxNumEvapPart> res;
    std::array<Float_t, MaxNumEvapPart> dEdxScaleRes{};
    std::array<Int_t, MaxNumEvapPart> colorRes{};
    std::array<TString, MaxNumEvapPart> evap;
    std::array<Float_t, MaxNumEvapPart> dEdxScaleEvap{};
    std::array<Int_t, MaxNumEvapPart> colorEvap{};
    // Per-step residue excitation: 0 = explicitly forced, 1 = ground,
    // 2 = uniform. -1 means not supplied and is rejected for reacted runs.
    std::array<Int_t, MaxNumEvapPart> residueExc{};
    // Per-step angular distribution: 0 = isotropic, 1 = Rutherford; -1 means
    // not supplied and is rejected for reacted runs.
    std::array<Int_t, MaxNumEvapPart> angularDist{};
    std::array<Double_t, MaxNumEvapPart> thetaCmMinDeg{};
    Double_t BeamEnergy = std::numeric_limits<Double_t>::quiet_NaN();
    Double_t KbFWHM = 0; // MeV FWHM at the accelerator
    TString entranceMaterial = "Ti";
    TString exitMaterial = "Ti";
    // Layer thickness, stored as either mg/cm² (areal density) or μm (linear
    // length along beam axis), keyed off the matching *ByLength flag.
    // Each layer's dedx_scale multiplies its stopping power; it is applied as
    // an equivalent thickness, since the energy lost in a layer is the
    // integral of dE/dx over its length (see BuildWindows).
    Double_t entranceThickness = -1.0;
    Bool_t entranceByLength = false;
    Double_t entranceScale = 1.0;
    Double_t exitThickness = -1.0;
    Double_t exitScale = 1.0;
    Double_t degraderScale = 1.0;
    Bool_t exitByLength = false;
    TString degraderMaterial = "";
    Double_t degraderLength = 0.0;
    Bool_t degraderByLength = true;
    // Reaction strip selection. Either 'strip' (a single value; -1 = unreacted
    // beam) OR both 'stripFirst' and 'stripLast' (a range). kStripUnset marks
    // "not set" so a missing key is caught explicitly rather than silently
    // defaulting.
    static constexpr Int_t kStripUnset = -99999;
    Int_t strip = kStripUnset;
    Int_t stripFirst = kStripUnset;
    Int_t stripLast = kStripUnset;
    // Exactly one detector-noise schema may be selected: 1 = absolute
    // Gaussian sigma [MeV], 2 = relative FWHM [%]. Values are per electrode;
    // a negative value disables noise. Scalar keys broadcast to anodes only.
    Int_t noiseMode = 0;
    std::array<Double_t, kNumElectrodes> Noise;
    Double_t NoiseCathode = -1;
    Int_t NEvents = 10; // events per selected strip
    Int_t Wait = 0;
    Int_t Update = 0;
    Double_t MaxTime = 2000.0;
    Double_t SimStep = 0.001;
    Int_t Method = 0;
    Int_t Threads = 1;
    TString FileName = "";
    TString FileOpt = "recreate";
    Int_t PrintOpt = 0;
    ULong64_t Seed = 1;

    controlFileParams() {
      Noise.fill(-1.0);
      dEdxScaleRes.fill(1.0f);
      dEdxScaleEvap.fill(1.0f);
      colorRes.fill(416);
      colorEvap.fill(616);
      residueExc.fill(-1);
      angularDist.fill(-1);
      thetaCmMinDeg.fill(1.0);
    }
  };

  controlFileParams ctf;
};

#endif

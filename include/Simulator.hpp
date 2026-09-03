#ifndef SIMULATOR_HPP
#define SIMULATOR_HPP

#include <algorithm>
#include <array>
#include <cstdlib>
#include <fstream>
#include <future>
#include <iostream>
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
  Int_t loadCtrlFile(char *fileName);
  Int_t run();

  // Used by the MT driver to repurpose a per-worker instance.
  void OverrideNEvents(Int_t n) { ctf.NEvents = n; }
  void OverrideOutputFile(const TString &f) { ctf.FileName = f; }
  void OverrideThreads(Int_t t) { ctf.Threads = t; }
  void DisableVisualization() {
    ctf.Update = 0;
    ctf.Wait = 0;
  }
  void SeedRandom(ULong_t s);

  void CalculateCMEnergyRange();
  void CalculateExcEnergyRange();
  void GenerateTraceDatabase(TString FileName, Double_t ThCMMin,
                             Double_t ThCMMax, Int_t ThSteps, Double_t PhiCMMin,
                             Double_t PhiCMMax, Int_t PhiSteps,
                             Double_t MaxTime, Double_t UserStep,
                             Int_t UpdateEnabled = 0, Int_t Wait = 0);
  Int_t SetAnode(Short_t Trans /*0..100*/, Int_t ELossBins = 400,
                 Float_t MaxELoss = 5);
  void SetBeamParticle(TString Name, Int_t Color, Float_t dEdxScale = 1.0);
  void SetCompoundParticle(TString Name);
  void SetDecayDaughter1(TString Name, Int_t Color);
  void SetDecayDaughter2(TString Name, Int_t Color);
  void SetEvapResAndPart(TString ResName, Int_t ResColor, TString ParName,
                         Int_t ParColor, Float_t dEdxScaleRes = 1.0,
                         Float_t dEdxScalePar = 1.0);
  void SetHeavyParticle(TString Name, Int_t Color, Int_t NEexc = 0,
                        Double_t *Eexc = 0 /*MeV*/);
  void SetLightParticle(TString Name, Int_t Color);
  void SetPrintLevel(Int_t PrintLevel);
  void SetROOTSystemPointer(TSystem *gSystem);
  void SetTargetParticle(TString Name);
  void Simulate(Int_t StpID, Int_t NEvents, Double_t MaxTime, Double_t UserStep,
                Int_t UpdateVis = 0, Int_t Wait = 0, TFile *ROOTfile = 0);
  void Simulate(Int_t StpID, Double_t ThCMMin, Double_t ThCMMax, Int_t ThSteps,
                Double_t PhiCMMin, Double_t PhiCMMax, Int_t PhiSteps,
                Double_t MaxTime, Double_t UserStep, Int_t Wait = 0);
  void WriteTraces(char *FileName);

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
  void DrawMUSIC(TEveManager *gEve, Short_t Transparency /*0..100*/);
  // Map a z position (cm along the beam axis) to the readout strip at that
  // depth: 0..17, or -1 for the dead layers and z outside the active volume.
  Int_t StripAtZ(Double_t z);

  // Event lifecycle.
  void FinalizeEvent(Int_t eventIndex);
  void ComputeExitEnergies();
  void ComputeDetectorResponse(Int_t event, Int_t reacStp, Int_t UpdateVis);
  // Propagate `PO` step-by-step through the gas, accumulating per-strip energy
  // deposits in `DE`. If `endZ > 0`, the loop exits when the particle crosses
  // that z (forward sweep up to a reaction vertex); otherwise the bound is the
  // full anode depth. If `reset_DE == false` the caller-provided `DE` array is
  // *added to* rather than zeroed first — used to continue an interrupted
  // sweep past a rejected reaction vertex.
  Int_t PropagateParticle(Particle *PO, Int_t Event, Double_t MaxTime,
                          Double_t UserStep, Double_t **DE,
                          Double_t endZ = -1.0, Bool_t reset_DE = true);

  // I/O and visualization.
  void CreateTracesAndTrajectories();
  TTree *InitTree(TFile *ROOTfile, TString FileOpt);
  void ResetBranches();
  void UpdateVisuals(Int_t event, Double_t Kbr, Double_t zr, Double_t TOF,
                     Int_t Wait = 0);

  // Kinematics.
  void SetInitialKinematics(Double_t Kbi);
  Int_t SetReactionKinematics(Double_t Kbr, Double_t zr, Double_t tof,
                              Double_t theta_CM = -1, Double_t phi_CM = -1);
  void PrintCompoundEexc(Double_t Kb, Double_t **DeltaEB);
  void PrintEnergetics(Double_t Kb, Double_t **DeltaEB);

  // Lifecycle / control.
  void InitCTF();
  void SetupRun();
  Int_t CheckMemoryUsage(Int_t Print = 0);
  Int_t runMultiThreaded();

  TRandom *Rdm;

  // Particles.
  Particle *Beam;
  Particle *Target;
  Particle *Compound;
  Particle *Heavy;
  Particle *DeDau1;
  Particle *DeDau2;
  Particle *Light;
  Particle **EvaP;
  Particle **EvaR;
  Int_t maxEvaporations;
  Int_t numEvaporations;
  Double_t Kb_after_window;
  Double_t *minEx;

  // Map (stpid, col) to flat index in ctf.Eres / per-electrode arrays.
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
  Double_t **DeltaEL;
  Double_t **DeltaEH;
  Double_t **DeltaED1;
  Double_t **DeltaED2;
  Double_t ***DeltaE_EvaP;
  Double_t ***DeltaE_EvaR;

  TString Name;
  Int_t PrintLevel;

  Double_t *SegLength;
  Double_t *SegEexcRange;
  Double_t *SegCMERange;
  Double_t CMEMax;
  Double_t CMEMin;
  Double_t EexcMax;
  Double_t EexcMin;

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
  static const Int_t N_STRIPS = 18;
  // Mirrors the experimental "events" tree: Left_0_17_dE[s] holds the left end
  // of strips 1..16 plus the full energy of the single-ended guard strips 0/17
  // (RightdE is 0 there). The strip total is L+R everywhere; no TotaldE branch,
  // matching the data tree.
  Float_t Left_0_17_dE[N_STRIPS];
  Float_t RightdE[N_STRIPS];
  Float_t Cathode;

  // MC truth branches (live on MCTree).
  // Energy sentinels: -1.0 = stopped in the gas, -2.0 = N/A for this event
  // (e.g. beam_energy_exit on a reacted event, evap slots of a disallowed
  // step).
  // *_stop_strip: 0..17 = stopped in that readout strip, -1 = did not stop
  // in a readout strip (exited the gas or stopped in a dead layer —
  // disambiguate via *_stop_z and the exit-energy branch), -2 = N/A.
  // *_stop_x/y/z hold the stop point, or the point where the particle left
  // the active volume if it exited.
  Int_t n_steps;        // configured reaction steps; on-disk array length
  Int_t reaction_strip; // strip selected for the vertex (-1 = unreacted run)
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

  std::ofstream Log;
  std::ofstream EnergeticsLog;

  TSystem *gSystem;
  Float_t MaxMemory;

  // Worker id (0 = master / single-threaded). Used to name per-worker
  // TGeoManagers and per-worker output files.
  Int_t workerId_ = 0;
  // Master (workerId 0) is chatty; MT workers stay quiet so only the
  // runMultiThreaded summary lines reach stdout.
  Bool_t verbose_ = true;
  TString ctrlFilePath_;

  static constexpr Double_t c = 29.9792458; // speed of light, cm/ns
  static constexpr Double_t pi = 3.14159265359;

  struct controlFileParams {
    TString gas = "4He";
    Float_t pressure = 760.0;    // Torr
    Float_t temperature = 293.0; // K
    Int_t ELossBins;
    Float_t MaxELoss;
    TString beamName;
    Float_t dEdxScaleBeam = 1.0;
    TString target;
    TString compound;
    Int_t NumEvapPart;
    // How much of the available energy the residue is left holding as internal
    // excitation. "forced" is the historical behaviour: Ex uniform on
    // [2/3, 1] x EneAvail, which favours energetically-allowed evaporation
    // chains but is wrong for a single-step reaction -- it locks away ~83% of
    // the available energy, so the product is too slow, stops too early, and
    // deposits too sharply just after the vertex.
    //   0 = forced (default, unchanged)  1 = ground state  2 = uniform
    // Stopping-power model for the gas: 0 = catima (built in), 1 = SRIM tables
    // read from disk. SRIM tables are per (ion, gas, pressure, temperature) and
    // are generated once by the srim-cache tool; a missing table is an error
    // rather than a silent fallback, since the two models disagree by ~10% and
    // a quiet substitution would invalidate any dedx_scale calibrated on one.
    Int_t stoppingModel = 0;

    Int_t residueExc = 0;
    // CM angular distribution of the two-body exit channel.
    //   0 = isotropic (default): cos(theta_CM) uniform on [-1, 1]. Correct for
    //       a compound-nucleus channel, where the residue's energy barely
    //       depends on angle because the ejectile is light.
    //   1 = Rutherford: dsigma/dOmega ~ 1/sin^4(theta_CM/2), the right shape
    //       for elastic/inelastic scattering off the gas. With an alpha
    //       ejectile the residue's energy swings ~20% across the angular
    //       range, so isotropic sampling makes the scattered beam lose ~19%
    //       of its energy where the real, forward-peaked process loses ~2%.
    Int_t angularDist = 0;
    // Small-angle cutoff for the Rutherford draw, in degrees. The cross
    // section diverges as theta -> 0; events below this transfer too little
    // energy to be visible anyway. Only used when angularDist == 1.
    Double_t thetaCmMinDeg = 1.0;
    static const Int_t MaxNumEvapPart = 10;
    TString *res = new TString[MaxNumEvapPart];
    Float_t *dEdxScaleRes = new Float_t[MaxNumEvapPart];
    Int_t *colorRes = new Int_t[MaxNumEvapPart];
    TString *evap = new TString[MaxNumEvapPart];
    Float_t *dEdxScaleEvap = new Float_t[MaxNumEvapPart];
    Int_t *colorEvap = new Int_t[MaxNumEvapPart];
    Double_t BeamEnergy; // MeV — KE at the accelerator (before windows)
    Double_t KbFWHM = 0; // MeV FWHM at the accelerator
    TString entranceMaterial = "Ti";
    TString exitMaterial = "Ti";
    // Layer thickness, stored as either mg/cm² (areal density) or μm (linear
    // length along beam axis), keyed off the matching *ByLength flag.
    // Each layer's dedx_scale multiplies its stopping power; it is applied as
    // an equivalent thickness, since the energy lost in a layer is the
    // integral of dE/dx over its length (see BuildWindows).
    Double_t entranceThickness = 0.9;
    Bool_t entranceByLength = false;
    Double_t entranceScale = 1.0;
    Double_t exitThickness = 0.9;
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
    // Per-electrode anode noise as relative resolution in % FWHM of the
    // electrode's energy deposit. 34 entries indexed by
    // ElectrodeIndex(stpid, col). -1 means "no noise on this electrode".
    // Scalar TOML eres = X broadcasts X to all 34 entries (anodes only).
    std::array<Double_t, kNumElectrodes> Eres;
    // Independent cathode-readout noise (% FWHM of the summed cathode
    // energy). Set via the `Cathode` key inside the [detector.eres] table
    // (scalar broadcast does not touch it).
    Double_t EresCathode = -1;
    Int_t NEvents;
    Int_t Wait;       // 1: canvas waits for user click; 0: no wait
    Int_t Update;     // 1: update visuals per event; 0: don't
    Double_t MaxTime; // ns
    Double_t SimStep; // cm
    Int_t Method;     // 0: Simulate; 1: GenerateTraceDatabase
    Int_t Threads = 1;
    TString FileName;
    TString FileOpt;
    Int_t reacClass;
    Int_t PrintOpt = 0;
  };

  controlFileParams ctf;
};

#endif

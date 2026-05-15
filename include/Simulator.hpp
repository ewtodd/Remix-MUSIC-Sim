#ifndef SIMULATOR_HPP
#define SIMULATOR_HPP

#include <algorithm>
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
  void SetStripEnergyResolution(Float_t Sigma /*MeV*/);
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
  Double_t EneSigma;
  Double_t *minEx;

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
  Float_t LeftdE[N_STRIPS];
  Float_t RightdE[N_STRIPS];
  Float_t TotaldE[N_STRIPS];
  Float_t Cathode;

  // MC truth branches (live on MCTree). Exit-energy sentinels:
  //   -1.0 = stopped in gas, -2.0 = N/A for this event (e.g. Kbeam_exit on
  //   a reacted event).
  Int_t reacStp;
  Float_t BeamEnergyAccel; // KE at the accelerator (= ctf.BeamEnergy)
  Float_t Kbi;             // KE at the gas surface (after entrance window)
  Float_t Kbr;             // KE at the reaction point
  Float_t Kbeam_exit;
  Float_t DeadUS_dE; // total dE in upstream dead gas layer (not read out)
  Float_t DeadDS_dE; // total dE in downstream dead gas layer
  Float_t *Kl;
  Float_t *Kh;
  Float_t *Kl_exit;
  Float_t *Kh_exit;
  Float_t *theta_CM;
  Float_t *phi_CM;
  Float_t *theta_l;
  Float_t *phi_l;
  Float_t *theta_h;
  Float_t *phi_h;
  Float_t *xfl;
  Float_t *yfl;
  Float_t *zfl;
  Float_t xr;
  Float_t yr;
  Float_t zr;
  Float_t xfe;
  Float_t yfe;
  Float_t zfe;
  Int_t resID;

  std::ofstream Log;

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
    Double_t entranceThickness = 0.9;
    Bool_t entranceByLength = false;
    Double_t exitThickness = 0.9;
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
    Double_t Eres = -1; // MeV — ad-hoc post-physics noise σ per strip
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

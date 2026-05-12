//===============================================================================
//--- MUSIC simulator ----------------------------------------------------------|
// Written by Daniel Santiago-Gonzalez                                          |
// To get the latest version type:                                              |
//  git clone https://gitlab.phy.anl.gov/music/sim.git                          |
// Compile with make                                                            |
//===============================================================================

#ifndef MUSIC_Simulator_hpp_INCLUDED
#define MUSIC_Simulator_hpp_INCLUDED

//C and C++ libraries
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <future>
#include <vector>
#include <algorithm>
// #include <sys/types.h>
// #include <sys/sysinfo.h>

//ROOT libraries
#include <TCanvas.h>
#include <TDirectory.h>
#include <TFileMerger.h>
#include <TROOT.h>
#include <TEveArrow.h>
#include <TEveGeoNode.h>
#include <TEveLine.h>
#include <TEveManager.h>
#include <TEveStraightLineSet.h>
#include <TF1.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>
#include <TGraph.h>
#include <TGLViewer.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TPolyLine.h>
#include <TRandom3.h>
#include <TSlider.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TTree.h>

// Useful libraries
//#include "EnergyLoss.hpp"
//#include "FourVector.hpp"
#include "Particle.hpp"
#include "NuclideFinder.hpp"

// NOTE: The inclusion of the header files in this .hpp instead of the .cpp file is
//       justified by "Headers and Includes: Why and How"  (see section 5 of 
//       http://www.cplusplus.com/forum/articles/10627/).

class MUSIC_Simulator {

public:
  MUSIC_Simulator(int workerId = 0);
  int loadCtrlFile(char* fileName);
  int run();

  // Overrides used by the MT driver to repurpose a per-worker instance.
  void OverrideNEvents(int n)             { ctf.NEvents = n; }
  void OverrideOutputFile(const std::string& f) { ctf.FileName = f; }
  void OverrideThreads(int t)             { ctf.Threads = t; }
  void DisableVisualization()             { ctf.Update = 0; ctf.Wait = 0; }
  void SeedRandom(unsigned long s);

  void CalculateCMEnergyRange();   // <- Do we need this?
  void CalculateExcEnergyRange();  // <- Do we need this?
  void GenerateTraceDatabase(std::string FileName,
			     double ThCMMin, double ThCMMax, int ThSteps,
			     double PhiCMMin, double PhiCMMax, int PhiSteps,
			     double MaxTime, double UserStep, int UpdateEnabled=0, int Wait=0);
  int SetAnode(short Trans/*From 0 to 100*/, int ELossBins=400, float MaxELoss=5);
  void SetBeamParticle(std::string Name, int Color, float dEdxScale=1.0);
  void SetCompoundParticle(std::string Name);
  void SetDecayDaughter1(std::string Name, int Color);
  void SetDecayDaughter2(std::string Name, int Color);
  void SetEvapResAndPart(std::string ResName, int ResColor,
			 std::string ParName, int ParColor,
			 float dEdxScaleRes=1.0, float dEdxScalePar=1.0);
  void SetHeavyParticle(std::string Name, int Color, int NEexc=0,
			double* Eexc=0/*MeV*/);
  void SetLightParticle(std::string Name, int Color);

  void SetPrintLevel(int PrintLevel);
  void SetROOTSystemPointer(TSystem* gSystem);
  void SetStripEnergyResolution(float Sigma/*MeV*/);
  void SetTargetParticle(std::string Name);
  void Simulate(int StpID, int NEvents, double MaxTime, double UserStep, int UpdateVis=0, int Wait=0,
		TFile* ROOTfile=0);
  void Simulate(int StpID, double ThCMMin, double ThCMMax, int ThSteps, double PhiCMMin, 
		double PhiCMMax, int PhiSteps, double MaxTime, double UserStep, int Wait=0);
  void WriteTraces(char* FileName);


private:
  void BuildGasMaterial();
  void BuildWindows();
  void BuildDegrader();
  catima::Material LookupMaterial(const std::string& name);
  catima::Material BuildSolidMaterial(const std::string& name, double thickness_mg_per_cm2);
  catima::Material BuildBulkMaterial(const std::string& name, double thickness_um);
  double EnergyOutOfMaterial(int A, int Z, double Ein_MeV, const catima::Material& mat);
  // Propagate through a material with catima-derived Gaussian straggling.
  // Returns Ein unchanged for neutrals (Z<=0).
  double EnergyThroughWithStraggling(int A, int Z, double Ein_MeV, const catima::Material& mat);
  void LoadHardcodedAnodeGeometry();
  void FinalizeEvent(int eventIndex);
  void ComputeExitEnergies();
  int CheckMemoryUsage(int Print=0);
  void ComputeDetectorResponse(int event, int reacStp, int UpdateVis);
  void CreateTracesAndTrajectories();
  void DrawMUSIC(TEveManager* gEve, short Transparency /*From 0 to 100*/);
  void InitCTF();
  TTree* InitTree(TFile* ROOTfile, std::string FileOpt);
  void PrintCompoundEexc(double Kb, double** DeltaEB);
  void PrintEnergetics(double Kb, double** DeltaEB);
  int PropagateParticle(Particle* PO, int Event, double MaxTime, double UserStep, double ** DE);
  void ResetBranches();
  void SetInitialKinematics(double Kbi);
  int SetReactionKinematics(double Kbr, double zr, double tof, double theta_CM=-1, double phi_CM=-1);
  void UpdateVisuals(int event, double Kbr, double zr, double TOF, int Wait=0);
  
  // Useful random number.
  //  TRandom3* Rdm;
  TRandom* Rdm;
  

  // Particle related stuff.
  Particle* Beam;
  Particle* Target;
  Particle* Compound;
  Particle* Heavy;
  Particle* DeDau1;
  Particle* DeDau2;
  Particle* Light;
  Particle** EvaP;
  Particle** EvaR;
  int maxEvaporations;
  int numEvaporations;
  double Kb_after_window;
  double EneSigma;
  double* minEx;

  // Arrays where the energy loss in each strip will be saved.
  double** DeltaEB_ave;// average beam energy loss
  double** DeltaEB;    // beam
  double** DeltaEL;    // light
  double** DeltaEH;    // heavy
  double** DeltaED1;   // decay daughter1
  double** DeltaED2;   // decay daughter2
  double*** DeltaE_EvaP;
  double*** DeltaE_EvaR;

  std::string Name;

  int PrintLevel;
 
  double* SegLength;
  double* SegEexcRange;
  double* SegCMERange;
  double CMEMax;
  double CMEMin;
  double EexcMax;
  double EexcMin;

  TH2F* HCT;
  TH2F* HCTB;
  TH2F* HPT;
  TGraph** Trace;      // Detector trece
  TGraph** TraceUB;    // Unreacted beam trace
  TGraph** TraceB;     // Beam trace
  TGraph*** TraceER;   // Evaporation residue traces
  TGraph*** TraceEP;   // Evaporated particle traces
  bool tracesCreated;

  TH1F* TraceMult;
  int NTraces;

  TF1* Gaussian;   // For randomizing the detector response

  int NEvents;

  // Visual stuff
  TEveArrow* TrackBeam;
  TEveArrow** TrackEvaP;
  TEveArrow** TrackEvaR;
  TEveManager* Eve;
  TCanvas* TraceCan;
  TLegend* LegCol;
  TLegend* LegPart;
  TPaveText* LabelKine;


  // Geometry stuff
  int AnodeRows;
  int AnodeCols;
  double AnodeDepth;  // distance along z-axis
  double AnodeLength; // distance along x-axis
  double AnodeHeight; // distance along y-axis
  short** AnodeColor;
  double** AnodeDX;
  double** AnodeDY;
  double** AnodeDZ;
  std::string** AnodeSegName;
  int** AnodeStpID;
  TGeoManager* Geo;
  TGeoMaterial* MatVacuum;
  TGeoMedium* Vacuum;
  TGeoVolume*** VolAnode;
  TGeoVolume* VolTop;
  TEveGeoTopNode* TopNode;
  
  // Nuclide finder
  NuclideFinder* NuF;

  // Gas material for stopping-power calculations (catima)
  catima::Material gas_;
  // Entrance / exit windows (configurable material + thickness in mg/cm^2)
  catima::Material entranceWindow_;
  catima::Material exitWindow_;
  // Optional upstream degrader (material + length in microns along beam axis).
  catima::Material degrader_;
  bool hasDegrader_ = false;
  // Per-layer enable flags. Disabled = the layer is skipped in all chains
  // (straggling, ComputeExitEnergies, PreWarmCatima). Triggered by
  // Pressure <= 0 for the gas and Thickness <= 0 for the windows.
  bool gasEnabled_ = true;
  bool entranceWindowEnabled_ = true;
  bool exitWindowEnabled_ = true;
  // Beam KE at the gas surface (after entrance window). Derived from ctf.BeamEnergy.
  double Kb_at_gas;

  // TTree stuff. Layout matches the upstream EventBuilderNearestGrid output:
  // an "event" tree carrying detector-level branches plus a sibling "MC" tree
  // carrying truth-only branches (friended; rows correspond 1:1).
  TTree* SimTree;   // stored as "event"
  TTree* MCTree;    // stored as "MC"
  static const int N_STRIPS = 18;
  static const int N_CHAN   = 36;
  Float_t LeftdE[N_STRIPS];
  Float_t RightdE[N_STRIPS];
  Float_t TotaldE[N_STRIPS];
  ULong64_t AllTimestamps[N_CHAN];
  UInt_t AllFlags[N_CHAN];
  Int_t Hits[N_CHAN];
  Float_t Cathode;
  Float_t Grid;
  Bool_t IsComplete;
  // MC truth branches (live on MCTree).
  // Sentinel values for exit-energy branches: -1.0 = particle stopped in gas,
  // -2.0 = particle does not exist for this event (e.g. Kbeam_exit on a reacted event).
  int reacStp;
  float BeamEnergyAccel;  // KE at the accelerator (= ctf.BeamEnergy)
  float Kbi;              // KE at the gas surface (= after entrance window)
  float Kbr;              // KE at the reaction point
  float Kbeam_exit;       // unreacted-beam KE after exit window (sentinel for reaction events)
  float* Kl;
  float* Kh;
  float* Kl_exit;         // light-product KE after exit window
  float* Kh_exit;         // heavy-product KE after exit window
  float* theta_CM;
  float* phi_CM;
  float* theta_l;
  float* phi_l;
  float* theta_h;
  float* phi_h;
  float* xfl;
  float* yfl;
  float* zfl;
  float xr;
  float yr;
  float zr;
  float xfe;
  float yfe;
  float zfe;
  int resID;

  
  // Log file
  std::ofstream Log;

  // To be used in CheckMemory
  TSystem* gSystem;
  float MaxMemory;

  // Worker id (0 for the main thread or single-threaded mode). Used to name
  // per-worker TGeoManagers and per-worker output files.
  int workerId_ = 0;
  // Gate for chatty stdout. True for the master (workerId 0), false for
  // every MT worker. Single-threaded mode runs as master so verbosity is
  // unchanged; in MT mode workers stay quiet and only the runMultiThreaded
  // summary lines reach stdout.
  bool verbose_ = true;
  // Ctrl-file path remembered so the MT driver can re-load it per worker.
  std::string ctrlFilePath_;

  int runMultiThreaded();
  void PreWarmCatima();

  // Use these two lines when compiling with root-config ver5
  //  static const double c = 29.9792458;  // Speed of light in cm/ns.
  //  static const double pi = 3.14159265359;
  // Use these two lines when compiling with root-config ver6
  static constexpr double c = 29.9792458;  // Speed of light in cm/ns.
  static constexpr double pi = 3.14159265359;



  struct controlFileParams {
    std::string gas = "4He";    // gas species name; see BuildGasMaterial for supported values
    float pressure = 760.0;     // Torr
    float temperature = 293.0;  // K
    int ELossBins;
    float MaxELoss;
    // beam
    std::string beamName;
    float dEdxScaleBeam=1.0;
    std::string target;
    std::string compound;
    int NumEvapPart;
    static const int MaxNumEvapPart=10;
    std::string* res = new std::string[MaxNumEvapPart];
    float* dEdxScaleRes = new float[MaxNumEvapPart];
    int* colorRes = new int[MaxNumEvapPart];
    std::string* evap = new std::string[MaxNumEvapPart];
    float* dEdxScaleEvap = new float[MaxNumEvapPart];
    int* colorEvap = new int[MaxNumEvapPart];
    double BeamEnergy;   // MeV - Kinetic energy of the beam at the accelerator (before windows)
    double KbFWHM=0;     // MeV - Beam energy spread (full-width half maximum), at the accelerator
    std::string entranceMaterial = "Ti";
    std::string exitMaterial     = "Ti";
    double entranceThickness = 0.9;  // mg/cm^2 - upstream window areal density
    double exitThickness     = 0.9;  // mg/cm^2 - downstream window areal density
    // Optional bulk degrader between accelerator and entrance window.
    // Empty material disables it; length is along the beam axis in microns.
    std::string degraderMaterial = "";
    double degraderLength = 0.0;     // microns
    // Reaction strip selection. Either set 'strip' to a single value (use
    // -1 for unreacted beam, 0..17 for a reaction strip), OR set both
    // 'stripFirst' and 'stripLast' to span a range. The sentinel below
    // marks "unset" so a missing key is caught explicitly rather than
    // silently falling back to a default.
    static constexpr int kStripUnset = -99999;
    int strip      = kStripUnset;
    int stripFirst = kStripUnset;
    int stripLast  = kStripUnset;
    double Eres=-1;  // MeV - Ad-hoc post-physics electronic noise sigma per strip.
                     // Negative = off (default). Gas straggling is now sampled from
                     // catima per integration step; Eres is reserved for tuning to data.
    int NEvents;     // Number of simulated events (recommendation: keep it <1000)
    int Wait;        // 1 - canvas waits for user's double click, 0 - no wait
    int Update;      // 1 - update visuals for every event, 0 - don't
    double MaxTime;  // ns - max time for an event
    double SimStep;  // cm - simulation steps size
    int Method;      // Select the simulation method: 0 - Simulate, 1 - GenerateTraceDatabase
    int Threads = 1; // Number of worker threads for the event loop (default: single-threaded)
    std::string FileName;
    std::string FileOpt;
    int reacClass;
    int PrintOpt=0;
  };
  
  controlFileParams ctf;

};

#endif

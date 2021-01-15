//===============================================================================
//--- MUSIC simulator ----------------------------------------------------------|
// Written by Daniel Santiago-Gonzalez                                          |
// ver 2.0 (2016/4)                                                             |
// To get the latest version type:                                              |
//  git clone https://gitlab.phy.anl.gov/music/sim.git                          |
//                                                                              |
// Compile by typing: root MakeAll.C                                            |
//                                                                              |
// May need libroot-math-physics-dev. Install using synaptic package manager.   |
// http://packages.ubuntu.com/trusty/i386/libroot-math-physics-dev/filelist     |
//===============================================================================

#ifndef MUSIC_Simulator_hpp_INCLUDED
#define MUSIC_Simulator_hpp_INCLUDED

//C and C++ libraries
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
// #include <sys/types.h>
// #include <sys/sysinfo.h>

//ROOT libraries
#include <TCanvas.h>
#include <TDirectory.h>
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
#include "EnergyLoss.hpp"
#include "FourVector.hpp"
#include "Particle.hpp"
#include "NuclideFinder.hpp"

// NOTE: The inclusion of the header files in this .hpp instead of the .cpp file is
//       justified by "Headers and Includes: Why and How"  (see section 5 of 
//       http://www.cplusplus.com/forum/articles/10627/).

class MUSIC_Simulator {

public:
  MUSIC_Simulator();
  int loadCtrlFile(char* fileName);
  int run();

  void CalculateCMEnergyRange();   // <- Do we need this?
  void CalculateExcEnergyRange();  // <- Do we need this?
  void GenerateTraceDatabase(std::string FileName,
			     double ThCMMin, double ThCMMax, int ThSteps,
			     double PhiCMMin, double PhiCMMax, int PhiSteps,
			     double MaxTime, double UserStep, int UpdateEnabled=0, int Wait=0);
  void SetAnode(std::string AnodeGeomFile, short Trans/*From 0 to 100*/, int ELossBins=400, 
		float MaxELoss=5);
  void SetBeamParticle(std::string Name, int Color, std::string ELossFile, double KineticE/*MeV*/);
  //  void SetBeamSpot(double diameter);
  void SetCompoundParticle(std::string Name);
  void SetDecayDaughter1(std::string Name, int Color, std::string ELossFile);
  void SetDecayDaughter2(std::string Name, int Color, std::string ELossFile);
  void SetEvapResAndPart(std::string ResName, std::string ResELossFile, int ResColor,
			 std::string ParName, std::string ParELossFile, int ParColor);
  void SetHeavyParticle(std::string Name, int Color, std::string ELossFile, int NEexc=0,
			double* Eexc=0/*MeV*/);
  void SetLightParticle(std::string Name, int Color, std::string ELossFile);

  void SetPrintLevel(int PrintLevel);
  void SetROOTSystemPointer(TSystem* gSystem);
  void SetStripEnergyResolution(float Sigma/*MeV*/);
  void SetTargetParticle(std::string Name);
  void Simulate(int StpID, int NEvents, double MaxTime, double UserStep, int UpdateVis=0, int Wait=0,
		std::string FileName="", std::string FileOpt="recreate");
  void Simulate(int StpID, double ThCMMin, double ThCMMax, int ThSteps, double PhiCMMin, 
		double PhiCMMax, int PhiSteps, double MaxTime, double UserStep, int Wait=0);
  void WriteTraces(char* FileName);


private:
  int CheckMemoryUsage(int Print=0);
  void ComputeDetectorResponse(int event, int reacStp, int UpdateVis);
  void CreateTracesAndTrajectories(int NEvents);
  void DrawMUSIC(TEveManager* gEve, short Transparency /*From 0 to 100*/);
  TTree* InitTree(TFile* ROOTfile, std::string FileOpt);
  void PrintCompoundEexc(double Kb, double** DeltaEB);
  void PrintEnergetics(double Kb, double** DeltaEB);
  int PropagateParticle(Particle* PO, int Event, double MaxTime, double UserStep, double ** DE);
  void SetInitialKinematics(double Kbi);
  int SetReactionKinematics(double Kbr, double zr, double tof, double theta_CM=-1, double phi_CM=-1);
  void UpdateVisuals(int event, double Kbr, double zr, double TOF, int Wait=0);
  
  // Useful random number.
  TRandom3* Rdm;

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
  int MaxEva;
  int CurEva;
  double Kb_after_window;
  double EneSigma;

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

  TH1F* TraceMult;
  int NTraces;

  TF1* Gaussian;   // For randomizing the detector response

  int NEvents;

  // Visual stuff
  TEveStraightLineSet** TrajH;
  TEveStraightLineSet** TrajD1;
  TEveStraightLineSet** TrajD2;
  TEveStraightLineSet** TrajL;
  TEveStraightLineSet*** TrajEvaP;
  TEveStraightLineSet*** TrajEvaR;
  TEveArrow* TrackBeam;
  TEveArrow** TrackEvaP;
  TEveArrow** TrackEvaR;
  TEveManager* Eve;
  TCanvas* TraceCan;
  TLegend* LegCol;
  TLegend* LegPart;
  TPaveText* LabelKine;


  // Geometry stuff
  int AnodeStps;
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

  // TTree stuff similar to the one used for experimental data.
  TTree* SimTree;
  float strip0;
  static const int ExpAnodeStps = 16;
  float* de_l;
  float* de_r;
  int* seg;
  float strip17;
  float cathode;
  int reacStp;
  float Kb;
  float* Kl;
  float* Kh;
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


  
  // Log file
  std::ofstream Log;

  // To be used in CheckMemory
  TSystem* gSystem;
  float MaxMemory;

  // Use these two lines when compiling with root-config ver5
  //  static const double c = 29.9792458;  // Speed of light in cm/ns.
  //  static const double pi = 3.14159265359;
  // Use these two lines when compiling with root-config ver6
  static constexpr double c = 29.9792458;  // Speed of light in cm/ns.
  static constexpr double pi = 3.14159265359;



  struct controlFileParams {
    std::string SRIMdir;
    int pressure; // Torr
    std::string AnodeGeom;
    int ELossBins = 300;
    float MaxELoss = 25;
    // beam
    std::string beam;
    std::string SRIMbeam;
    std::string target;
    std::string compound;
    static const int NumEvapPart = 10;
    std::string* res = new std::string[NumEvapPart];
    std::string* SRIMres = new std::string[NumEvapPart];
    std::string* evap = new std::string[NumEvapPart];
    std::string* SRIMevap = new std::string[NumEvapPart];
    int* color = new int[NumEvapPart];
    double Kb;  // MeV - Energy of the beam after the Ti window and degrader (if any)
    int strip;      // Strip where reaction takes place
    float Eres; // MeV - Strip energy resolution (larger values increase signal randomness)
    int NEvents;  // Number of simulated events (recommendation: keep it <1000)
    int Wait;       // 1 - canvas waits for user's double click, 0 - no wait
    int Update;     // 1 - update visuals for every event, 0 - don't
    double MaxTime;   // ns - max time for an event
    double SimStep;  // cm - simulation steps size
    int Method;     // Select the simulation method: 0 - Simulate, 1 - GenerateTraceDatabase
    std::string FileName;
    std::string FileOpt;
  };
  
  controlFileParams ctf;

};

#endif

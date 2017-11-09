//===============================================================================
//--- MUSIC simulator ----------------------------------------------------------|
// Written by Daniel Santiago-Gonzalez                                          |
// ver 2.0 (2016/4)                                                             |
// To get the latest version type:                                              |
// git clone https://dasago@bitbucket.org/music_anl/music_simulator.git         |
//                                                                              |
// Compile by typing: root -l MakeAll.C                                         |
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

//ROOT libraries
#include <TCanvas.h>
#include <TEveArrow.h>
#include <TEveGeoNode.h>
#include <TEveManager.h>
#include <TF1.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>
#include <TGraph.h>
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
#include <TTree.h>

// Useful libraries
#include "../physicstools/EnergyLoss.hpp"
#include "../physicstools/FourVector.hpp"
#include "../physicstools/Particle.hpp"
#include "../physicstools/NuclideFinder.hpp"

// NOTE: The inclusion of the header files in this .hpp instead of the .cpp file is
//       justified by "Headers and Includes: Why and How"  (see section 5 of 
//       http://www.cplusplus.com/forum/articles/10627/).

class MUSIC_Simulator {
public:
  MUSIC_Simulator();
  void CalculateCMEnergyRange();   // <- Do we need this?
  void CalculateExcEnergyRange();  // <- Do we need this?
  void GenerateTraceDatabase(std::string FileName,
			     double ThCMMin, double ThCMMax, int ThSteps,
			     double PhiCMMin, double PhiCMMax, int PhiSteps,
			     double MaxTime, double UserDT, int Wait=0);
  void SetAnode(std::string AnodeGeomFile, short Trans/*From 0 to 100*/);
  void SetBeamParticle(std::string Name, int Color, std::string ELossFile, double KineticE/*MeV*/);
  //  void SetBeamSpot(double diameter);
  void SetCompoundParticle(std::string Name);
  void SetHeavyParticle(std::string Name, int Color, std::string ELossFile, int NEexc=0,
			double* Eexc=0/*MeV*/);
  void SetLightParticle(std::string Name, int Color, std::string ELossFile);
  void SetPrintLevel(int PrintLevel);
  void SetStripEnergyResolution(float Sigma/*MeV*/);
  void SetTargetParticle(std::string Name);
  void Simulate(int StpID, int NEvents, double MaxTime, double UserDT, int Wait=0);
  void Simulate(int StpID, double ThCMMin, double ThCMMax, int ThSteps, double PhiCMMin, 
		double PhiCMMax, int PhiSteps, double MaxTime, double UserDT, int Wait=0);
  void WriteTraces(char* FileName);

private:
  void ComputeDetectorResponse(int event);
  void CreateTracesAndTrajectories(int NEvents);
  void DrawMUSIC(TEveManager* gEve, short Transparency /*From 0 to 100*/);
  void PrintCompoundEexc(double Kb, double** DeltaEB);
  double** PropagateParticle(Particle* PO, int Event, double MaxTime, double UserDT);
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
  Particle* Light;
  double Kb_after_window;
  double EneSigma;

  // Arrays where the energy loss in each strip will be saved.
  double** DeltaEB_ave;// average beam energy loss
  double** DeltaEB;    // beam
  double** DeltaEL;    // light
  double** DeltaEH;    // heavy

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
  TGraph*** Trace;
  TGraph*** TraceH;
  TGraph*** TraceL;
  TGraph** TraceB;
  int NTraces;

  TF1* Gaussian;   // For randomizing the detector response

  int NEvents;

  // Visual stuff
  TEveStraightLineSet** TrajH;
  TEveStraightLineSet** TrajL;
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
  float* andl;
  float* andr;
  int* seg;
  float strip17;
  float cathode;
  float Kl;
  float Kh;
  float theta_l;
  float theta_h;
  float phi_l;
  float phi_h;

  // Log file
  ofstream Log;

  static const double c = 29.9792458;  // Speed of light in cm/ns.
  static const double pi = 3.14159265359;
};

#endif

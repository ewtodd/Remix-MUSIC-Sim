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
#include "../../PhysicsTools/EnergyLoss.hpp"
#include "../../PhysicsTools/FourVector.hpp"
#include "../../PhysicsTools/Particle.hpp"
#include "../../NuclideFinder/NuclideFinder.hpp"

// NOTE: The inclusion of the header files in this .hpp instead of the .cpp file is
//       justified by "Headers and Includes: Why and How"  (see section 5 of 
//       http://www.cplusplus.com/forum/articles/10627/).

class MUSIC_Simulator {
public:
  MUSIC_Simulator();
  void CalculateCMEnergyRange();   // <- Do we need this?
  void CalculateExcEnergyRange();  // <- Do we need this?
  void DrawMUSIC(TEveManager* gEve, short Transparency /*From 0 to 100*/);
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
  void Simulate(int StpNum, int NEvents, double MaxTime, double UserDT, int Wait=0);
  void WriteTraces(char* FileName);

private:
  void SetInitialKinematics(double Kbi);
  int SetReactionKinematics(double Kbr, double zr, double tof);
  double** PropagateParticle(Particle* PO, int Event, double MaxTime, double UserDT);
  void PrintCompoundEexc(double Kb, double** DeltaEB);

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

  int NEvents;
  TEveStraightLineSet** TrajH;
  TEveStraightLineSet** TrajL;
  TEveManager* Eve;
  TEveGeoTopNode* TopNode;

  // Geometry stuff.
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
  TGeoManager* Geo;
  TGeoMaterial* MatVacuum;
  TGeoMedium* Vacuum;
  TGeoVolume*** VolAnode;
  TGeoVolume* VolTop;
 
  // Nuclide finder
  NuclideFinder* NuF;

  static const double c = 29.9792458;  // Speed of light in cm/ns.
  static const double pi = 3.14159265359;
};

#endif

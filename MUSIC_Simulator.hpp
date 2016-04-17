// May need libroot-math-physics-dev. Install using synaptic package manager.
// http://packages.ubuntu.com/trusty/i386/libroot-math-physics-dev/filelist
#ifndef MUSIC_Simulator_hpp_INCLUDED   
#define MUSIC_Simulator_hpp_INCLUDED   

//C and C++ libraries
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

//ROOT libraries
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2.h>
#include <TLine.h>
#include <TMath.h>
#include <TPolyLine.h>
#include <TRandom3.h>
#include <TEveArrow.h>
#include <TEveGeoNode.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>
#include <TGeoXtru.h>
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
  double* CalculateELoss(Particle* P, int Event);
  void DrawMUSIC(TEveManager* gEve, short Transparency /*From 0 to 100*/);
  void DrawTrajecotries(TEveManager* gEve);
  void SetAnode(std::string AnodeGeomFile, short Trans/*From 0 to 100*/);
  void SetBeamParticle(std::string Name, int Color, std::string ELossFile, double KineticE/*MeV*/);
  //  void SetBeamSpot(double diameter);
  void SetCompoundParticle(std::string Name);
  void SetHeavyParticle(std::string Name, int Color, std::string ELossFile, int NEexc=0, 
			double* Eexc=0/*MeV*/);
  void SetLightParticle(std::string Name, int Color, std::string ELossFile); 
  void SetStripEnergyResolution(float Sigma/*MeV*/);
  void SetTargetParticle(std::string Name);
  void Simulate(int StpNum, int NEvents, double MaxTime, double UserDT);
  void WriteTraces(char* FileName);

private:
  void SetInitialKinematics(double Kbi);
  void SetReactionKinematics(double Kbr, double zr);
  double* CalcAverageBeamELoss();
  double** PropagateParticle(Particle* PO, int Event, double MaxTime, double UserDT);

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
 
  double* SegLength;
  double* SegEexcRange;
  double* SegCMERange;
  double CMEMax;
  double CMEMin;
  double EexcMax;
  double EexcMin;

  TH2F* HELoss;
  TH2F* HELossB;
  TH2F** HELossC;
  TGraph*** Trace;
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
  TGeoMaterial* MatAl;
  TGeoMaterial* MatSi;
  TGeoMaterial* MatVacuum;
  TGeoMedium* Al;
  TGeoMedium* CD2;
  TGeoMedium* CF4;
  TGeoMedium* D2;
  TGeoMedium* He3;
  TGeoMedium* He4;
  TGeoMedium* Kapton;
  TGeoMedium* LiF;
  TGeoMedium* Mylar;
  TGeoMedium* Si;
  TGeoMedium* Ti;
  TGeoMedium* Vacuum;
  TGeoVolume* VolIC;
  TGeoVolume* VolICBkFlange;
  TGeoVolume* VolICFlange;
  TGeoVolume* VolICPSGFrame;
  TGeoVolume** VolICSec;
  TGeoVolume* VolICWin;
  TGeoVolume* VolSD;
  TGeoVolume* VolSolDSDoor;
  TGeoVolume*** VolAnode;
  TGeoVolume* VolSolUSDoor;
  TGeoVolume* VolTgt;
  TGeoVolume* VolTgtFrame;
  TGeoVolume* VolTgtWinDS;
  TGeoVolume* VolTgtWinUS;
  TGeoVolume* VolTop;
 
  // Nuclide finder
  NuclideFinder* NuF;


  static const double c = 29.9792458;  // Speed of light in cm/ns.
  static const double pi = 3.14159265359;

};

#endif

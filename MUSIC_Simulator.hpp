// May need libroot-math-physics-dev. Install using synaptic package manager.
// http://packages.ubuntu.com/trusty/i386/libroot-math-physics-dev/filelist
#ifndef MUSIC_Simulator_hpp_INCLUDED   
#define MUSIC_Simulator_hpp_INCLUDED   

//C and C++ libraries
#include <iostream>
#include <fstream>
#include <string>

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
#include "../../PhysicsTools/SRIM_Table_Maker.hpp"
#include "../../NuclideFinder/NuclideFinder.hpp"


// NOTE: The inclusion of the header files in this .hpp instad of the .cpp file is
//       justified by "Headers and Includes: Why and How"  (see section 5 of 
//       http://www.cplusplus.com/forum/articles/10627/).


class MUSIC_Simulator {
public:
  MUSIC_Simulator();
  void CalculateCMEnergyRange();
  void CalculateExcEnergyRange();
  double* CalculateELoss(Particle* P, EnergyLoss* PInTgt);
  void CreateMUSIC();
  void DrawMUSIC(TEveManager* gEve, short Transparency /*From 0 to 100*/);
  void DrawTrajecotries(TEveManager* gEve);
  void SetBeamParticle(string ParticleName, double M, int Q, int Color, double KineticE);
  bool SetEnergyLossFile(string ParticleName, string TgtELossFile);
  void SetFusedParticle(string ParticleName, double M, int Q, int NEexc=0, double* Eexc=0);
  void SetHeavyParticle(string ParticleName, double M, int Q, int Color, int NEexc=0, double* Eexc=0);
  void SetLightParticle(string ParticleName, double M, int Q, int Color); 
  void SetParamDirectory(string Dir);
  void SetPreviousCrossSection(string FileName, string Format, short Marker, short Color);
  void SetSegmentLength(int NSegments, float* SegLength /*cm*/);
  void SetStripEnergyResolution(float Sigma);
  void SetTargetParticle(string ParticleName, double M, int Q);
  void ShowCMEnergyRange();
  void ShowPreviousCrossSection(float XMin, float XMax, float YMin, float YMax);
  void Simulate(int NEvents);
  void Simulate(int Reaction, int StpNum, int NEvents);
  void WriteTraces(char* FileName);

private:
  // Useful random number.
  TRandom3* Rdm;

  // Particle related stuff.
  Particle* Beam;
  Particle* Target;
  Particle** Fused;
  Particle** Heavy;
  Particle** Light;
  EnergyLoss* BeamInTgt;
  EnergyLoss** FusedInTgt;
  EnergyLoss** HeavyInTgt;
  EnergyLoss** LightInTgt;

  int MaxParticles;
  int MaxPrevCS;
  int NFusedParticles;
  int NHeavyParticles;
  int NLightParticles;
  int NPrevCS;
  int NSegments;

  double Kb_after_window;
  double EneSigma;

  string ParamDirectory;
  string Name;
 
  double* SegLength;
  double* SegEexcRange;
  double* SegCMERange;
  double CMEMax;
  double CMEMin;
  double EexcMax;
  double EexcMin;
  double VolW, VolH, VolL;

  TGraph** PrevCS;
  TH2F* HELoss;
  TGraph** Trace;
  int NTraces;


  // Geometry stuff.
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
  TGeoVolume* VolRDBody;// New recoil det
  TGeoVolume* VolRDFace;// New recoil det
  TGeoVolume* VolRDSi1;   // New recoil det
  TGeoVolume* VolRDSi2;   // New recoil det
  TGeoVolume* VolRDSi3;   // New recoil det
  TGeoVolume* VolSD;
  TGeoVolume* VolSolDSDoor;
  TGeoVolume*** VolAnode;
  TGeoVolume* VolSolUSDoor;
  TGeoVolume* VolTgt;
  TGeoVolume* VolTgtFrame;
  TGeoVolume* VolTgtWinDS;
  TGeoVolume* VolTgtWinUS;
  TGeoVolume* VolTop;

  // For energy loss tables.
  SRIM_Table_Maker* SRIM;
 
  // Nuclide finder
  NuclideFinder* NuF;


  static const double c = 29.9792458;  // Speed of light in cm/ns.

};

#endif

/*******************************************************************
Header file: EnergyLoss.hpp

Description: Simple class that calculates the energy loss of an ion
  in a gas target. The energy loss table can be loaded from an
  output SRIM file, with the energy of the ion, the electrical, 
  nuclear stopping powers (dE/dx), etc.  The units are assumed to be
  MeV and MeV/mm for the ion's energy and the stopping powers,
  respectively. 

Typical usage:
  EnergyLoss* Ne20InHe4 = new EnergyLoss();
  Ne20InHe4->LoadSRIMFile("N20_in_He4_500Torr_90K.srim");
  //                                         MeV   cm    cm
  double Efinal = Ne20InHe4->GetFinalEnergy(100.0, 1.5, 0.01);

Compile with: 
  g++ -shared -fPIC EnergyLoss.cpp `root-config --cflags --glibs` -o EnergyLoss.so

Author: Daniel Santiago-Gonzalez
2012-Sep
*******************************************************************/
#ifndef EnergyLoss_hpp_INCLUDED   
#define EnergyLoss_hpp_INCLUDED   

#include <iostream>
#include <cmath>
#include <fstream>
#include <string.h>
#include <TGraph.h>

class EnergyLoss{
public:
  EnergyLoss(std::string SRIM_file="", float IonMass=0/*MeV/c^2*/, float dEdxScale=1.0);
  void GetBraggCurve(float InitE, int NSteps, float* Dist, float StepSize);
  void GetBraggCurves(int NCurves, float* InitE, int NSteps, float FinalDist);
  void GetBraggCurves(int NCurves, float* InitE, int NSteps, float* Dist, float StepSize);
  void GetBraggCurves(int NCurves, float InitE, int* NSteps, float** Dist, float StepSize);
  void GetEnergyCurve(float InitE, int NSteps, float* Dist, float StepSize);
  void GetEnergyCurves(int NCurves, float* InitE, int NSteps, float* DistArray, float StepSize);
  void GetEvDCurve(float InitEne, float FinalDepth, int steps);
  double GetFinalEnergy(double InitialEnergy/*MeV*/, double PathLength/*cm*/, double StepSize/*cm*/);
  float GetFinalEnergy(float InitialEnergy/*MeV*/, float PathLength/*cm*/, float StepSize/*cm*/);
  double GetInitialEnergy(float FinalEnergy/*MeV*/, float PathLength/*cm*/, float StepSize/*cm*/);
  double GetOptimumStepSize(float Energy/*MeV*/);
  double GetPathLength(float InitialEnergy/*MeV*/, float FinalEnergy/*MeV*/, float DeltaT/*ns*/);
  double GetTimeOfFlight();
  double GetTimeOfFlight(float InitialEnergy, float PathLength, float StepSize);
  bool LoadSRIMFile(std::string FileName);
  bool LoadLISEFile(std::string FileName);
  void SetIonMass(float IonMass/*MeV/c^2*/);


  bool GoodELossFile;
  TGraph* EvD;
  TGraph** BraggCurve;
  TGraph** EnergyCurve;
  std::string FileName;

private:

  // These methods on their own sometimes do not converge. It is
  // better to use the GetFinalEnergy and GetInitialEnergy methods.
  double GetEnergyLoss(double initial_energy, double distance);
  float GetEnergyLoss(float initial_energy, float distance);

  // Private members
  bool Energy_in_range;
  double* dEdx_e;
  double* dEdx_n;
  double* IonEnergy;
  double IonMass;
  int last_point;
  int points;
  double TOF;
  float dEdxScale;

};

#endif

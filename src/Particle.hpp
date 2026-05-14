/*******************************************************************
Header file: Particle.hpp

Description: class used in MUSIC and HELIOS simulators.

Compile with: 
On linux
  g++ -shared -fPIC Particle.cpp `root-config --cflags --glibs` -o Particle.so
On MacOS (work in progress)
  g++ -shared -fPIC Particle.cpp FourVector.so EnergyLoss.so `root-config --cflags --glibs` -lEve -o Particle.so

Author: Daniel Santiago-Gonzalez
2014
*******************************************************************/

#ifndef Particle_hpp_INCLUDED   
#define Particle_hpp_INCLUDED   

// C++ libraries
#include <iostream>
#include <string.h>

// ROOT libraries
#include <TEveStraightLineSet.h>
#include <TRandom3.h>

// Useful for calculating energy losses in media.
#include "EnergyLoss.hpp"
// My version of a 4-vector class.
#include "FourVector.hpp"


class Particle{
public:
  // Constructors
  Particle(std::string Name, double M=0, int Q=0, bool SaveTrajectory=0);

  // Methods
  void Boost(double BetaX, double BetaY, double BetaZ);
  void Copy(Particle* Other);
  void CopyTrace(int& NumPts, float* t, float* x, float* y, float* z, float* K);
  void GetBeta(double& BetaX, double& BetaY, double& BetaZ);
  short GetColor() {return AttColor;};
  int GetCurrentExcState();
  double GetEexc();
  double GetEexc(int ExcState);
  double GetEnergyLoss(int MediumID, double InitE/*MeV*/, double PathLength/*cm*/);
  double GetFinalEnergy(int MediumID, double InitE, double PathLength);
  double GetFinalEnergyStraggled(int MediumID, double InitE, double PathLength, TRandom* rng);
  double GetInitialEnergy(int MediumID, double FinalE, double PathLength);
  double GetKE();
  double GetOptimumStepSize(int MediumID, double Energy);
  FourVector GetP();
  void GetP(double& P0, double& P1, double& P2, double& P3);
  double GetPathLength(int MediumID, double InitE, double FinalE, double DeltaT);
  double GetPhi();
  double GetPhiX();
  double GetTheta();
  double GetThetaX();
  double GetTimeOfFlight(int MediumID);
  double GetTimeOfFlight(int MediumID, float InitialEnergy, float PathLength);
  void GetTrajectoryAtt(short& Color, short& Style, short& Width);
  void GetX(double& X0, double& X1, double& X2, double& X3);
  void Print(std::ostream& log=std::cout);
  void ResetTrace();
  void ResetKinematics();
  void SetCurrentExcState(int ExcState);
  void SetExcEnergies(int N, double* Eexc, double* Prob=0);
  void SetExcEnergy(double Ex);
  void SetMedium(const catima::Material* gas, float dEdxScale=1.0);
  void SetP(FourVector P);
  void SetP(double P0, double P1, double P2, double P3);
  void SetReactionIndex(int RI);
  void SetTracePoint(float t, float x, float y, float z, float K);
  void SetX(double X0, double X1, double X2, double X3);
  void SetTrajectoryAtt(short Color, short Style=1, short Width=1);

  // Assignment operator 
  //  Particle & operator=(const Particle &rhs);

  // Members
  int A;
  double Mass;
  static const int MaxPoints = 50000;
  std::string Name;
  int NEexc;
  int Q;
  int RI;
  bool SaveTrajectory;
  bool DoNotPropagate;
  TEveStraightLineSet* Trajectory;
  int Z;

private:
  // For the traces (time, space coords and number of points).
  float* TrT;
  float* TrX;
  float* TrY;
  float* TrZ;
  float* TrK;
  int TrPts;

  FourVector P;               // Four-momentum
  FourVector X;               // Four-position

  // Excitation energy
  double* Eexc;
  // Probability for each excited states to be selected
  double* ProbExc;
  // Probadility distribution function
  TRandom3* PDF;
  
  // Energy loss stuff
  const catima::Material* gas_;
  float dEdxScale_;
  EnergyLoss** IonInMedium;
  int NumMedia;
  static const int MaxMedia = 40;

  int CurrentExcState;

  // For particle trajectories
  short AttColor;
  short AttStyle;
  short AttWidth;
};

#endif

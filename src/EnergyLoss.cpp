/*******************************************************************
Code: EnergyLoss.cpp

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

#include "EnergyLoss.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////
// Constructor. When called with default values, i.e. EnergyLoss(), the SRIM table will 
// need to be loaded by calling LoadSRIMFile and the ion mass by SetIonMass.
////////////////////////////////////////////////////////////////////////////////////////
EnergyLoss::EnergyLoss(string SRIM_file, float IonMass)
{
  BraggCurve = 0;
  dEdx_e = 0;
  dEdx_n = 0;
  Energy_in_range = 1;
  EvD = new TGraph();
  GoodELossFile = 0;
  IonEnergy = 0;
  this->IonMass = IonMass;  // In MeV/c^2
  last_point = 0;
  points = 0;
  TOF = 0;
  if (SRIM_file!="") 
    LoadSRIMFile(SRIM_file);
}


////////////////////////////////////////////////////////////////////////////////////////
// For the case in which only one curve is created.
////////////////////////////////////////////////////////////////////////////////////////
void EnergyLoss::GetBraggCurve(float InitE, int NSteps, float* Dist, float StepSize)
{
  float InitialEnergy[1] = {InitE};
  GetBraggCurves(1, InitialEnergy, NSteps, Dist, StepSize);
  return;
}



void EnergyLoss::GetBraggCurves(int NCurves, float* InitE, int NSteps, float FinalDist)
{
  float CurrentE=0, DeltaE=0, Dist=FinalDist/NSteps, IntegrationStep=0.05;
  BraggCurve = new TGraph*[NCurves];
  for (int c=0; c<NCurves; c++) {
    BraggCurve[c] = new TGraph();
    CurrentE = InitE[c];
    for (int s=0; s<NSteps; s++) {
      // Integrate the energy loss in this step (it may be a large step, like in the 
      // detector MUSIC)
      DeltaE = CurrentE - GetFinalEnergy(CurrentE, Dist, IntegrationStep);
      if (Energy_in_range) {
	BraggCurve[c]->SetPoint(s, s*Dist, DeltaE);
	CurrentE -= DeltaE;
      }
      else
	break;
    }
  }
  return;
}


////////////////////////////////////////////////////////////////////////////////////////
// This is another version of the GetBraggCurves method in which the user specifies the lenght
// of each step in the elements of Dist array.  The lenght of the steps is then Dist[0], Dist[1],
// etc as opposed to the previous version in which all the steps have a fixed length equal to
// FinalDist/NSteps.
////////////////////////////////////////////////////////////////////////////////////////
void EnergyLoss::GetBraggCurves(int NCurves, float* InitE, int NSteps, float* Dist, float StepSize)
{
  float CurrentE, DeltaE, CurrentPos;
  BraggCurve = new TGraph*[NCurves];
  for (int c=0; c<NCurves; c++) {
    BraggCurve[c] = new TGraph();
    CurrentE = InitE[c];
    CurrentPos = 0;
    for (int s=0; s<NSteps; s++) {
      // Integrate the energy loss in this step (it may be a large step, like in the 
      // detector MUSIC)
      DeltaE = CurrentE - GetFinalEnergy(CurrentE, Dist[s], StepSize);  // Integrating over steps.
      if (Energy_in_range) {
	BraggCurve[c]->SetPoint(2*s, CurrentPos, DeltaE);
	BraggCurve[c]->SetPoint(2*s+1, CurrentPos+Dist[s], DeltaE);
	CurrentE -= DeltaE;
	CurrentPos += Dist[s];
      }
      else {
	DeltaE = CurrentE;
	BraggCurve[c]->SetPoint(2*s, CurrentPos, DeltaE);
	BraggCurve[c]->SetPoint(2*s+1, CurrentPos+Dist[s], DeltaE);
	break;
      }
    }
  }
  return;
}


////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////
void EnergyLoss::GetBraggCurves(int NCurves, float InitE, int* NSteps, float** Dist, float StepSize)
{
  float CurrentE, DeltaE, CurrentPos;
  BraggCurve = new TGraph*[NCurves];
  for (int c=0; c<NCurves; c++) {
    BraggCurve[c] = new TGraph();
    CurrentE = InitE;
    CurrentPos=0;
    for (int s=0; s<NSteps[c]; s++) {
      // Integrate the energy loss in this step (it may be a large step, like in the 
      // detector MUSIC)
      DeltaE = CurrentE - GetFinalEnergy(CurrentE, Dist[c][s], StepSize);  // Integrating over steps.
      if (Energy_in_range) {
	BraggCurve[c]->SetPoint(s, CurrentPos+Dist[c][s], DeltaE);
	CurrentE -= DeltaE;
	CurrentPos += Dist[c][s];
      }
      else
	break;
    }
  }
  return;
}


////////////////////////////////////////////////////////////////////////////////////////
// Similarly as in the GetBraggCurves methods, here the EnergyCurves (TGraph pointers)
// are created but the y-axis is now the remaining energy at a given distance. 
////////////////////////////////////////////////////////////////////////////////////////
void EnergyLoss::GetEnergyCurves(int NCurves, float* InitE, int NSteps, float* Dist, float StepSize)
{
  float CurrentE, CurrentPos;
  EnergyCurve = new TGraph*[NCurves];
  for (int c=0; c<NCurves; c++) {
    EnergyCurve[c] = new TGraph();
    CurrentE = InitE[c];
    CurrentPos=0;
    for (int s=0; s<NSteps; s++) {
      // Integrate the energy loss in this step (it may be a large step, like in the 
      // detector MUSIC)
      CurrentE = GetFinalEnergy(CurrentE, Dist[s], StepSize);  // Integrating over steps.
      if (Energy_in_range) {
	EnergyCurve[c]->SetPoint(s, CurrentPos+Dist[s], CurrentE);
	CurrentPos += Dist[s];
      }
      else
	break;
    }
  }
  return;
}


////////////////////////////////////////////////////////////////////////////////////////
// For the (simpler) case in which only one curve is created.
////////////////////////////////////////////////////////////////////////////////////////
void EnergyLoss::GetEnergyCurve(float InitE, int NSteps, float* Dist, float StepSize)
{
  float InitialEnergy[1] = {InitE};
  GetEnergyCurves(1, InitialEnergy, NSteps, Dist, StepSize);
  return;
}


////////////////////////////////////////////////////////////////////////////////////////
// Get the energy loss of the ion in the gas for a given ion's energy and a differential
// distance through the gas target. Most important function in this class. 
////////////////////////////////////////////////////////////////////////////////////////
double EnergyLoss::GetEnergyLoss(double energy /*MeV*/, double distance /*cm*/)
{
  int i = -1;
  double dEdx = 0;
  if (Energy_in_range) {
    // Look for two points for which the initial energy lies between them.
    // This for-loop should find the points unless there was a big jump from
    // the energy used in the last point and the energy used now.
    for(int p=last_point-1; p<points-1; p++){
      if (last_point>=0)
	if(energy>=IonEnergy[p]  && energy<IonEnergy[p+1]){
	  i = p+1;
	  last_point = p;
	  break;
	}
    }
    // It is probable that if the point wasn't found could have been because of
    // a big jump in the energy (see above), so we need to look in the remaining
    // points.
    if (i==-1) {
      for(int p=0; p<last_point-1; p++){
	if(energy>=IonEnergy[p]  && energy<IonEnergy[p+1]){
	  i = p+1;
	  last_point = p;
	  break;
	}
      }
    }
    // If after the last two for-loops 'i' is still -1 it means the energy was out of range.
    if(i==-1){
      //cout << "*** EnergyLoss warning: energy not within range: " << energy << endl;
      Energy_in_range = 0;
      return 0;
    }
    // If the initial energy is within the range of the function, get the stopping power
    // for the initial energy.  
    double E1 = IonEnergy[i-1];        double E2 = IonEnergy[i];
    double dEdx_e1 = dEdx_e[i-1];      double dEdx_e2 = dEdx_e[i];  
    double dEdx_n1 = dEdx_n[i-1];      double dEdx_n2 = dEdx_n[i];  
    // Interpolating the electric stopping power (from point 1 to 'e').
    double dEdx_e1e = dEdx_e1 + (energy - E1)*(dEdx_e2 - dEdx_e1)/(E2 - E1);
    // Interpolating the nuclear stopping power (usually negligable).
    double dEdx_n1e = dEdx_n1 + (energy - E1)*(dEdx_n2 - dEdx_n1)/(E2 - E1);
    // The stopping power units are in MeV/mm so we multiply by 10 to convert to MeV/cm.
    dEdx = (dEdx_e1e+dEdx_n1e)*10*distance;
  }  
  return (1.15*dEdx);  // Warning: For 17F(alpha,p) only!!
}

////////////////////////////////////////////////////////////////////////////////////////
// Float version of this method
////////////////////////////////////////////////////////////////////////////////////////
float EnergyLoss::GetEnergyLoss(float energy /*MeV*/, float distance /*cm*/)
{
  float dEdx = (float)GetEnergyLoss((double)energy, (double) distance);
  return dEdx;
}



////////////////////////////////////////////////////////////////////////////////////////
// This function sets the points of the energy v. distance TGraph object (called EvD)
// for a given initial energy, going from distance 0 to a final distance (depth) in a
// given number of steps.  The EvD object can then be diplayed by 
//   EL->EvD->Draw("l");
// assuming EL is an EnergyLoss pointer.
////////////////////////////////////////////////////////////////////////////////////////
void EnergyLoss::GetEvDCurve(float InitEne/*MeV*/, float FinalDepth/*cm*/, int steps)
{
  double current_ene, current_depth;
  last_point = 0;
  Energy_in_range = 1;
  current_ene = InitEne;
  current_depth = 0;
  for(int s=0; s<steps; s++){
    current_ene -= GetEnergyLoss(current_ene, (double)(FinalDepth/steps));
    current_depth += FinalDepth/steps;
    EvD->SetPoint(s, current_depth, current_ene);
  }
}


////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////
double EnergyLoss::GetInitialEnergy(float FinalEnergy /*MeV*/, float PathLength /*cm*/, 
				    float StepSize/*cm*/)
{
  double Energy = FinalEnergy;
  int Steps = (int)floor(PathLength/StepSize);
  last_point = 0;
  // The function starts by assuming FinalEnergy is within the energy range, but
  // this could be changes in the GetEnergyLoss() function.
  Energy_in_range = 1;

  for (int s=0; s<Steps; s++) {
    Energy = Energy + GetEnergyLoss(Energy, (double)(PathLength/Steps));
    if (!Energy_in_range)
      break;
  } 
  Energy = Energy + GetEnergyLoss(Energy, (double)(PathLength-Steps*StepSize));
  
  // If the energy is out of the stopping power table range, return an
  // unrealistic energy value.
  if (!Energy_in_range)
    Energy = -1000;

  return Energy;
}


////////////////////////////////////////////////////////////////////////////////////////
// 
////////////////////////////////////////////////////////////////////////////////////////
double EnergyLoss::GetFinalEnergy(double InitialEnergy /*MeV*/, double PathLength /*cm*/, 
				  double StepSize/*cm*/)
{
  double c = 29.9792458;  // Speed of light in cm/ns
  double Energy = InitialEnergy;
  int Steps = (int)floor(PathLength/StepSize);
  last_point = 0;
  TOF = 0;

  // The function starts by assuming InitialEnergy is within the energy range, but
  // this could be changes in the GetEnergyLoss() function.
  Energy_in_range = 1;

  // Loop over the steps
  for (int s=0; s<Steps; s++) {
    if (IonMass>0 && Energy>0) {
      //      double vel = sqrt(2*Energy/IonMass)*c;    // non-relativistic
      double vel = c*sqrt(1 - pow(IonMass/(IonMass+Energy),2));
      double DeltaT = StepSize/vel;
      TOF += DeltaT;
    }
    Energy = Energy - GetEnergyLoss(Energy, (double)(PathLength/Steps));
    if (!Energy_in_range)
      break;
  }
  // After the loop, there may by a very small distance to be covered
  // by the ion. The following lines account for this distance.
  if (IonMass>0 && Energy>0) {
    double vel = sqrt(2*Energy/IonMass);
    double DeltaT = (PathLength-Steps*StepSize)/vel;
    TOF += DeltaT;
  }
  Energy = Energy - GetEnergyLoss(Energy, (double)(PathLength-Steps*StepSize));
 
  // If the energy is out of the stopping power table range, return an
  // unrealistic energy value.
  if (!Energy_in_range) 
    Energy = -1000;
  return Energy;
}

////////////////////////////////////////////////////////////////////////////////////////
// Float version of this method
////////////////////////////////////////////////////////////////////////////////////////
float EnergyLoss::GetFinalEnergy(float InitialEnergy /*MeV*/, float PathLength /*cm*/, 
				 float StepSize/*cm*/)
{
  float Energy = (float)GetFinalEnergy((double)InitialEnergy, (double)PathLength, (double)StepSize);
  return Energy;
}


////////////////////////////////////////////////////////////////////////////////////////
// For a given incident energy, this method returns the optimum step size (in cm).
////////////////////////////////////////////////////////////////////////////////////////
double EnergyLoss::GetOptimumStepSize(float Energy /*MeV*/)
{
  double StepSize = Energy/GetEnergyLoss((double)Energy, (double)1.0);
  return StepSize;
}


////////////////////////////////////////////////////////////////////////////////////////
// Calulates the ion's path length in cm.
////////////////////////////////////////////////////////////////////////////////////////
double EnergyLoss::GetPathLength(float InitialEnergy /*MeV*/, float FinalEnergy /*MeV*/,
				 float DeltaT /*ns*/)
{
  double c = 29.9792458;  // Speed of light in cm/ns
  double L = 0, DeltaX = 0;
  double Kn = InitialEnergy;
  int n = 0;
  int max_n = (int)(10E6);
  last_point = 0;
  Energy_in_range = 1;

  if (IonMass==0)
    cout << "*** EnergyLoss Error: Path length cannot be calculated for IonMass = 0." << endl;
  else {
    // The path length (L) is proportional to sqrt(Kn). After the sum, L will be multiplied by
    // the proportionality factor.
    while (Kn>FinalEnergy && n<max_n) {
      L += sqrt(Kn);                       // DeltaL going from point n to n+1.
      DeltaX = sqrt(2*Kn/IonMass)*DeltaT*c;// dx = v*dt
      Kn -= GetEnergyLoss(Kn, DeltaX);     // After increasing L the energy at n+1 is calculated.
      n++;    
    }
    if (n>=max_n) {
      cout << "*** EnergyLoss Warning: Full path length wasn't reached after " << max_n 
	   << " iterations." << endl;
      L = 0;
    }
    else
      L *= sqrt(2/IonMass)*DeltaT*c;
  }
  return L;
}



////////////////////////////////////////////////////////////////////////////////////////
// Calulates the ion's time of flight in ns.
////////////////////////////////////////////////////////////////////////////////////////
double EnergyLoss::GetTimeOfFlight(float InitialEnergy /*MeV*/, float PathLength /*cm*/, 
				   float StepSize /*cm*/)
{
  double c = 29.9792458;  // Speed of light in cm/ns
  double TOF = 0;
  double Kn = InitialEnergy;
  int Steps = PathLength/StepSize;
  last_point = 0;
  Energy_in_range = 1;

  if (IonMass==0)
    cout << "*** EnergyLoss Error: Time of flight cannot be calculated for IonMass = 0." << endl;
  else {
    // The TOF is proportional to 1/sqrt(Kn). After the sum, TOF will be multiplied by
    // the proportionality factor.
    for (int n=0; n<Steps; n++) {
      TOF += 1/sqrt(Kn);                 // DeltaT going from point n to n+1
      Kn -= GetEnergyLoss(Kn, (double)StepSize); // Then, the kinetic energy at point n+1 is calculated
      if (Kn<0)
	break;
    }
    TOF *= sqrt(IonMass/2)*StepSize/c;
  }
  return TOF;
}

////////////////////////////////////////////////////////////////////////////////////////
// Returns the time of flight computed in the GetFinalEnergy method.
////////////////////////////////////////////////////////////////////////////////////////
double EnergyLoss::GetTimeOfFlight()
{
  return TOF;
}



////////////////////////////////////////////////////////////////////////////////////////
// Reads the stopping power as function of incident ion energy from a table generated
// with the program SRIM.
////////////////////////////////////////////////////////////////////////////////////////
bool EnergyLoss::LoadSRIMFile(string FileName)
{
  double IonEnergy, dEdx_e, dEdx_n;
  string aux, unit;
  int str_counter=0;
  ifstream Read(FileName.c_str());
  this->FileName = FileName;
  last_point = 0;
  if(!Read.is_open()) {
    cout << "*** EnergyLoss Error: File " << FileName << " was not found." << endl;
    GoodELossFile = 0;
  } 
  else {
    GoodELossFile = 1;
    Energy_in_range = 1;        
    // Read all the string until you find "Straggling", then read the next 7 strings.
    // (this method is not elegant at all but there is no time to make it better)
    do 
      Read >> aux;
    while (aux!="Straggling");
    for (int i=0; i<7; i++)
      Read >> aux;

    do {
      Read >> aux;
      str_counter++;
    } while (aux!="-----------------------------------------------------------");
    Read.close();
   
    str_counter--;
    points = str_counter/10;
   
    // Create the arrays depending on the number rows in the file.
    this->IonEnergy = new double[points];
    this->dEdx_e = new double[points];
    this->dEdx_n = new double[points];    
 
    // Go to the begining of the file and read it again to now save the info in the
    // newly created arrays.
    Read.open(FileName.c_str());
    do 
      Read >> aux;
    while (aux!="Straggling");

    for (int i=0; i<7; i++)
      Read >> aux;
    
    for (int p=0; p<points; p++) {
      Read >> IonEnergy >> unit >> dEdx_e >> dEdx_n >> aux >> aux >> aux >> aux >> aux >> aux;
      if (unit=="eV")
	IonEnergy *= 1E-6;
      else if (unit=="keV")
	IonEnergy *= 0.001;
      else if (unit=="GeV")
	IonEnergy *= 1000;
      //cout << p << " " << IonEnergy << " " << unit << " " << dEdx_e << " " << dEdx_n << endl;
      this->IonEnergy[p] = IonEnergy;
      this->dEdx_e[p] = dEdx_e;
      this->dEdx_n[p] = dEdx_n;
    }    
  }
  return GoodELossFile;
}


////////////////////////////////////////////////////////////////////////////////////////
// IonMass is used in methods GetPathLength and GetTimeOfFlight.
////////////////////////////////////////////////////////////////////////////////////////
void EnergyLoss::SetIonMass(float IonMass)
{
  this->IonMass = IonMass;
  return;  
}

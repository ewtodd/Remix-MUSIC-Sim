// Methods for the Particle class. See header file for class description.
// Compile with: 
// On linux
//  g++ -shared -fPIC Particle.cpp `root-config --cflags --glibs` -o Particle.so
// On MacOS (work in progress)
// g++ -shared -fPIC Particle.cpp FourVector.so EnergyLoss.so `root-config --cflags --glibs` -lEve -o Particle.so
// By Daniel Santiago-Gonzalez
// 2014

#include "Particle.hpp"
#include <cmath>

using namespace std;

///////////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////////
Particle::Particle(string Name, double M, int Q, bool SaveTrajectory)
{
  this->Name = Name;
  this->Q = Q;
  Mass = M;
  A = 0;  // Not used at this moment.
  Z = Q;
  NEexc = 1;
  this->SaveTrajectory = SaveTrajectory;
  DoNotPropagate = false;
  
  // Variables
  AttColor = 1;
  AttStyle = 1;
  AttWidth = 2;
  CurrentExcState = 0;
  NumMedia = 0;
  P.SetName("P_"+Name);
  P.SetCoords(M,0,0,0);
  RI = 0;
  TrPts = 0;
  X.SetName("X_"+Name);
  X.SetCoords(0,0,0,0);
  // Pointers
  Eexc = new double[NEexc];
  Eexc[0] = 0.0;
  PDF = 0;
  ProbExc = 0;
  IonInMedium = new EnergyLoss*[MaxMedia];
  for (int m=0; m<MaxMedia; m++)
    IonInMedium[m] = 0;
  gas_ = nullptr;
  dEdxScale_ = 1.0;
  TrT = new float[MaxPoints];
  TrX = new float[MaxPoints];
  TrY = new float[MaxPoints];
  TrZ = new float[MaxPoints];
  TrK = new float[MaxPoints];
  // Particle trajectories
  Trajectory = new TEveStraightLineSet();
  AllTraj = new TEveStraightLineSet*[MaxEvents];
  for (int e=0; e<MaxEvents; e++) 
    AllTraj[e] = new TEveStraightLineSet();
}


///////////////////////////////////////////////////////////////////////////////////
// Lorentz boost for the position and momentum four-vectors.
///////////////////////////////////////////////////////////////////////////////////
void Particle::Boost(double BetaX, double BetaY, double BetaZ)
{
  X.Boost(BetaX, BetaY, BetaZ);
  P.Boost(BetaX, BetaY, BetaZ);
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Copy members of another Particle object.  It does not modify Name. Needs to be
// updated.
///////////////////////////////////////////////////////////////////////////////////
//Particle& Particle::operator=(const Particle& rhs) {
void Particle::Copy(Particle* rhs) {
  if (this != rhs) {
    // Do the assignment operation.
    A = rhs->A;
    SetExcEnergies(rhs->NEexc, rhs->Eexc);
    Mass = rhs->Mass;
    Q = rhs->Q;
    Z = rhs->Z;
    X = rhs->X;
    P = rhs->P;
    if (rhs->gas_) SetMedium(rhs->gas_, rhs->dEdxScale_);
  }
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Copy the trace information to individual arrays.
///////////////////////////////////////////////////////////////////////////////////
void Particle::CopyTrace(int& NumPts, float* t, float* x, float* y, float* z, float* K)
{
  int TotPoints = TrPts;
  if (TrPts>MaxPoints)
    TotPoints = MaxPoints;
  for (int p=0; p<TotPoints; p++) {
    t[p] = TrT[p];
    x[p] = TrX[p];
    y[p] = TrY[p];
    z[p] = TrZ[p];
    K[p] = TrK[p];
  }
  NumPts = TotPoints;
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Get the cartesian components of the velocity in units of c.
///////////////////////////////////////////////////////////////////////////////////
void Particle::GetBeta(double& BetaX, double& BetaY, double& BetaZ)
{
  double E = P.GetX0();
  double px = P.GetX1();
  double py = P.GetX2();
  double pz = P.GetX3();
  BetaX = px/E;
  BetaY = py/E;
  BetaZ = pz/E;
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Used only in ANASEN Simulator and HELIOS Simulator.
///////////////////////////////////////////////////////////////////////////////////
int Particle::GetCurrentExcState()
{
  return CurrentExcState;
}


///////////////////////////////////////////////////////////////////////////////////
// Get an excitation energy based on the probaility distribution, ExcProb, which
// was specified by the user in SetExcEnergies().
///////////////////////////////////////////////////////////////////////////////////
double Particle::GetEexc()
{
  double Eexc = 0;
  if (NEexc==1)
    Eexc = this->Eexc[0];
  else if (PDF!=0 && ProbExc!=0 && this->Eexc!=0) {
    double Rdm = PDF->Uniform();
    for (int ne=0; ne<NEexc; ne++) {
      if (Rdm>=ProbExc[ne] && Rdm<ProbExc[ne+1]) {
	CurrentExcState = ne;
	Eexc = this->Eexc[ne];
	break;
      }
    }
  }
  // Working to make this part of the method obsolete and change to
  // use ExcProb only.
  // if (this->Eexc!=0 && CurrentExcState>=0 && CurrentExcState<NEexc)
  //   Eexc = this->Eexc[CurrentExcState];
  return Eexc;
}

double Particle::GetEexc(int ExcState)
{
  double Eexc = 0;
  if (this->Eexc!=0 && ExcState>=0 && ExcState<NEexc)
    Eexc = this->Eexc[ExcState];
  return Eexc;
}


///////////////////////////////////////////////////////////////////////////////////
// OBSOLETE!! 
// Returns the energy loss for a given medium (MediumID) and initial energy (InitE) 
// over a small (differential) path length.
///////////////////////////////////////////////////////////////////////////////////
double Particle::GetEnergyLoss(int MediumID, double InitE/*MeV*/, double PathLength/*cm*/)
{
  double DE = 0;
  //  int m = MediumID;
  // if (m>=0 && m<NumMedia) {
  //   DE = IonInMedium[m]->GetEnergyLoss(InitE, PathLength);   
  //   if (DE<0)
  //     DE = 0;
  // }
  return DE;
}



///////////////////////////////////////////////////////////////////////////////////
// Compute final energy from a given medium ID, initial energy, path length and 
// step size. Wrapped method of EnergyLoss class.
///////////////////////////////////////////////////////////////////////////////////
double Particle::GetFinalEnergy(int MediumID, double InitE/*MeV*/, double PathLength/*cm*/,
				double StepSize/*cm*/)
{
  double FinalE = InitE;
  int m = MediumID;
  if (m>=0 && m<NumMedia) {
    FinalE = IonInMedium[m]->GetFinalEnergy(InitE, PathLength, StepSize);   
    if (FinalE<0)
      FinalE = 0;
  }
  return FinalE;
}


///////////////////////////////////////////////////////////////////////////////////
// Compute initial energy from a given medium ID, final energy, path length and 
// step size. Wrapped method of EnergyLoss class
///////////////////////////////////////////////////////////////////////////////////
double Particle::GetInitialEnergy(int MediumID, double FinalE/*MeV*/, double PathLength/*cm*/,
				  double StepSize/*cm*/)
{
  double InitE = FinalE;
  int m = MediumID;
  if (m>=0 && m<NumMedia) {
    InitE = IonInMedium[m]->GetInitialEnergy(FinalE, PathLength, StepSize);   
    if (InitE<0)
      InitE = 0;
  }
  return InitE;
}


///////////////////////////////////////////////////////////////////////////////////
// Return the kinetic energy (total energy minus mass minus excitation energy).
///////////////////////////////////////////////////////////////////////////////////
double Particle::GetKE()
{
  double KE = P.GetX0() - Mass - GetEexc();
  if (KE<0.0) {
    KE = 0.0;
    cout << "Warning(" << Name << "): unrealistic negative kinetic energy (KE) -> KE set equal to 0."
	 << endl;
  }
  return KE;
}


///////////////////////////////////////////////////////////////////////////////////
// For a given incident energy, this method returns the optimum step size (in cm).
///////////////////////////////////////////////////////////////////////////////////
double Particle::GetOptimumStepSize(int MediumID, double Energy/*MeV*/)
{
  double StepSize = 0.1;
  int m = MediumID;
  if (m>=0 && m<NumMedia) 
    StepSize = 0.01*(IonInMedium[m]->GetOptimumStepSize(Energy));
  return StepSize;
}




///////////////////////////////////////////////////////////////////////////////////
// Return four momentum objects.
///////////////////////////////////////////////////////////////////////////////////
FourVector Particle::GetP()
{
  FourVector P = this->P;
  return P;
}


///////////////////////////////////////////////////////////////////////////////////
// Return coordinates of four momentum.
///////////////////////////////////////////////////////////////////////////////////
void Particle::GetP(double& P0, double& P1, double& P2, double& P3)
{
  P0 = P.GetX0();
  P1 = P.GetX1();
  P2 = P.GetX2();
  P3 = P.GetX3();
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////////
double Particle::GetPathLength(int MediumID, double InitE/*MeV*/, double FinalE/*MeV*/,
			       double DeltaT/*ns*/)
{
  double PL = 0;
  int m = MediumID;
  if (m>=0 && m<NumMedia) 
    PL = IonInMedium[m]->GetPathLength((float)InitE, (float)FinalE, (float)DeltaT);
  return PL;
}


///////////////////////////////////////////////////////////////////////////////////
// Get the azimuthal angle of the momentum vector in radians (from 0 to 2*pi).
///////////////////////////////////////////////////////////////////////////////////
double Particle::GetPhi()
{
  double pi = 3.14159265359;
  double phi = 0;
  double px = P.GetX1();
  double py = P.GetX2();
  if (px>=0 && py>0) 
    phi = atan(py/px);
  else if (px<0 && py>0)
    phi = pi + atan(py/px);
  else if (px<0 && py<0)
    phi = pi + atan(py/px);
  else if (px>0 && py<0)
    phi = 2*pi + atan(py/px);

  return phi;
}

///////////////////////////////////////////////////////////////////////////////////
// Get the azimuthal angle of the position vector in radians (from 0 to 2*pi).
///////////////////////////////////////////////////////////////////////////////////
double Particle::GetPhiX()
{
  double pi = 3.14159265359;
  double phi = 0;
  double x = X.GetX1();
  double y = X.GetX2();
  if (x>=0 && y>0) 
    phi = atan(y/x);
  else if (x<0 && y>0)
    phi = pi + atan(y/x);
  else if (x<0 && y<0)
    phi = pi + atan(y/x);
  else if (x>0 && y<0)
    phi = 2*pi + atan(y/x);

  return phi;
}


///////////////////////////////////////////////////////////////////////////////////
// Get the polar angle of the momentum vector in radians.
///////////////////////////////////////////////////////////////////////////////////
double Particle::GetTheta()
{
  double pi = 3.14159265359;
  double theta;
  double px = P.GetX1();
  double py = P.GetX2();
  double pz = P.GetX3();
  theta = atan(sqrt(px*px+py*py)/pz);
  if (pz<0)
    theta += pi;
  return theta;
}

///////////////////////////////////////////////////////////////////////////////////
// Get the polar angle of the position vector in radians.
///////////////////////////////////////////////////////////////////////////////////
double Particle::GetThetaX()
{
  double pi = 3.14159265359;
  double theta;
  double x = X.GetX1();
  double y = X.GetX2();
  double z = X.GetX3();
  theta = atan(sqrt(x*x+y*y)/z);
  if (z<0)
    theta += pi;
  return theta;
}



///////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////
double Particle::GetTimeOfFlight(int MediumID)
{
  double TOF = 0;
  int m = MediumID;
  if (m>=0 && m<NumMedia) 
    TOF = IonInMedium[m]->GetTimeOfFlight();
  return TOF;
}

///////////////////////////////////////////////////////////////////////////////////
// Calculate the particle's time of flight (in ns) based on the medium stopping 
// power, initial energy and track lenght.
///////////////////////////////////////////////////////////////////////////////////
double Particle::GetTimeOfFlight(int MediumID, float InitialEnergy/*MeV*/, float PathLength/*cm*/,
				 float StepSize/*cm*/)
{
  double TOF = 0;
  int m = MediumID;
  if (m>=0 && m<NumMedia) 
    TOF = IonInMedium[m]->GetTimeOfFlight(InitialEnergy, PathLength, StepSize);
  return TOF;
}


///////////////////////////////////////////////////////////////////////////////////
// Return the trajectory attributes.
///////////////////////////////////////////////////////////////////////////////////
void Particle::GetTrajectoryAtt(short& Color, short& Style, short& Width)
{
  Color = AttColor;
  Style = AttStyle;
  Width = AttWidth;
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Return coordinates of position four vector.
///////////////////////////////////////////////////////////////////////////////////
void Particle::GetX(double& X0, double& X1, double& X2, double& X3)
{
  X0 = X.GetX0();
  X1 = X.GetX1();
  X2 = X.GetX2();
  X3 = X.GetX3();
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Simple function that prints some variables of this class.
///////////////////////////////////////////////////////////////////////////////////
void Particle::Print(ostream& log)
{
  log << "|== Particle " << Name << " =================================|"<< endl;
  log << "| mass = " << Mass << " MeV/c^2   Z = " << Q << " e\n"
      << "| Eexc = ";
  if (NEexc>0) {
    for (int n=0; n<NEexc; n++) {
      log << Eexc[n];
      if (n<NEexc-1)
	log << ", ";
      else
	log << "\n";
    }
  }
  else
    log << " 0 MeV" << endl;
  log << "| KE = " << GetKE() << " Mev" << endl;
  log << "| ";
  X.Print(log);
  log << "| ";
  P.Print(log);
  if (SaveTrajectory) {
    log << "| Trajectory: " << Trajectory  << " C=" << Trajectory->GetLineColor() 
	<< " W=" << Trajectory->GetLineWidth() << " S=" << Trajectory->GetLineStyle() << endl;
  } 
  else 
    log << "| Trajectory object not saved." << endl;
  if (NumMedia>0 && gas_) {
    log << "| Stopping-power via catima; gas density = "
        << gas_->density() << " g/cm³, dE/dx scale = " << dEdxScale_ << endl;
  }
  log << "|==================================================|"<< endl;
  return;
}


///////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////
void Particle::ResetTrace()
{
  int TotPoints = TrPts;
  if (TrPts==0)
    TotPoints = MaxPoints;
  for (int p=0; p<TotPoints; p++) {
    TrT[p] = -1000;
    TrX[p] = 0;
    TrY[p] = 0;
    TrZ[p] = -1000;
    TrK[p] = 0;
  }
  TrPts = 0;
  return;
}


///////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////
void Particle::ResetKinematics()
{
  P.SetCoords(Mass,0,0,0);
  X.SetCoords(0,0,0,0);
  Eexc[0] = 0.0;
  return;
}



///////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////
void Particle::SetCurrentExcState(int ExcState)
{
  CurrentExcState = ExcState;
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Set the number of excited states for this particle and obtain their values from 
// the array Eexc[]. The user has the option to provide the probability for each
// excitation energy to be selected. This can be specified in the Prob array. If
// not specified (i.e. Prob=0, default value), all excited states will have the 
// same probability to be selected when calling the GetEexc() method.
// Prob example: if Prob[3] = {1.5, 3.5, 2.5}, the cumulative probability is given 
// by, ProbExc[4] = {0, 1.5, 5.0, 7.5}. Then, the final normalized array is
// ProbExc[4] = {0, 0.2, 0.666, 1}
///////////////////////////////////////////////////////////////////////////////////
void Particle::SetExcEnergies(int N, double* Eexc, double* Prob)
{
  double Norm = 0;
  double Cumulative = 0;
  if (N>0 && Eexc!=0) {
    NEexc = N;
    this->Eexc = new double[N];
    //    this->Prob = new double[N];
    ProbExc = new double[N+1];
    // Set the excitation energies and get the normalization factor
    // for ProbExc
    for (int n=0; n<N; n++) {
      this->Eexc[n] = Eexc[n];
      if (Prob!=0)
	Norm += Prob[n];
      else
	Norm += 1.0/N;
    }
    // Get the cumulative probability
    ProbExc[0] = 0.0;
    for (int n=0; n<N; n++) {
      if (Prob!=0)
	Cumulative += Prob[n];
      else
	Cumulative += 1.0/N;
      ProbExc[n+1] = Cumulative;
    }
    // Normalize ProbExc
    for (int n=0; n<N+1; n++) 
      ProbExc[n] = ProbExc[n]/Norm;

    // Pseudo-random number generator
    PDF = new TRandom3();
    PDF->SetSeed();
  }
  return;
}



///////////////////////////////////////////////////////////////////////////////////
// This method is used when the user wants to specify a single excitation energy.
///////////////////////////////////////////////////////////////////////////////////
void Particle::SetExcEnergy(double Ex) {
  Eexc[0] = Ex;
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Configure energy loss in the (single) gas medium using catima.
// A is derived from Mass (MeV/c²) / atomic mass unit; Z is the particle's Z.
///////////////////////////////////////////////////////////////////////////////////
void Particle::SetMedium(const catima::Material* gas, float dEdxScale)
{
  if (gas == nullptr) return;
  if (Z <= 0) {
    // Expected for neutral particles (e.g. neutrons); skip silently. They
    // are propagated by ExitWindow / kinematics-only paths, not catima.
    return;
  }
  const double amu_MeV = 931.49410242;
  int A_derived = (Mass > 0.0) ? int(std::round(Mass / amu_MeV)) : 0;
  if (A_derived <= 0) {
    cout << Name << ": cannot derive mass number from Mass=" << Mass << " MeV/c²" << endl;
    return;
  }
  if (A == 0) A = A_derived;
  NumMedia = 1;
  gas_ = gas;
  dEdxScale_ = dEdxScale;
  if (IonInMedium[0]) delete IonInMedium[0];
  IonInMedium[0] = new EnergyLoss(A_derived, Z, Mass, gas, dEdxScale);
  return;
}

///////////////////////////////////////////////////////////////////////////////////
// Copy just the coordinates of the four vector V to the class member P.
///////////////////////////////////////////////////////////////////////////////////
void Particle::SetP(FourVector V) 
{
  double P0, P1, P2, P3;
  P0 = V.GetX0();
  P1 = V.GetX1();
  P2 = V.GetX2();
  P3 = V.GetX3();
  P.SetCoords(P0, P1, P2, P3);
  return;
}


///////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////
void Particle::SetP(double P0, double P1, double P2, double P3) 
{
  P.SetCoords(P0, P1, P2, P3);
  return;
}


///////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////
void Particle::SetReactionIndex(int RI)
{
  this->RI = RI;
  return;
}


///////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////
void Particle::SetTracePoint(float t, float x, float y, float z, float K)
{
  int p = TrPts;
  if (p<MaxPoints) {
    TrT[p] = t;
    TrX[p] = x;
    TrY[p] = y;
    TrZ[p] = z;
    TrK[p] = K;
  }
  else 
    cout << "Warning: " << Name << " reached maximum number of trace points." << endl;
  TrPts++;
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// This method is needed because the attributes of TEveStraightLineSet elements
// cannot be set until lines have been added to it (via the AddLine method). 
// Otheriwse, it gives a 'segmentation violation' error.
///////////////////////////////////////////////////////////////////////////////////
void Particle::SetTrajectoryAtt(short Color, short Style, short Width)
{
  AttColor = Color;
  AttStyle = Style;
  AttWidth = Width;
  return;
}


///////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////
void Particle::SetX(double X0, double X1, double X2, double X3) 
{
  X.SetCoords(X0, X1, X2, X3);
  return;
}

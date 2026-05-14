

#include "FourVector.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////
// Empty constructor
///////////////////////////////////////////////////////////////////////////////////
FourVector::FourVector()
{
  Name = "";
  x[0] = x[1] = x[2] = x[3] = 0;
}


///////////////////////////////////////////////////////////////////////////////////
// Full constructor
///////////////////////////////////////////////////////////////////////////////////
FourVector::FourVector(string Name, double x0, double x1, double x2, double x3)
{
  SetName(Name);
  SetCoords(x0, x1, x2, x3);
}


///////////////////////////////////////////////////////////////////////////////////
// Lorentz boost
///////////////////////////////////////////////////////////////////////////////////
void FourVector::Boost(double BetaX, double BetaY, double BetaZ)
{
  double A[4][4];
  double new_x[4];
  const double B[3] = {BetaX, BetaY, BetaZ};
  const double Beta = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
  const double Gamma = 1/sqrt(1-Beta*Beta);
  // Transformation matrix
  A[0][0] = Gamma;
  for (int i=1; i<4; i++) {
    A[0][i] = A[i][0] = -Gamma*B[i-1];
    for (int j=1; j<4; j++)
      if (Beta==0)
	A[j][i] = A[i][j] = Delta(i,j);
      else
	A[j][i] = A[i][j] = (Gamma-1)*B[i-1]*B[j-1]/(Beta*Beta) + Delta(i,j);
  }

  for (int i=0; i<4; i++) {
    new_x[i] = 0;
    for (int j=0; j<4; j++)
      new_x[i] += A[i][j]*x[j];
  }

  for (int i=0; i<4; i++) 
    x[i] = new_x[i];
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Kronecker delta
///////////////////////////////////////////////////////////////////////////////////
double FourVector::Delta(int i, int j) 
{
  double d = 0;
  if (i==j)
    d = 1;
  return d;
}


///////////////////////////////////////////////////////////////////////////////////
// One-liners
///////////////////////////////////////////////////////////////////////////////////
double FourVector::GetX0() { return x[0]; }
double FourVector::GetX1() { return x[1]; }
double FourVector::GetX2() { return x[2]; }
double FourVector::GetX3() { return x[3]; }
string FourVector::GetName() { return Name; }


///////////////////////////////////////////////////////////////////////////////////
// Get the polar angle of the momentum vector in radians.
///////////////////////////////////////////////////////////////////////////////////
double FourVector::GetTheta()
{
  double px = GetX1();
  double py = GetX2();
  return atan2(sqrt(px*px + py*py), GetX3());
}


///////////////////////////////////////////////////////////////////////////////////
// Set coordinates of four-vector
///////////////////////////////////////////////////////////////////////////////////
void FourVector::SetCoords(double x0, double x1, double x2, double x3)
{
  x[0] = x0;
  x[1] = x1;
  x[2] = x2;
  x[3] = x3;
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Give it a name
///////////////////////////////////////////////////////////////////////////////////
void FourVector::SetName(string Name)
{
  this->Name = Name;
}


///////////////////////////////////////////////////////////////////////////////////
// Print basic info
///////////////////////////////////////////////////////////////////////////////////
void FourVector::Print(ostream& log) 
{
  log << Name << " = (" << x[0] << ", " << x[1] << ", " << x[2] << ", " << x[3] << ")" << endl;
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Overloading operators (cool stuff)
///////////////////////////////////////////////////////////////////////////////////
// Assignment.
// Take a const-reference to the right-hand side of the assignment.
// Return a non-const reference to the left-hand side.
FourVector& FourVector::operator=(const FourVector& rhs) {
  if (this != &rhs) {
    // Do the assignment operation.
    //cout << "Operation =" << endl;
    x[0] = rhs.x[0];
    x[1] = rhs.x[1];
    x[2] = rhs.x[2];
    x[3] = rhs.x[3];
  }
  return *this;  // Return a reference to myself.
}

// Increment.
FourVector& FourVector::operator+=(const FourVector& rhs) 
{
  //cout << "Operation +=" << endl;
  x[0] += rhs.x[0];
  x[1] += rhs.x[1];
  x[2] += rhs.x[2];
  x[3] += rhs.x[3];
  return *this;  // Return a reference to myself.
}

// Addition.
// Add this instance's value to other, and return a new instance with the result.
const FourVector FourVector::operator+(const FourVector& other) const
{
  //cout << "Operation +" << endl;
  FourVector result;
  result.SetCoords(x[0]+other.x[0], x[1]+other.x[1], x[2]+other.x[2], x[3]+other.x[3]);
  return result;
}

// Decrement.
FourVector& FourVector::operator-=(const FourVector& rhs) 
{
  //cout << "Operation -=" << endl;
  x[0] -= rhs.x[0];
  x[1] -= rhs.x[1];
  x[2] -= rhs.x[2];
  x[3] -= rhs.x[3];
  return *this;  // Return a reference to myself.
}

// Subtraction.
// Add this instance's value to other, and return a new instance with the result.
const FourVector FourVector::operator-(const FourVector& other) const 
{
  //cout << "Operation -" << endl;
  FourVector result;
  result.SetCoords(x[0]-other.x[0], x[1]-other.x[1], x[2]-other.x[2], x[3]-other.x[3]);
  return result;
}

// Dot product.
double FourVector::operator*(const FourVector& P) {
  //cout << "Operation * (dot product)" << endl;
  return(x[0]*P.x[0] - x[1]*P.x[1] - x[2]*P.x[2] - x[3]*P.x[3]);
}

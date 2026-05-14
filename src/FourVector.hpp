// References for overloading operators:
// http://courses.cms.caltech.edu/cs11/material/cpp/donnie/cpp-ops.html
// Compile with: 
// g++ -shared -fPIC FourVector.cpp -o FourVector.so
#ifndef FourVector_hpp_INCLUDED   
#define FourVector_hpp_INCLUDED   

#include <iostream>
#include <cmath>
#include <string.h>

class FourVector{
  
public:
  FourVector();
  FourVector(std::string Name, double x0=0, double x1=0, double x2=0, double x3=0);
  
  //  FourVector( const FourVector& other);
  void Boost(double BetaX, double BetaY, double BetaZ);
  double GetX0();
  double GetX1();
  double GetX2();
  double GetX3();
  std::string GetName();
  double GetTheta();
  void SetCoords(double x0, double x1, double x2, double x3);
  void SetName(std::string Name);
  void Print(std::ostream& log=std::cout);
  
  FourVector & operator=(const FourVector &rhs);
  FourVector & operator+=(const FourVector &rhs) ;
  const FourVector operator+(const FourVector &other) const;
  FourVector & operator-=(const FourVector &rhs);
  const FourVector operator-(const FourVector &other) const;
  double operator*(const FourVector &P);

private:
  double Delta(int i, int j);

  std::string Name;
  double x[4];
};

#endif

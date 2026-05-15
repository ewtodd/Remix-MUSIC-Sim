#ifndef FOURVECTOR_HPP
#define FOURVECTOR_HPP

#include <cmath>
#include <iostream>

#include <TROOT.h>
#include <TString.h>

class FourVector {
public:
  FourVector();
  FourVector(TString Name, Double_t x0 = 0, Double_t x1 = 0, Double_t x2 = 0,
             Double_t x3 = 0);

  void Boost(Double_t BetaX, Double_t BetaY, Double_t BetaZ);
  Double_t GetX0();
  Double_t GetX1();
  Double_t GetX2();
  Double_t GetX3();
  TString GetName();
  Double_t GetTheta();
  void SetCoords(Double_t x0, Double_t x1, Double_t x2, Double_t x3);
  void SetName(TString Name);
  void Print(std::ostream &log = std::cout);

  FourVector &operator=(const FourVector &rhs);
  FourVector &operator+=(const FourVector &rhs);
  const FourVector operator+(const FourVector &other) const;
  FourVector &operator-=(const FourVector &rhs);
  const FourVector operator-(const FourVector &other) const;
  Double_t operator*(const FourVector &P);

private:
  Double_t Delta(Int_t i, Int_t j);

  TString Name;
  Double_t x[4];
};

#endif

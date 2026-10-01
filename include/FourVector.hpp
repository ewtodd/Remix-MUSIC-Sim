#ifndef FOURVECTOR_HPP
#define FOURVECTOR_HPP

#include <cmath>
#include <iostream>

#include <TROOT.h>
#include <TString.h>

class FourVector {
public:
  FourVector();
  FourVector(TString vectorName, Double_t x0 = 0, Double_t x1 = 0,
             Double_t x2 = 0, Double_t x3 = 0);

  void Boost(Double_t BetaX, Double_t BetaY, Double_t BetaZ);
  Double_t GetX0() const;
  Double_t GetX1() const;
  Double_t GetX2() const;
  Double_t GetX3() const;
  TString GetName() const;
  Double_t GetTheta() const;
  void SetCoords(Double_t x0, Double_t x1, Double_t x2, Double_t x3);
  void SetName(TString vectorName);
  void Print(std::ostream &log = std::cout) const;

  FourVector &operator=(const FourVector &rhs);
  FourVector &operator+=(const FourVector &rhs);
  const FourVector operator+(const FourVector &other) const;
  FourVector &operator-=(const FourVector &rhs);
  const FourVector operator-(const FourVector &other) const;
  Double_t operator*(const FourVector &P) const;

private:
  static Double_t Delta(Int_t i, Int_t j);

  TString Name;
  Double_t x[4];
};

#endif

#include "FourVector.hpp"

FourVector::FourVector() {
  Name = "";
  x[0] = x[1] = x[2] = x[3] = 0;
}

FourVector::FourVector(TString Name, Double_t x0, Double_t x1, Double_t x2,
                       Double_t x3) {
  SetName(Name);
  SetCoords(x0, x1, x2, x3);
}

void FourVector::Boost(Double_t BetaX, Double_t BetaY, Double_t BetaZ) {
  Double_t A[4][4];
  Double_t new_x[4];
  const Double_t B[3] = {BetaX, BetaY, BetaZ};
  const Double_t Beta = std::sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
  const Double_t Gamma = 1 / std::sqrt(1 - Beta * Beta);
  A[0][0] = Gamma;
  for (Int_t i = 1; i < 4; i++) {
    A[0][i] = A[i][0] = -Gamma * B[i - 1];
    for (Int_t j = 1; j < 4; j++)
      if (Beta == 0)
        A[j][i] = A[i][j] = Delta(i, j);
      else
        A[j][i] = A[i][j] =
            (Gamma - 1) * B[i - 1] * B[j - 1] / (Beta * Beta) + Delta(i, j);
  }

  for (Int_t i = 0; i < 4; i++) {
    new_x[i] = 0;
    for (Int_t j = 0; j < 4; j++)
      new_x[i] += A[i][j] * x[j];
  }
  for (Int_t i = 0; i < 4; i++)
    x[i] = new_x[i];
}

Double_t FourVector::Delta(Int_t i, Int_t j) { return (i == j) ? 1.0 : 0.0; }

Double_t FourVector::GetX0() { return x[0]; }
Double_t FourVector::GetX1() { return x[1]; }
Double_t FourVector::GetX2() { return x[2]; }
Double_t FourVector::GetX3() { return x[3]; }
TString FourVector::GetName() { return Name; }

Double_t FourVector::GetTheta() {
  Double_t px = GetX1();
  Double_t py = GetX2();
  return std::atan2(std::sqrt(px * px + py * py), GetX3());
}

void FourVector::SetCoords(Double_t x0, Double_t x1, Double_t x2, Double_t x3) {
  x[0] = x0;
  x[1] = x1;
  x[2] = x2;
  x[3] = x3;
}

void FourVector::SetName(TString Name) { this->Name = Name; }

void FourVector::Print(std::ostream &log) {
  log << Name << " = (" << x[0] << ", " << x[1] << ", " << x[2] << ", " << x[3]
      << ")" << std::endl;
}

FourVector &FourVector::operator=(const FourVector &rhs) {
  if (this != &rhs) {
    x[0] = rhs.x[0];
    x[1] = rhs.x[1];
    x[2] = rhs.x[2];
    x[3] = rhs.x[3];
  }
  return *this;
}

FourVector &FourVector::operator+=(const FourVector &rhs) {
  x[0] += rhs.x[0];
  x[1] += rhs.x[1];
  x[2] += rhs.x[2];
  x[3] += rhs.x[3];
  return *this;
}

const FourVector FourVector::operator+(const FourVector &other) const {
  FourVector result;
  result.SetCoords(x[0] + other.x[0], x[1] + other.x[1], x[2] + other.x[2],
                   x[3] + other.x[3]);
  return result;
}

FourVector &FourVector::operator-=(const FourVector &rhs) {
  x[0] -= rhs.x[0];
  x[1] -= rhs.x[1];
  x[2] -= rhs.x[2];
  x[3] -= rhs.x[3];
  return *this;
}

const FourVector FourVector::operator-(const FourVector &other) const {
  FourVector result;
  result.SetCoords(x[0] - other.x[0], x[1] - other.x[1], x[2] - other.x[2],
                   x[3] - other.x[3]);
  return result;
}

// Minkowski dot product (+, −, −, −).
Double_t FourVector::operator*(const FourVector &P) {
  return x[0] * P.x[0] - x[1] * P.x[1] - x[2] * P.x[2] - x[3] * P.x[3];
}

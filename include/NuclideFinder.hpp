#ifndef NUCLIDEFINDER_HPP
#define NUCLIDEFINDER_HPP

// Look up nuclide masses, spins, and parities by name ("4He"), by (Z, A), or
// by index. Mass units: "MeV/c^2", "u", or "micro-u".

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <TROOT.h>

class NuclideFinder {
public:
  NuclideFinder();

  Int_t GetA(Int_t Index);
  Int_t GetA(std::string Nuclide);
  Int_t GetArraySize();
  Float_t GetGSspin(Int_t Index);
  Int_t GetGSparity(Int_t Index);
  Double_t GetMass(Int_t Index, std::string Units);
  Double_t GetMass(Int_t atomicNumber, Int_t massNumber, std::string Units);
  Double_t GetMass(std::string Nuclide, std::string Units);
  Int_t GetN(std::string Nuclide);
  std::string GetName(Int_t Index);
  std::string GetName(Int_t atomicNumber, Int_t massNumber);
  Int_t GetZ(Int_t Index);
  Int_t GetZ(std::string Nuclide);
  Int_t IsMeasured(Int_t Index);
  Bool_t Contains(std::string Nuclide);
  Int_t LoadAMEFile(std::string file);
  Int_t LoadNubaseFile(std::string file);

private:
  std::string GetProperName(std::string Nuclide);
  void Init();

  Int_t Size;
  std::vector<Int_t> A;
  std::vector<Int_t> N;
  std::vector<Int_t> Z;
  std::vector<Int_t> Measured;
  std::vector<std::string> Name;
  std::vector<Double_t> Mass; // micro-u
  std::vector<Float_t> GSspin;
  std::vector<Int_t> GSparity;
};

#endif

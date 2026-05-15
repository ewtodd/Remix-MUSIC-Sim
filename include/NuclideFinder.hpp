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

#include <TROOT.h>

class NuclideFinder {
public:
  NuclideFinder();

  Int_t GetA(Int_t Index);
  Int_t GetArraySize();
  Float_t GetGSspin(Int_t Index);
  Int_t GetGSparity(Int_t Index);
  Double_t GetMass(Int_t Index, std::string Units);
  Double_t GetMass(Int_t Z, Int_t A, std::string Units);
  Double_t GetMass(std::string Nuclide, std::string Units);
  Int_t GetN(std::string Nuclide);
  std::string GetName(Int_t Index);
  std::string GetName(Int_t Z, Int_t A);
  Int_t GetZ(Int_t Index);
  Int_t GetZ(std::string Nuclide);
  Int_t IsMeasured(Int_t Index);
  Int_t LoadAMEFile(std::string file);
  Int_t LoadNubaseFile(std::string file);

private:
  std::string GetProperName(std::string Nuclide);
  void Init();
  void SetGSProperties();

  Int_t Size;
  Int_t *A;
  Int_t *N;
  Int_t *Z;
  Int_t *Measured;
  std::string *Name;
  Double_t *Mass; // micro-u
  Float_t *GSspin;
  Int_t *GSparity;

  Int_t nubaseLines;
  // nubase2020.asc row.
  struct nubaseEntry {
    Int_t massNumber;
    Int_t atomicNumber;
    Int_t stateIndex; // 0 = ground state; 1..3 isomers
    std::string Aelem;
    char Ssym;
    Double_t massExcess; // keV
    Int_t massMeas;      // 1 = measured, 0 = systematics (#)
    Double_t massUnc;
    Int_t massUncMeas;
    Float_t exc; // isomer excitation energy [keV]
    Int_t excMeas;
    Float_t excUnc;
    Int_t excUncMeas;
    Int_t uncertIGorder; // 1 = order of isomer or g.s. uncertain
    Float_t halfLife;
    Int_t halfLifeMeas;
    Int_t stable;
    std::string unitTime;
    Float_t halfLifeUnc;
    Int_t halfLifeUncMeas;
    Float_t twoJ;  // 2J
    Int_t spinUnc; // 0 = uncertain
    Int_t parity;  // ±1
    Int_t parityUnc;
    char JPiExtra; // '*' measured, '#' systematics, 'T' isospin
    Int_t yearENSDF;
    Int_t discovery;
    std::string decayBR;
    Double_t mass; // micro-u
    Double_t timeToSec;
    Double_t halfLifeSec;
  };
};

#endif

/***********************************************************************
Header file: SRIM_Table_Maker.hpp

Description: A class that creates input files for SRIM and then runs 
             SRIM (actually, SRModule.exe) via the program 'wine'.  The 
             SRIM outputs are limited to the following predefined cases: 
             0 - CD2
             1 - CF4 (gas)
             2 - Helium-3 (gas)
             3 - Helium-4 (gas)
             4 - Kapton
             5 - Silicon
             6 - LiF
             7 - Havar
             8 - Oxygen-16 (gas, diatomic molecule)
             9 - Oxygen-18 (gas, diatomic molecule)
            10 - Deuterium (gas, diatomic molecule)
            11 - Titanium-48
            12 - Gold-197
            13 - Mylar
            14 - Si3N4
            15 - CH2
            16 - 6Li
            17 - CH4 (gas)

Author: Daniel Santiago-Gonzalez
2014-07

***********************************************************************/
#ifndef SRIM_Table_Maker_hpp_INCLUDED   
#define SRIM_Table_Maker_hpp_INCLUDED   

#include <iostream>
#include <fstream>
#include <string.h>

// To excecute an external program.
#include <unistd.h>
#include <sys/wait.h>
// Contains 'chdir' function.
#include <unistd.h>
// Contains 'remove' function.
#include <cstdio>

class SRIM_Table_Maker{
public:

  SRIM_Table_Maker(std::string SRModulePath);
  int GetMediumIndex(std::string Medium);
  int GetMaxCases();
  double GetDensity(std::string Medium);
  void MakeTable(std::string SRIM_output, int Case, int Charge/*e*/, double Mass/*u*/);
  void MakeTable(std::string SRIM_output, std::string TgtMed, int Charge/*e*/, double Mass/*u*/);
  void SetGasDensity(int Case, float Pressure /*Torr*/, float Temperature /*K*/);
  void SetGasDensity(std::string TgtMed, float Pressure /*Torr*/, float Temperature /*K*/);

private:

  void MakeSRIMInputFile(std::string SRIM_output, int Case, int Charge, double Mass /*u*/);
  int RunSRModule();

  std::string SRModulePath;

  int MaxCases, MaxNumElem;
  int* Phase;
  double* Density;
  std::string* CompCorr;
  int* NumElem;
  std::string** TgtComp;
  std::string* TgtMed;

};

#endif

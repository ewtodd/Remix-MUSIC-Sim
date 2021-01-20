///////////////////////////////////////////////////////////////////////////////////////
// Methods for SRIM_Table_Maker class
// Author: Daniel Santiago-Gonzalez
// 2014-07
//
// See header file for class description.
// Compile with: g++ -shared -fPIC SRIM_Table_Maker.cpp -o SRIM_Table_Maker.so
///////////////////////////////////////////////////////////////////////////////////////

#include "SRIM_Table_Maker.hpp"

using namespace std;


///////////////////////////////////////////////////////////////////////////////////////
// Constructor.  It requires the path of the SRModule.exe program given as a string.
///////////////////////////////////////////////////////////////////////////////////////
SRIM_Table_Maker::SRIM_Table_Maker(string SRModulePath)
{
  this->SRModulePath = SRModulePath;
  MaxCases = 19;              // <-- Increase this number if you add a new case
  MaxNumElem = 10;
  Phase = new int[MaxCases];
  Density = new double[MaxCases];
  CompCorr = new string[MaxCases];
  NumElem = new int[MaxCases];
  TgtComp = new string*[MaxCases];
  TgtMed = new string[MaxCases];
  for (int c=0; c<MaxCases; c++)
    TgtComp[c] = new string[MaxNumElem];

  // Case 0: CD2
  // Density ref: http://homepages.rpi.edu/~danony/Papers/2010/Deuterated%20Target.pdf
  TgtMed[0] = "CD2";
  Phase[0] = 0;
  Density[0] = 1.06;
  CompCorr[0] = ".9843957";
  NumElem[0] = 2;
  TgtComp[0][0] = "1   \"Hydrogen\"               4             2.014";
  TgtComp[0][1] = "6   \"Carbon\"                 2             12.011";

  // Case 1: CF4 (gas)
  TgtMed[1] = "CF4";
  Phase[1] = 1;
  Density[1] = 0.00125;  // Default SRIM density. Change with SetGasDensity(P, T).
  CompCorr[1] = ".9585996";
  NumElem[1] = 2;
  TgtComp[1][0] = "6   \"Carbon\"                 1             12.011";
  TgtComp[1][1] = "9   \"Fluorine\"               4             18.998";

  // Case 2: Helium-3 (gas)
  TgtMed[2] = "3He";
  Phase[2] = 1;
  Density[2] = 0.00012545; // LISE++ density at 273 K and 760 Torr. Change with SetGasDensity(P, T).
  CompCorr[2] = "1";
  NumElem[2] = 1;
  TgtComp[2][0] = "2   \"Helium\"                 1             3.016";

  // Case 3: Helium-4 (gas)
  TgtMed[3] = "4He";
  Phase[3] = 1;
  Density[3] = 0.0001665;  // LISE++ density at 273 K and 760 Torr. Change with SetGasDensity(P, T).
  CompCorr[3] = "1";
  NumElem[3] = 1;
  TgtComp[3][0] = "2   \"Helium\"                 1             4.003";

  // Case 4: Kapton
  TgtMed[4] = "Kapton";
  Phase[4] = 0;
  Density[4] = 1.42;
  CompCorr[4] = "1";
  NumElem[4] = 4;
  TgtComp[4][0] = "1   \"Hydrogen\"               25.64101      1.008";
  TgtComp[4][1] = "6   \"Carbon\"                 56.41008      12.011";
  TgtComp[4][2] = "7   \"Nitrogen\"               5.128087      14.007";
  TgtComp[4][3] = "8   \"Oxygen\"                 12.82082      15.999";

  // Case 5: Silicon
  TgtMed[5] = "Si";
  Phase[5] = 0;
  Density[5] = 2.3212;
  CompCorr[5] = "1";
  NumElem[5] = 1;
  TgtComp[5][0] = "14   \"Silicon\"               1             28.086";
  
  // Case 6: Lithium6-Fluoride
  // Refs: http://www.sigmaaldrich.com/catalog/product/aldrich/601411?lang=en&region=US
  TgtMed[6] = "LiF";
  Phase[6] = 0;
  Density[6] = 2.55;
  CompCorr[6] = "1";
  NumElem[6] = 2;
  TgtComp[6][0] = "3   \"Lithium\"                1             6.015";
  TgtComp[6][1] = "9   \"Fluorine\"               1             18.998";

  // Case 7: Havar
  TgtMed[7] = "Havar";
  Phase[7] = 0;
  Density[7] = 8.3;
  CompCorr[7] = "1";
  NumElem[7] = 8;
  TgtComp[7][0] = "6   \"Carbon\"                 .9647583      12.011";
  TgtComp[7][1] = "24   \"Chromium\"              22.28578      51.996";
  TgtComp[7][2] = "25   \"Manganese\"             1.687388      54.938";
  TgtComp[7][3] = "26   \"Iron\"                  18.11391      55.847";
  TgtComp[7][4] = "27   \"Cobalt\"                41.78285      58.933";
  TgtComp[7][5] = "28   \"Nickel\"                12.83355      58.69";
  TgtComp[7][6] = "42   \"Molybdenum\"            1.44937       95.94";
  TgtComp[7][7] = "74   \"Tungsten\"              .8823931      183.85";

  // Case 8: 16O (gas, diatomic molecule)
  TgtMed[8] = "16O";
  Phase[8] = 1;
  Density[8] = 0.00066514; // LISE++ density at 273 K and 760 Torr. Change with SetGasDensity(P, T).
  CompCorr[8] = "1";
  NumElem[8] = 1;
  TgtComp[8][0] = "8   \"Oxygen\"                 2            15.999";
  
  // Case 9: 18O (gas, diatomic molecule)
  TgtMed[9] = "18O";
  Phase[9] = 1;
  Density[9] = 0.00066514; // Same as 16O case.
  CompCorr[9] = "1";
  NumElem[9] = 1;
  TgtComp[9][0] = "8   \"Oxygen\"                 2            17.999";

  // Case 10: D2 (gas, diatomic molecule)
  TgtMed[10] = "D2";
  Phase[10] = 1;
  Density[10] = 0.0001; 
  CompCorr[10] = "1";
  NumElem[10] = 1;
  TgtComp[10][0] = "1   \"Hydrogen\"                 2            2.014";

  // Case 11: Titanium-48
  TgtMed[11] = "48Ti";
  Phase[11] = 0;
  Density[11] = 4.5189; 
  CompCorr[11] = "1";
  NumElem[11] = 1;
  TgtComp[11][0] = "22   \"Titanium\"              1             47.9";

  // Case 12: Gold-197
  TgtMed[12] = "197Au";
  Phase[12] = 0;
  Density[12] = 19.311; 
  CompCorr[12] = "1";
  NumElem[12] = 1;
  TgtComp[12][0] = "79   \"Gold\"                  1             196.97";

  // Case 13: Mylar
  TgtMed[13] = "Mylar";
  Phase[13] = 0;
  Density[13] = 1.397; 
  CompCorr[13] = "0.9570329";
  NumElem[13] = 3;
  TgtComp[13][0] = "1   \"Hydrogen\"                8             1.008";
  TgtComp[13][1] = "6   \"Carbon\"                  10            12.011";
  TgtComp[13][2] = "8   \"Oxygen\"                  4             15.999";

  // Case 14: Si3N4
  TgtMed[14] = "Si3N4";
  Phase[14] = 0;
  Density[14] = 3.44; 
  CompCorr[14] = "1";
  NumElem[14] = 2;  
  TgtComp[14][0] = "14   \"Silicon\"               4             28.086";
  TgtComp[14][1] = "7   \"Nitrogen\"               3             14.007";

  // Case 15: CH2
  TgtMed[15] = "CH2";
  Phase[15] = 0;
  Density[15] = 0.93;
  CompCorr[15] = ".9843957";
  NumElem[15] = 2;
  TgtComp[15][0] = "1   \"Hydrogen\"               4             1.008";
  TgtComp[15][1] = "6   \"Carbon\"                 2             12.011";

  // Case 16: 6Li
  TgtMed[16] = "6Li";
  Phase[16] = 0;
  Density[16] = 0.534;
  CompCorr[16] = "1";
  NumElem[16] = 1;
  TgtComp[16][0] = "3   \"Lithium\"                1             6.015";

  // Case 17: CF4 (gas)
  TgtMed[17] = "CH4";
  Phase[17] = 1;
  Density[17] = 0.000717;  // Default SRIM density. Change with SetGasDensity(P, T).
  CompCorr[17] = "1.004811";
  NumElem[17] = 2;
  TgtComp[17][0] = "1   \"Hydrogen\"               4             1.008";
  TgtComp[17][1] = "6   \"Carbon\"                 1             12.011";

  // Case 18: 40Ar (gas)
  TgtMed[18] = "40Ar";
  Phase[18] = 1;
  Density[18] = 1.6504;    // Default SRIM density. Change with SetGasDensity(P, T).
  CompCorr[18] = "1";
  NumElem[18] = 1;
  TgtComp[18][0] = "18   \"Argon\"                 1             39.9481";

}


///////////////////////////////////////////////////////////////////////////////////////
// Returns the density of the medium specified as a string in the argument. Returns 1
// if the medium is not found.
///////////////////////////////////////////////////////////////////////////////////////
double SRIM_Table_Maker::GetDensity(string Medium)
{
  double density = 1;
  for (int i=0; i<MaxCases; i++)
    if (Medium==TgtMed[i]) {
      density = Density[i];
      break;
    }
  return density;
}


///////////////////////////////////////////////////////////////////////////////////////
// Returns the index of the medium specified as a string in the argument. Returns -1
// if the medium is not found.
///////////////////////////////////////////////////////////////////////////////////////
int SRIM_Table_Maker::GetMediumIndex(string Medium)
{
  int MedID = -1;
  for (int i=0; i<MaxCases; i++)
    if (Medium==TgtMed[i]) {
      MedID = i;
      break;
    }
  return MedID;
}


///////////////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////////////
int SRIM_Table_Maker::GetMaxCases()
{
  return MaxCases;
}


///////////////////////////////////////////////////////////////////////////////////////
// Function that writes the input file (SR.IN) for SRModule.exe.
// Since SRModule.exe is a program for windows the input text file must have the DOS 
// end-of-line signal composed by two chars (0D0A, in hex).  Simple endl or \n do not
// work, hence the two char EOL1,2.
///////////////////////////////////////////////////////////////////////////////////////
void SRIM_Table_Maker::MakeSRIMInputFile(string SRIM_output, int Case, int Charge, double Mass /*u*/)
{
  // These 2 lines get the current working directory.
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));

  // If the SRIM_output string does not contain the root directory (/)
  // then we need to add the cwd. Otherwise SRModule will try to
  // create the SRIM table in a directory that may not exists and will
  // give an error message.
  if (SRIM_output[0]!='/') {
    string temp = cwd;
    temp.append("/");
    temp.append(SRIM_output);
    SRIM_output = temp;
  }
  
  // Energy range in keV
  float EMin = 0.5, EMax = 100000.0*Mass;
  char EOL1 = 0x0D;
  char EOL2 = 0x0A;
  string InputFile = SRModulePath + "SR.IN";
  ofstream SR_IN(InputFile.c_str());
  SR_IN << "---Stopping/Range Input Data (Number-format: Period = Decimal Point)" << EOL1 << EOL2
	<< "---Output File Name" << EOL1 << EOL2
	<< "\"" << SRIM_output << "\"" << EOL1 << EOL2
	<< "---Ion(Z), Ion Mass(u)" << EOL1 << EOL2
	<< Charge << "   " << Mass << EOL1 << EOL2
	<< "---Target Data: (Solid=0,Gas=1), Density(g/cm3), Compound Corr." << EOL1 << EOL2
	<< Phase[Case] << "    " << Density[Case] << "    " << CompCorr[Case] << EOL1 << EOL2
	<< "---Number of Target Elements" << EOL1 << EOL2
	<< " " << NumElem[Case] << EOL1 << EOL2
	<< "---Target Elements: (Z), Target name, Stoich, Target Mass(u)" << EOL1 << EOL2;
  for (int elem=0; elem<NumElem[Case]; elem++) 
    SR_IN << TgtComp[Case][elem] << EOL1 << EOL2;
  SR_IN << "---Output Stopping Units (1-8)" << EOL1 << EOL2
	<< " 3" << EOL1 << EOL2
	<< "---Ion Energy : E-Min(keV), E-Max(keV)" << EOL1 << EOL2
	<< " " << EMin << "    " << EMax << EOL1 << EOL2;
  SR_IN.close();
  return;
}


///////////////////////////////////////////////////////////////////////////////////////
// This method makes calls the function that creates the SRIM input file and then runs
// the SRModule program. It uses the Case number to determine which table to make.
///////////////////////////////////////////////////////////////////////////////////////
void SRIM_Table_Maker::MakeTable(string SRIM_output, int Case, int Charge, double Mass /*u*/)
{
  if (Charge>0 && Mass>0) {
    MakeSRIMInputFile(SRIM_output, Case, Charge, Mass);
    RunSRModule();
  }
  return;
}


///////////////////////////////////////////////////////////////////////////////////////
// This method makes calls the function that creates the SRIM input file and then runs
// the SRModule program. It uses the target medium name (TgtMed) to determine which
// table to make.
///////////////////////////////////////////////////////////////////////////////////////
void SRIM_Table_Maker::MakeTable(string SRIM_output, string TgtMed, int Charge, double Mass /*u*/)
{
  int Case = -1;
  for (int c=0; c<MaxCases; c++) 
    if (this->TgtMed[c]==TgtMed) {
      Case = c;
      break;
    }
  if (Case>=0) 
    MakeTable(SRIM_output, Case, Charge, Mass);
  else
    cout << "SRIM_Table_Maker warning: " << TgtMed << " not a valid medium" << endl;
  return;
}


///////////////////////////////////////////////////////////////////////////////////////
// References: 
// 1) http://www.cplusplus.com/forum/lounge/17684/
// 2) http://stackoverflow.com/questions/15750068/fork-and-execl-in-c-execl-is-failing
// 3) http://www.srim.org/SRIM/Tutorials/SRIM%20Special%20-%20SR-Module.doc
// 4) http://linux.die.net/man/3/chdir
// 5) https://www.systutorials.com/docs/linux/man/3p-execlp/#lbAH
///////////////////////////////////////////////////////////////////////////////////////
int SRIM_Table_Maker::RunSRModule()
{
  string TempDir;
  pid_t pid = fork(); // Create a child process
  switch (pid) {
  case -1: // Error 
    std::cerr << "ERROR: fork() failed.\n";
    _exit(1);

  case 0: // Child process

    // SRModule requires the system to change dir.
    if (chdir(SRModulePath.c_str())==-1)      
      cout << "SRIM_Table_Maker error: directory \'" << SRModulePath << "\' does not exist" << endl;

    // Execute the program
    execlp("/usr/bin/wine", "wine", "SRModule.exe", NULL); 
    std::cerr << "ERROR: execl() failed! Check path of SRModule.exe: \n"
	      << SRModulePath << "\n";
    // execl doesn't return unless there's an error
    _exit(1);
    
  default: // Parent process 
    // cout << "\nSRModule created with process id " << pid << "\n";
    int* status = 0;
    while (wait(status) != pid) {
      // Wait until child process exits.
    }
    //cout << "Process exited with " << WEXITSTATUS(status) << "\n";
    return 0;
  }
}


///////////////////////////////////////////////////////////////////////////////////////
// Calculates the gas density for selected cases.  The density can be obtained from
// the ideal gas law 
//   PV = (m/M)RT,
// where m (in g) is the total mass of the gas and M (in g/mol) is the molar mass of a
// gas molecule. The density, rho=m/V (in g/cm^3), is then
//   rho = PM/RT
// We now convert the pressure units in the gas constant, R, from atm to Torr
//   R = 82.05736 cm^3*atm/(K*mol) = 62363.59 cm^3*Torr/(K*mol)
///////////////////////////////////////////////////////////////////////////////////////
void SRIM_Table_Maker::SetGasDensity(int Case, float Pressure /*Torr*/, float Temperature /*K*/)
{
  double MolarMass = 0; // g/mol
  if (Case==1 || Case==2 || Case==3 || Case==8 || Case==9 || Case==10 || Case==17 || Case==18) {
    if (Case==1)
      MolarMass = 88.0043;
    else if (Case==2)
      MolarMass = 3.016;
    else if (Case==3)
      MolarMass = 4.003;
    else if (Case==8)
      MolarMass = 31.999;  // http://encyclopedia.airliquide.com/encyclopedia.asp?GasID=48#GeneralData
    else if (Case==9)
      MolarMass = 35.998;  // Just from twice the mass of 18O
    else if (Case==10)
      MolarMass = 4.0282;  // http://encyclopedia.airliquide.com/Encyclopedia.asp?GasID=20#GeneralData
    else if (Case==17)
      MolarMass = 16.04;   // https://en.wikipedia.org/wiki/Methane
    else if (Case==18)
      MolarMass = 39.948;   

    Density[Case] = Pressure*MolarMass/Temperature/62363.59;  // g/cm^3
  }
  else 
    cout << "SRIM_Table_Maker warning: Case " << Case << " is not a gas." << endl;
  return;
}


///////////////////////////////////////////////////////////////////////////////////////
// Same as above but the user specifies the name of the target medium and not the 
// Case number.
///////////////////////////////////////////////////////////////////////////////////////
void SRIM_Table_Maker::SetGasDensity(string TgtMed, float Pressure /*Torr*/, float Temperature /*K*/)
{
  int Case = -1;
  for (int c=0; c<MaxCases; c++) 
    if (this->TgtMed[c]==TgtMed) {
      Case = c;
      break;
    }
  if (Case>=0) 
    SetGasDensity(Case, Pressure, Temperature);
  else
    cout << "SRIM_Table_Maker warning: " << TgtMed << " not a valid medium" << endl;
  return;
}

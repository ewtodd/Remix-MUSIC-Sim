

#include <iostream>
#include <cstdio>
//#include <unistd.h>  // needed for getcwd()
#include "NuclideFinder.hpp"
#include "SRIM_Table_Maker.hpp"

using namespace std;


struct controlFileParams {
  string SRModuleDir;
  string OutputDir;
  float GasP; // Torr
  float GasT; // K
  int GasIndex;
  int NumParticles;
  const int MaxNumPart = 10;
  string* particle = new string[MaxNumPart];
};


int loadCtrlFile(char* fileName, controlFileParams& ctr);


int main(int argc, char* argv[])
{
  cout << "==========================================================================" << endl;
  cout << "|--- make-srim-table ----------------------------------------------------|" << endl;
  cout << "| Created by Daniel Santiago-Gonzalez (Argonne National Laboratory)      |" << endl;
  cout << "| ver 1.0 (2020/1)                                                       |" << endl;
  cout << "| Usage: ./make-srim-table control.file                                  |" << endl;
  cout << "| For documentation, software updates and license details, please visit: |" << endl;
  cout << "| https://gitlab.phy.anl.gov/music/sim                                   |" << endl;
  cout << "==========================================================================" << endl;

  controlFileParams ctr;
  
  if (argc==2) {
    cout << "make-srim-table# Loading control file: " << argv[1] << endl;
    if (loadCtrlFile(argv[1], ctr)==0)
      cout << "make-srim-table warning: invalid control file (check address)." << endl;
  }
  else {
    cout << "make-srim-table error: only one argument is expected." << endl;
    return 0;    
  }
  
  char* srimFile;
  
  auto NuF = new NuclideFinder();
  auto SRIM = new SRIM_Table_Maker(ctr.SRModuleDir);
  SRIM->SetGasDensity(ctr.GasIndex, ctr.GasP, ctr.GasT);
  string SRIMFile[ctr.NumParticles];
  for (int p=0; p<ctr.NumParticles; p++) {
    cout << "make-srim-table# Creating SRIM files " << ctr.particle[p] << " ..." << endl;
    sprintf(srimFile,
	    "%s%s_in_GasIndex%i_%.0fTorr_%.0fK.srim",
	    ctr.OutputDir.c_str(), ctr.particle[p].c_str(), ctr.GasIndex, ctr.GasP, ctr.GasT);
    cout << "make-srim-table# " << srimFile << endl;
    // Before making the tables we need the masses (in u) and atomic numbers (in e).
    double m = NuF->GetMass(ctr.particle[p], "u");
    int Z = NuF->GetZ(ctr.particle[p]);
    // Now we have all the information to make the energy loss tables.
    SRIM->MakeTable(srimFile, ctr.GasIndex, Z, m);
  }
  return 1;
}


////////////////////////////////////////////////////////////////////////////////////
// Method to load the parameters from the control file
// Control file example below horizontal line (---), first two lines are skipped
// and can be used for comments. Third column can be used for comments:
// ------------------------------------------------------
// Control file for make-srim-table                     | File description
// Parameter      Value     Comment                     | Column description
// DirSRModule    SRModule/ Directory of SRModule       | first data line
// ...                                                  | 
////////////////////////////////////////////////////////////////////////////////////
int loadCtrlFile(char* fileName, controlFileParams& ctr)
{
  int Status = 0;
  ifstream Ctrl(fileName);
  string aux, ParName, ParVal;
  if (Ctrl.is_open()) {
    getline(Ctrl,aux);  // skipping first line
    getline(Ctrl,aux);  // skipping second line
    while (!Ctrl.eof()) {
      Ctrl >> ParName >> ParVal;
      getline(Ctrl,aux);

      if (ParName=="SRModuleDir")
	ctr.SRModuleDir = ParVal;
      else if (ParName=="OutputDir")
	ctr.OutputDir = ParVal;

      else if (ParName=="GasIndex")
	ctr.GasIndex = atoi(ParVal.c_str());
      else if (ParName=="GasP")
	ctr.GasP = atof(ParVal.c_str());
      else if (ParName=="GasT")
	ctr.GasT = atof(ParVal.c_str());

      // Number of particles for which the tables will be generated for the given medium
      else if (ParName=="NumParticles")
      	ctr.NumParticles = atoi(ParVal.c_str());
      else if (ParName=="p0")
	ctr.particle[0] = ParVal;
      else if (ParName=="p1")
	ctr.particle[1] = ParVal;
      else if (ParName=="p2")
	ctr.particle[2] = ParVal;
      else if (ParName=="p3")
	ctr.particle[3] = ParVal;
      else if (ParName=="p4")
	ctr.particle[4] = ParVal;
      else if (ParName=="p5")
	ctr.particle[5] = ParVal;
      else if (ParName=="p6")
	ctr.particle[6] = ParVal;
      else if (ParName=="p7")
	ctr.particle[7] = ParVal;
      else if (ParName=="p8")
	ctr.particle[8] = ParVal;
      else if (ParName=="p9")
	ctr.particle[9] = ParVal;
      else
	cout << "make-srim-table warning: control file parameter \'" << ParName << "\' not recognized."
	     << endl;
      
      
#if 0
      else if (ParName=="str")
      	ctr.var = ParVal;
      else if (ParName=="int")
      	ctr.var = atoi(ParVal.c_str());
      else if (ParName=="float")
      	ctr.var = atof(ParVal.c_str());
#endif
    }
    Status = 1;
  } 
  return Status;
}


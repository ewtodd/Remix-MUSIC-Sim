/////////////////////////////////////////////////////////////////////////////
// Description: ROOT macro that creates SRIM tables for the MUSIC Simulator.
//
// Usage: root -l MakeSRIMTables.C
//
// Created by: Daniel Santiago-Gonzalez
// Date: Dec 2018
/////////////////////////////////////////////////////////////////////////////
#include <unistd.h>  // needed for getcwd()

void MakeSRIMTables()
{
  // These 3 lines get the current working directory (to load the SRIM files)
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  string SRIMdir = cwd;
  SRIMdir += "/SRIM_files/";

  //=======================================================================
  // CONTROL PANEL
  //
  string OutputDir = SRIMdir;
  string SRModPath = "/home/dasago/.wine/drive_c/Program\ Files\ \(x86\)/SRIM/SR\ Module/";
  //
  // Particles for which the Stopping Power tables will be generated
  const int NumParticles = 7;
  string particle[NumParticles];
  particle[0] = "59Fe";
  particle[1] = "1H"; 
  particle[2] = "4He";
  particle[3] = "62Ni";
  particle[4] = "61Ni";
  particle[5] = "62Co";
  particle[6] = "59Co";
  //
  // Target gas related stuff
  float GasP = 600; // Torr
  float GasT = 293; // K
  int GasIndex = 3;
  // You have to use the same index numbers as in the SRIM_Table_Maker class.
  //   0 - CD2
  //   1 - CF4 (gas)
  //   2 - Helium-3 (gas)
  //   3 - Helium-4 (gas)
  //   4 - Kapton
  //   5 - Silicon
  //   6 - LiF
  //   7 - Havar
  //   8 - Oxygen-16 (gas)
  //   9 - Oxygen-18 (gas)
  //  10 - Deuterium (gas)
  //  17 - Methane (CH4 gas)
  //=======================================================================


  NuclideFinder* NuF = new NuclideFinder();
  SRIM_Table_Maker* SRIM;
  SRIM = new SRIM_Table_Maker(SRModPath);
  SRIM->SetGasDensity(GasIndex, GasP, GasT);
  string SRIMFile[NumParticles ];
  for (int p=0; p<NumParticles; p++) {
    cout << "Creating SRIM files for " << particle[p] << " ..." << endl;
    SRIMFile[p] = OutputDir + particle[p] + 
      Form("_in_GasIndex%i_%.0fTorr_%.0fK.srim", GasIndex, GasP, GasT);
    cout << SRIMFile[p] << endl;
    // Before making the tables we need the masses (in u) and atomic numbers (in e).
    double m = NuF->GetMass(particle[p], "u");
    int Z = NuF->GetZ(particle[p]);
    // Now we have all the information to make the energy loss tables.
    SRIM->MakeTable(SRIMFile[p], GasIndex, Z, m);
  }
}

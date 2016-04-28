/////////////////////////////////////////////////////////////////////////////
// Description: Script to tests the MUSIC_Simulator class.
//
// Usage: root -l RunMUSICSim.C
//
// Created by: Daniel Santiago-Gonzalez
// Date: March 2016
/////////////////////////////////////////////////////////////////////////////
{
  //=======================================================================
  // CONTROL PANEL
  string SRIM_dir = "/home/dasago/Dropbox/Codes/MUSIC/Simulator/SRIM_files/";

  // Geometrical parameters (all distances in cm)
  string AnodeGeom = "AnodeGeometry";

  string beam = "20Ne";
  string target = "4He";
  string compound = "24Mg";
  string light = "n";
  string heavy = "23Mg";
  // Energy of the beam after the window.
  double Kb = 55;
  
  int Strip = 3;
  int NEvents = 50;
  double MaxTime = 1000; // ns
  double UserDT = 0.1;     // ns

  // Target gas related stuff
  float GasP = 400; // Torr
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
  //=======================================================================


  /////////////////////////////////////////////////////////////////////////////
  // Load the necessary libraries for the script to run.
  /////////////////////////////////////////////////////////////////////////////
  gStyle->SetOptStat("");
  gSystem->Load("../../PhysicsTools/EnergyLoss.so"); 
  gSystem->Load("../../PhysicsTools/FourVector.so"); 
  gSystem->Load("../../PhysicsTools/Particle.so"); 
  gSystem->Load("../../PhysicsTools/SRIM_Table_Maker_cpp.so");
  gSystem->Load("../../NuclideFinder/NuclideFinder_cpp.so"); 
  gSystem->Load("MUSIC_Simulator_cpp.so"); 




  /////////////////////////////////////////////////////////////////////////////
  // In this part, energy loss tables are created using a SRIM_Table_Maker
  // object. Alternatively, the tables can be generated 'by hand'.
  /////////////////////////////////////////////////////////////////////////////
  string SRModPath = "/home/dasago/.wine/drive_c/Program\ Files\ \(x86\)/SRIM/SR\ Module/";
  SRIM_Table_Maker* SRIM = new SRIM_Table_Maker(SRModPath);
  SRIM->SetGasDensity(GasIndex, GasP, GasT);
  string particle[4];
  particle[0] = beam;
  particle[1] = compound;
  particle[2] = light;
  particle[3] = heavy;
  string SRIMFile[4];
  for (int p=0; p<4; p++) {
    cout << "Creating SRIM files for " << particle[p] << " ..." << endl;
    SRIMFile[p] = SRIM_dir + particle[p] + Form("_in_4He_%.0fTorr_%.0fK.srim", GasP, GasT);
  }
  // Before making the tables we need the masses (in u) and atomic numbers (in e).
  NuclideFinder* NuF = new NuclideFinder();
  double mb = NuF->GetMass(beam, "u");
  int Zb = NuF->GetZ(beam);
  double mt = NuF->GetMass(target, "u");
  int Zt = NuF->GetZ(target);
  double mc = NuF->GetMass(compound, "u");
  int Zc = NuF->GetZ(compound);
  double ml = NuF->GetMass(light, "u");
  int Zl = NuF->GetZ(light);
  double mh = NuF->GetMass(heavy, "u");
  int Zh = NuF->GetZ(heavy);
  // Now we have all the information to make the energy loss tables.
  SRIM->MakeTable(SRIMFile[0], GasIndex, Zb, mb);
  SRIM->MakeTable(SRIMFile[1], GasIndex, Zc, mc);
  SRIM->MakeTable(SRIMFile[2], GasIndex, Zl, ml);
  SRIM->MakeTable(SRIMFile[3], GasIndex, Zh, mh);



  /////////////////////////////////////////////////////////////////////////////
  // Simulator stuff
  /////////////////////////////////////////////////////////////////////////////
  // gRandom->SetSeed();
  MUSIC_Simulator* MUSIC = new MUSIC_Simulator();
  MUSIC->SetStripEnergyResolution(0.025);
  // Geometry
  MUSIC->SetAnode(AnodeGeom, 90);
  // Beam
  MUSIC->SetBeamParticle(beam, kBlack, SRIMFile[0], Kb);
  // Target
  MUSIC->SetTargetParticle(target);
  // Compound particle
  MUSIC->SetCompoundParticle(compound);
  // Light evaporation residue (e.g. proton)
  MUSIC->SetLightParticle(light, kRed, SRIMFile[2]);
  // Heavy evaporation residue (e.g. 23Na)
  MUSIC->SetHeavyParticle(heavy, kBlue, SRIMFile[3]);

  // Release the Kracken!!
  MUSIC->Simulate(Strip, NEvents, MaxTime, UserDT, 1);
  MUSIC->WriteTraces(Form("Traces_Stp%d_%s_%s.root",Strip,target.c_str(),light.c_str()));

}

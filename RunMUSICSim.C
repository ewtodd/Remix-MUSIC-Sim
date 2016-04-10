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
  string ParamDir = "/home/dasago/Dropbox/Codes/MUSIC/Simulator/SRIM_files/";

  // Geometrical parameters (all distances in cm)
  string AnodeGeom = "AnodeGeometry";

  string beam = "20Ne";
  string target = "4He";
  string fused = "24Mg";
  string light = "4He";
  string heavy = "20Ne";
  // Energy of the beam after the window.
  double Kb = 55;
  
  int Strip = 4;
  int NEvents = 100;

  // Target gas related stuff
  float GasP = 350; // Torr
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


  NuclideFinder* NuF = new NuclideFinder();
  // Masses (must be in MeV/c^2) and charges (in e)
  double mb = NuF->GetMass(beam, "MeV/c^2");
  int Zb = NuF->GetZ(beam);
  double mt = NuF->GetMass(target, "MeV/c^2");
  int Zt = NuF->GetZ(target);
  double mf = NuF->GetMass(fused, "MeV/c^2");
  int Zf = NuF->GetZ(fused);
  double ml = NuF->GetMass(light, "MeV/c^2");
  int Zl = NuF->GetZ(light);
  double mh = NuF->GetMass(heavy, "MeV/c^2");
  int Zh = NuF->GetZ(heavy);

  /////////////////////////////////////////////////////////////////////////////
  // In this part, energy loss tables are created using a SRIM_Table_Maker
  // object. Alternatively, the tables can be generated 'by hand'.
  /////////////////////////////////////////////////////////////////////////////
  string SRModPath = "/home/dasago/.wine/drive_c/Program\ Files\ \(x86\)/SRIM/SR\ Module/";
  SRIM_Table_Maker* SRIM = new SRIM_Table_Maker(SRModPath);
  SRIM->SetGasDensity(GasIndex, GasP, GasT);
  string particle[4];
  particle[0] = beam;
  particle[1] = fused;
  particle[2] = light;
  particle[3] = heavy;
  string SRIMFile[4];
  for (int p=0; p<4; p++) {
    cout << "Creating SRIM files for " << particle[p] << " ..." << endl;
    SRIMFile[p] = particle[p] + Form("_in_4He_%.0fTorr_%.0fK.srim", GasP, GasT);
  }
  double Conv = 931.494061; // from u to MeV/c^2
  SRIM->MakeTable(ParamDir+SRIMFile[0], GasIndex, Zb, mb/Conv);
  SRIM->MakeTable(ParamDir+SRIMFile[1], GasIndex, Zf, mf/Conv);
  SRIM->MakeTable(ParamDir+SRIMFile[2], GasIndex, Zl, ml/Conv);
  SRIM->MakeTable(ParamDir+SRIMFile[3], GasIndex, Zh, mh/Conv);



  /////////////////////////////////////////////////////////////////////////////
  // Simulator stuff
  /////////////////////////////////////////////////////////////////////////////
  gRandom->SetSeed();
  MUSIC_Simulator* MUSIC = new MUSIC_Simulator();
  MUSIC->SetStripEnergyResolution(0.05);
  MUSIC->SetParamDirectory(ParamDir);
  // Geometry
  MUSIC->SetAnode(AnodeGeom, 90);
  // Beam
  MUSIC->SetBeamParticle("beam", mb, Zb, kBlack, Kb);
  MUSIC->SetEnergyLossFile("beam", SRIMFile[0]);
  // Target
  MUSIC->SetTargetParticle("target", mt, Zt);
  // Fused particle
  MUSIC->SetFusedParticle(fused, mf, Zf);
  MUSIC->SetEnergyLossFile(fused, SRIMFile[1]);
  // Light outgoing particle (e.g. proton)
  MUSIC->SetLightParticle(light, ml, Zl, kRed);
  MUSIC->SetEnergyLossFile(light, SRIMFile[2]);
  // Heavy outgoing particle (e.g. 23Na)
  MUSIC->SetHeavyParticle(heavy, mh, Zh, kBlue);
  MUSIC->SetEnergyLossFile(heavy, SRIMFile[3]);
  int Reaction = 0;
  MUSIC->Simulate(Reaction, Strip, NEvents);
  MUSIC->WriteTraces(Form("TracesR%d_Stp%d_%s_%s.root",Reaction,Strip,target.c_str(),light.c_str()));

  TEveManager* Eve = new TEveManager(1000, 1000, kTRUE, "V");
  // Axes
  TEveArrow* Xaxis = new TEveArrow(20,0,0,-10,0,0);
  Xaxis->SetMainColor(kGreen); 
  Eve->AddElement(Xaxis);
  TEveArrow* Yaxis = new TEveArrow(0,20,0,0,-10,0);
  Yaxis->SetMainColor(kYellow);
  Eve->AddElement(Yaxis);
  TEveArrow* Zaxis = new TEveArrow(0,0,30,0,0,0);
  Eve->AddElement(Zaxis);
  Zaxis->SetMainTransparency(65);

  MUSIC->DrawTrajecotries(Eve);
  MUSIC->DrawMUSIC(Eve, 85);
}

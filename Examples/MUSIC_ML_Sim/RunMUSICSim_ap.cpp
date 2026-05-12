/////////////////////////////////////////////////////////////////////////////
// Description: ROOT macro to run the MUSIC_Simulator.
//
// Usage: root -l RunMUSICSim_ap.C
//
// Created by: Daniel Santiago-Gonzalez
// Date: Oct 2019
/////////////////////////////////////////////////////////////////////////////
#include <unistd.h>  // needed for getcwd()

//void RunMUSICSim_ap()
{
  // These 3 lines get the current working directory (to load the SRIM files)
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  string SRIMdir = cwd;

  //=======================================================================
  // CONTROL PANEL (units: dist=cm, energy=MeV, angle=deg)
  string AnodeGeom = "AnodeGeometry";
  int Pressure = 450; // Torr
  int ELossBins = 300;
  float MaxELoss = 5;
  SRIMdir += "/SRIM_files/";
  // beam
  string beam = "17F";
  string SRIMbeam = SRIMdir + Form("17F_in_GasIndex3_%dTorr_293K.srim", Pressure);
  // target (no SRIM file, assuming it does not propagate)
  string target = "4He";
  // compound nucleus (no SRIM file, assuming it does not propagate)
  string compound = "21Na";
  
  // Reac1: 17F(4He,p)20Ne
  // Chan1: 21Na*->p+20Ne*
  const int NumEvapPart = 1;
  string res[NumEvapPart];
  string SRIMres[NumEvapPart];
  string evap[NumEvapPart];
  string SRIMevap[NumEvapPart];
  evap[0] = "p";    SRIMevap[0] = SRIMdir + Form("1H_in_GasIndex3_%dTorr_293K.srim",Pressure);
  res[0] = "20Ne";   SRIMres[0] = SRIMdir + Form("20Ne_in_GasIndex3_%dTorr_293K.srim",Pressure);

  // Cosmetics
  int color[NumEvapPart];
  for (int er=0; er<NumEvapPart; er++) {
    if (evap[er]=="n")
      color[er] = kBlue;
    if (evap[er]=="p")
      color[er] = kRed;
  }

  double Kb = 37.5;//42.59;   // MeV - Energy of the beam after the Ti window and Al degrader
  int strip = 4;     // Strip where reaction takes place
  float Eres = 0.1;  // MeV - Strip energy resolution (larger values increase signal randomness)
  int NEvents = 100;   // Number of simulated events
  int Wait = 1;       // 1 - canvas waits for user's double click, 0 - no wait
  int Update = 1;     // 1 - update visuals for every event, 0 - don't
  double MaxTime = 2000;   // ns - max time for an event
  double SimStep = 0.001;     // cm - simulation steps size
  int Method = 0;    // Select the simulation method: 0 - Simulate, 1 - GenerateTraceDatabase
  string FileName = Form("Traces_17F_%.1fMeV_%dTorr_alpha_p_Noise100keV_Stp%d_test_pts.root", Kb, Pressure, strip);
  string FileOpt = "recreate"; // recreate or update

  // The following control variables only apply for GenerateTraceDatabase (Method=1)
  double ThCMMin = 0.1;
  double ThCMMax = 179.9;
  int ThSteps = 6;
  double PhiCMMin = 0.1;
  double PhiCMMax = 359.9;
  int PhiSteps = 6;
  //=======================================================================


  /////////////////////////////////////////////////////////////////////////////
  // Simulator stuff (typically not modified)
  /////////////////////////////////////////////////////////////////////////////

  MUSIC_Simulator* MUSIC = new MUSIC_Simulator();
  MUSIC->SetROOTSystemPointer(gSystem);
  MUSIC->SetPrintLevel(1);
  MUSIC->SetStripEnergyResolution(Eres);
  // Geometry
  MUSIC->SetAnode(AnodeGeom, 90/*transparency 0-100*/, ELossBins, MaxELoss);
  // Beam
  MUSIC->SetBeamParticle(beam, kBlack, SRIMbeam, Kb);
  // Target
  MUSIC->SetTargetParticle(target);
  // Compound particle
  MUSIC->SetCompoundParticle(compound);
  // Evaporation residues and particles
  for (int i=0; i<NumEvapPart; i++)
    MUSIC->SetEvapResAndPart(res[i], SRIMres[i], kGreen, evap[i], SRIMevap[i], color[i]);

  if (Method==0) {
    // Simulate events for one strip or generate trace data base (see below)
    MUSIC->Simulate(strip, NEvents, MaxTime, SimStep, Update, Wait, FileName, FileOpt);
  }
  else if (Method==1) {
    MUSIC->GenerateTraceDatabase("TraceDB_ap.root", 
				 ThCMMin, ThCMMax, ThSteps, 
				 PhiCMMin, PhiCMMax, PhiSteps,
				 MaxTime, SimStep, Update, Wait);
  }
}

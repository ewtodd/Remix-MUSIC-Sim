/////////////////////////////////////////////////////////////////////////////
// Description: ROOT macro to run the MUSIC_Simulator.
//
// Usage: root -l RunMUSICSim_ap.C
//
// Created by: Daniel Santiago-Gonzalez
// Date: Oct 2019
/////////////////////////////////////////////////////////////////////////////
#include <unistd.h>  // needed for getcwd()

void RunMUSICSim_an()
{
  // These 3 lines get the current working directory (to load the SRIM files)
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  string SRIMdir = cwd;

  //=======================================================================
  // CONTROL PANEL (units: dist=cm, energy=MeV, angle=deg)
  string AnodeGeom = "AnodeGeometry";
  int ELossBins = 300;
  float MaxELoss = 2;
  SRIMdir += "/SRIM_files/";
  // beam
  string beam = "10B";
  string SRIMbeam = SRIMdir + "10B_in_GasIndex3_130Torr_293K.srim";
  // target (no SRIM file, assuming it does not propagate)
  string target = "4He";
  // compound nucleus (no SRIM file, assuming it does not propagate)
  string compound = "14N";
  
  // Reac1: 12C(16C,4n2p)22Ne
  // Chan1: 28Mg*->n+27Mg*->n+26Mg*->n+25Mg*->n+24Mg*->p+23Na*->p+22Ne
  const int NumEvapPart = 1;
  string res[NumEvapPart];
  string SRIMres[NumEvapPart];
  string evap[NumEvapPart];
  string SRIMevap[NumEvapPart];
  evap[0] = "n";    SRIMevap[0] = "";
  res[0] = "13N";   SRIMres[0] = SRIMdir + "13N_in_GasIndex3_130Torr_293K.srim";
  int color[NumEvapPart];
  for (int er=0; er<NumEvapPart; er++) {
    if (evap[er]=="n")
      color[er] = kBlue;
    if (evap[er]=="p")
      color[er] = kRed;
  }

  double Kb = 6.7;   // MeV - Energy of the beam after the Ti window and Al degrader
  int Strip = 2;     // Strip where reaction takes place
  float Eres = 0.01;  // MeV - Strip energy resolution (larger values increase signal randomness)
  int NEvents = 30;   // Number of simulated events (recommendation: keep it <1000)
  int Wait = 0;       // 1 - canvas waits for user's double click, 0 - no wait
  int Update = 0;     // 1 - update visuals for every event, 0 - don't
  double MaxTime = 1000;   // ns - max time for an event
  double UserDT = 0.1;     // ns - simulation time steps
  double ThCMMin = 0.5;
  double ThCMMax = 179.5;
  int ThSteps = 36;
  double PhiCMMin = 0.1;
  double PhiCMMax = 359.9;
  int PhiSteps = 72;
  //=======================================================================


  /////////////////////////////////////////////////////////////////////////////
  // Simulator stuff (typically not modified)
  /////////////////////////////////////////////////////////////////////////////

  MUSIC_Simulator* MUSIC = new MUSIC_Simulator();
  MUSIC->SetPrintLevel(0);
  MUSIC->SetStripEnergyResolution(Eres);
  // Geometry
  MUSIC->SetAnode(AnodeGeom, 90, ELossBins, MaxELoss);
  // Beam
  MUSIC->SetBeamParticle(beam, kBlack, SRIMbeam, Kb);
  // Target
  MUSIC->SetTargetParticle(target);
  // Compound particle
  MUSIC->SetCompoundParticle(compound);
  // Evaporation residues and particles
  for (int i=0; i<NumEvapPart; i++)
    MUSIC->SetEvapResAndPart(res[i], SRIMres[i], kGreen, evap[i], SRIMevap[i], color[i]);

  // Simulate events for one strip or generate trace data base (see below)
  //  MUSIC->Simulate(Strip, NEvents, MaxTime, UserDT, Update, Wait);
  //  MUSIC->WriteTraces(Form("Traces_Stp%d_ap.root", Strip));

  // Uncomment lines below to generate trace database (scanning thorugh angles in CM)
  MUSIC->GenerateTraceDatabase("TraceDB_an.root", 
  			       ThCMMin, ThCMMax, ThSteps, 
  			       PhiCMMin, PhiCMMax, PhiSteps,
  			       MaxTime, UserDT, Update, Wait);
 
}

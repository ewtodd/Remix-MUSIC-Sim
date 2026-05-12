/////////////////////////////////////////////////////////////////////////////
// Description: ROOT macro to run the MUSIC_Simulator.
//
// Usage: root -l RunMUSICSim_16C.C
//
// Created by: Daniel Santiago-Gonzalez
// Date: Dec 2018
/////////////////////////////////////////////////////////////////////////////
#include <unistd.h>  // needed for getcwd()


void RunMUSICSim_16C()
{
  // These 3 lines get the current working directory (to load the SRIM files)
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  string SRIMdir = cwd;

  //=======================================================================
  // CONTROL PANEL (units: dist=cm, energy=MeV, angle=deg)
  string AnodeGeom = "AnodeGeometry";
  int ELossBins = 300;
  float MaxELoss = 25;
  SRIMdir += "/SRIM_files/";
  // beam
  string beam = "16C";
  string SRIMbeam = SRIMdir + "16C_in_GasIndex17_400Torr_293K.srim";
  // target (no SRIM file, assuming it does not propagate)
  string target = "12C";
  // compound nucleus (no SRIM file, assuming it does not propagate)
  string compound = "28Mg";
  
  // Reac1: 12C(16C,4n2p)22Ne
  // Chan1: 28Mg*->n+27Mg*->n+26Mg*->n+25Mg*->n+24Mg*->p+23Na*->p+22Ne
  const int NumEvapPart = 6;
  string res[NumEvapPart];
  string SRIMres[NumEvapPart];
  string evap[NumEvapPart];
  string SRIMevap[NumEvapPart];
  evap[0] = "n";    SRIMevap[0] = SRIMdir + "";
  res[0] = "27Mg";  SRIMres[0] = SRIMdir + "27Mg_in_GasIndex17_400Torr_293K.srim";
  evap[1] = "n";    SRIMevap[1] = SRIMdir + "";
  res[1] = "26Mg";  SRIMres[1] = SRIMdir + "26Mg_in_GasIndex17_400Torr_293K.srim";
  evap[2] = "n";    SRIMevap[2] = SRIMdir + "";
  res[2] = "25Mg";  SRIMres[2] = SRIMdir + "25Mg_in_GasIndex17_400Torr_293K.srim";
  evap[3] = "n";    SRIMevap[3] = SRIMdir + "";
  res[3] = "24Mg";  SRIMres[3] = SRIMdir + "24Mg_in_GasIndex17_400Torr_293K.srim";
  evap[4] = "p";    SRIMevap[4] = SRIMdir + "1H_in_GasIndex17_400Torr_293K.srim";
  res[4] = "23Na";  SRIMres[4] = SRIMdir + "23Na_in_GasIndex17_400Torr_293K.srim";
  evap[5] = "p";    SRIMevap[5] = SRIMdir + "1H_in_GasIndex17_400Torr_293K.srim";
  res[5] = "22Ne";  SRIMres[5] = SRIMdir + "22Ne_in_GasIndex17_400Torr_293K.srim";
  int color[NumEvapPart];
  for (int er=0; er<NumEvapPart; er++) {
    if (evap[er]=="n")
      color[er] = kBlue;
    if (evap[er]=="p")
      color[er] = kRed;
  }

  double Kb = 61.6;   // MeV - Energy of the beam after the Ti window and Al degrader
  int Strip = 11;     // Strip where reaction takes place
  float Eres = 0.05;  // MeV - Strip energy resolution (larger values increase signal randomness)
  int NEvents = 30;   // Number of simulated events (recommendation: keep it <1000)
  int Wait = 1;       // 1 - canvas waits for user's double click, 0 - no wait
  int Update = 1;     // 1 - update visuals for every event, 0 - don't
  double MaxTime = 1000;   // ns - max time for an event
  double UserDT = 0.1;     // ns - simulation time steps
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
  // Simulate events for one strip
  MUSIC->Simulate(Strip, NEvents, MaxTime, UserDT, Update, Wait);
  MUSIC->WriteTraces(Form("Traces_Stp%d.root", Strip));

}

/////////////////////////////////////////////////////////////////////////////
// Description: ROOT macro to run the MUSIC_Simulator.
//
// Usage: root -l RunMUSICSim.C
//
// Created by: Daniel Santiago-Gonzalez
// Date: Oct 2019
/////////////////////////////////////////////////////////////////////////////
#include <unistd.h>  // needed for getcwd()

void RunMUSICSim()
{
  // These 3 lines get the current working directory (to load the SRIM files)
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  string SRIMdir = cwd;

  //=======================================================================
  // CONTROL PANEL (units: dist=cm, energy=MeV, angle=deg)
  string reaction = "(a,n)";  // Options (a,a), (a,n), (a,2n), (a,p)
  int pressure = 600; // Torr
  string AnodeGeom = "AnodeGeometry";
  int ELossBins = 300;
  float MaxELoss = 25;
  SRIMdir += "/SRIM_files/";
  // beam
  string beam = "59Fe";
  string SRIMbeam = SRIMdir + "59Fe_in_GasIndex3_600Torr_293K.srim";
  // target (no SRIM file, assuming it does not propagate)
  string target = "4He";
  // compound nucleus (no SRIM file, assuming it does not propagate)
  string compound = "63Ni";
  
  // Reac1: 12C(16C,4n2p)22Ne
  // Chan1: 28Mg*->n+27Mg*->n+26Mg*->n+25Mg*->n+24Mg*->p+23Na*->p+22Ne
  int NumEvapPart = 0;
  if (reaction=="(a,a)" || reaction=="(a,n)" ||reaction=="(a,p)")
    NumEvapPart = 1;
  else if (reaction=="(a,2n)")
    NumEvapPart = 2;
  
  string* res = new string[NumEvapPart];
  string* SRIMres = new string[NumEvapPart];
  string* evap = new string[NumEvapPart];
  string* SRIMevap = new string[NumEvapPart];
  if (reaction=="(a,a)") {
    evap[0] = "4He";    
    res[0] = beam;
  }
  else if (reaction=="(a,p)") {
    evap[0] = "1H";    
    res[0] = "62Co";    
  }
  else if (reaction=="(a,n)") {
    evap[0] = "n";
    res[0] = "62Ni";        
  }
  else if (reaction=="(a,2n)") {
    evap[0] = "n";
    res[0] = "62Ni";    
    evap[1] = "n";
    res[1] = "61Ni";
  }

  SRIMres[0] = SRIMdir + Form("%s_in_GasIndex3_%dTorr_293K.srim",res[0].c_str(), pressure);
  SRIMevap[0] = SRIMdir + Form("%s_in_GasIndex3_%dTorr_293K.srim",evap[0].c_str(), pressure);
  if (NumEvapPart>1) {
    SRIMres[1] = SRIMdir + Form("%s_in_GasIndex3_%dTorr_293K.srim",res[1].c_str(), pressure);
    SRIMevap[1] = SRIMdir + Form("%s_in_GasIndex3_%dTorr_293K.srim",evap[1].c_str(), pressure);
  }

  int color[NumEvapPart];
  for (int er=0; er<NumEvapPart; er++) {
    if (evap[er]=="n")
      color[er] = kBlue;
    else if (evap[er]=="1H")
      color[er] = kRed;
    else if (evap[er]=="4He")
      color[er] = kMagenta;
    else
      color[er] = kGray+2;
  }


  double Kb = 204;    // MeV - Energy of the beam after the Ti window and Al degrader
  int strip = 4;     // Strip where reaction takes place
  float Eres = 0.005;  // MeV - Strip energy resolution (larger values increase signal randomness)
  int NEvents = 300;   // Number of simulated events (recommendation: keep it <1000)
  int Wait = 0;       // 1 - canvas waits for user's double click, 0 - no wait
  int Update = 0;     // 1 - update visuals for every event, 0 - don't
  double MaxTime = 1000;   // ns - max time for an event
  double SimStep = 0.001;     // cm - simulation steps size
  int Method = 0;    // Select the simulation method: 0 - Simulate, 1 - GenerateTraceDatabase
  string FileName = Form("Traces_Stp%d_%sMeV.root", strip, reaction.c_str());
  //  string FileName = Form("Traces_59Fe_beam_230MeV_test.root");
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
 
  if (Method==0) {
    // Simulate events for one strip or generate trace data base (see below)
    MUSIC->Simulate(strip, NEvents, MaxTime, SimStep, Update, Wait, FileName, FileOpt);
  }
  else if (Method==1) {
    MUSIC->GenerateTraceDatabase("TraceDB.root", 
				 ThCMMin, ThCMMax, ThSteps, 
				 PhiCMMin, PhiCMMax, PhiSteps,
				 MaxTime, SimStep, Update, Wait);
  }
}

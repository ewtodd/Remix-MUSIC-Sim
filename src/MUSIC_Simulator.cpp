// Methods for MUSIC_Simulator class.
// See header file for class description and compilation instructions.

#include "MUSIC_Simulator.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////////
MUSIC_Simulator::MUSIC_Simulator()
{
  Name = "MUSIC_Simulator";

  // Initialize variables.
  CMEMax = 0;
  CMEMin = 0;
  EexcMax = 0;
  EexcMin = 0;
  Kb_after_window = 0;
  NTraces = 0;
  PrintLevel = 0;
  SegLength = 0;
  SegCMERange = 0;
  SegEexcRange = 0;

  // Initialize selected pointers to zero (non-zero pointers initialized below).
  Beam = 0;
  Gaussian = new TF1("Gaussian","gaus",-1e2,1e2);
  Gaussian->SetNpx(1000);
  Target = 0;
  Trace = 0;
  Compound = 0;
  Light = 0;
  Heavy = 0;
  DeDau1 = 0;
  DeDau2 = 0;

  // Evaporated particles and residues arrays
  MaxEva = 20;                       // currently only allowing 20 evaporated particles
  CurEva = 0;                        // This index will increase when a new evap particle is created
  EvaP = new Particle*[MaxEva];
  EvaR = new Particle*[MaxEva];
  Kl = new float[MaxEva];
  Kh = new float[MaxEva];
  theta_CM = new float[MaxEva];
  phi_CM = new float[MaxEva];
  theta_l = new float[MaxEva];
  phi_l = new float[MaxEva];
  theta_h = new float[MaxEva];
  phi_h = new float[MaxEva];
  xfl = new float[MaxEva];
  yfl = new float[MaxEva];
  zfl = new float[MaxEva];
  for (int er=0; er<MaxEva; er++) {
    Kl[er] = 0;
    Kh[er] = 0;
    phi_CM[er] = -1;
    theta_CM[er] = -1;
    phi_l[er] = -1;
    theta_l[er] = -1;
    phi_h[er] = -1;
    theta_h[er] = -1;
    xfl[er] = 0;
    yfl[er] = 0;
    zfl[er] = -1000;
  }
  xr = yr = 0;
  zr = -1000;
  xfe = yfe = 0;
  zfe = -1000;

  // Nuclide finder object
  NuF = new NuclideFinder();

  // Geometry manager
  Geo = new TGeoManager("Geo", "MUSIC geometry manager");

  // Define some materials and media
  MatVacuum = new TGeoMaterial("Vac", 0, 0, 0);
  // NOTE: Not sure about units of arguments
  Vacuum = new TGeoMedium("Vacuum", 1, MatVacuum);
  
  // Make the top container volume
  VolTop = Geo->MakeBox("VolTop", Vacuum, 300., 300., 300.);
  Geo->SetTopVolume(VolTop);
  // Zero other volume pointers
  VolAnode = 0;
  
  // Pseudo-random number generator.
  Rdm = new TRandom3();
  Rdm->SetSeed();     // Provide a seed that depends on the time.
  
  // Arrays for the TTree
  SimTree = 0;
  de_r = new float[ExpAnodeStps];
  de_l = new float[ExpAnodeStps];
  seg = new int[ExpAnodeStps];
  for (int stp=0; stp<ExpAnodeStps; stp++) {
    de_r[stp] = 0;
    de_l[stp] = 0;
    seg[stp] = -1;
  }

  InitCTF(); // Initialize control file parameters.
  
  gSystem = 0;
}


///////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::CalculateCMEnergyRange()
{
#if 0
  double pb, Eb;
  double CME_beg, CME_end, pCM;
  double mb = Beam->Mass;
  double mt = Target->Mass;
  FourVector Pb("Pb");
  FourVector Pt("Pt", mt, 0, 0, 0);
  FourVector Ptot("Total four-mom. in the lab");
  double Kb = ctf.Kb;
  double Kb_min;
  float TotalLength = 0;
  
  // Linear momentum and total energy of the beam particle in the lab with the current
  // value of the kinetic energy.
  pb = sqrt(2*mb*Kb*(1 + Kb/(2*mb)));
  Eb = sqrt(mb*mb + pb*pb);
  Pb.SetCoords(Eb, 0, 0, pb);
  // Total four momentum
  Ptot = Pb + Pt;
    // Calculate the initial momentum (of beam and target) in the CM.
    pCM = sqrt((Ptot*Ptot - pow(mt+mb,2))*(Ptot*Ptot - pow(mt-mb,2))/(4*(Ptot*Ptot)));
    // Center-of-mass energy at the beginning of MUSIC (the CM energy is the sum of the kinetic
    // energies).
    CMEMax = CME_beg = (sqrt(mb*mb+pCM*pCM) - mb) + (sqrt(mt*mt+pCM*pCM) - mt);
    //    cout << "Full CM energy range covered in " << TotalLength << "cm:\n CME(beg) = " 
    //	 << CME_beg << " MeV";
    // Now get the CM energy at the end of the segments.
    Kb_min = BeamInTgt->GetFinalEnergy(Kb, TotalLength, 0.001);
    pb = sqrt(2*mb*Kb_min*(1 + Kb_min/(2*mb)));
    Eb = sqrt(mb*mb + pb*pb);
    Pb.SetCoords(Eb, 0, 0, pb);
    Ptot = Pb + Pt;
    pCM = sqrt((Ptot*Ptot - pow(mt+mb,2))*(Ptot*Ptot - pow(mt-mb,2))/(4*(Ptot*Ptot)));
    CMEMin = CME_end = (sqrt(mb*mb+pCM*pCM) - mb) + (sqrt(mt*mt+pCM*pCM) - mt);
    //   cout << "   CME(end) = " << CME_end << " MeV" << endl;

    //   cout << "CM energy range in each segment:" << endl;
   

    for (int i=0; i<NSegments; i++) {
      Kb = BeamInTgt->GetFinalEnergy(Kb, SegLength[i], 0.001);
      // Linear momentum and total energy of the beam particle in the lab with the current
      // value of the kinetic energy.
      pb = sqrt(2*mb*Kb*(1 + Kb/(2*mb)));
      Eb = sqrt(mb*mb + pb*pb);
      // Four-momentum of the beam.
      Pb.SetCoords(Eb, 0, 0, pb);
      // Total four momentum
      Ptot = Pb + Pt;
      // Calculate the initial momentum (of beam and target) in the CM.
      pCM = sqrt((Ptot*Ptot - pow(mt+mb,2))*(Ptot*Ptot - pow(mt-mb,2))/(4*(Ptot*Ptot)));
      // In this case the CM energy is the sum of the kinetic energies.
      CME_end = (sqrt(mb*mb+pCM*pCM) - mb) + (sqrt(mt*mt+pCM*pCM) - mt);
      SegCMERange[i] = CME_beg - CME_end;
      // cout << i << "\t" << SegLength[i] << " cm \t" << SegCMERange[i] << " MeV \t(Kb=" 
      // 	   << Kb << " MeV)"<< endl;
      // The excitation energy at the end of this segment is the excitation energy at the beginnig
      // of the next segment.
      if (i+1<NSegments) 
	CME_beg = CME_end;
    }
  }
  else {
    cout << "Warning: energy loss file for beam has not been loaded.  Calculations were not made.\n";
  }
#endif
  return;
}


///////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::CalculateExcEnergyRange()
{
#if 1
  double pb, Eb;
  double Eexc_beg, Eexc_end;
  double mb = Beam->Mass;
  double mt = Target->Mass;
  double mf = Compound->Mass;
  FourVector Pb("Pb");
  FourVector Pt("Pt", mt, 0, 0, 0);
  FourVector Ptot("Total four-mom. in the lab");
  double Kb = ctf.Kb;
  double Kb_min;
  float TotalLength = 0;
  
 
  for (int i=0; i<AnodeRows; i++) 
    TotalLength += AnodeDZ[i][0];
    
  // Linear momentum and total energy of the beam particle in the lab with the current
  // value of the kinetic energy.
  pb = sqrt(2*mb*Kb*(1 + Kb/(2*mb)));
  Eb = sqrt(mb*mb + pb*pb);
  Pb.SetCoords(Eb, 0, 0, pb);
  // Total four momentum
  Ptot = Pb + Pt;
  // Energy of excitation at the begining of MUSIC.
  EexcMax = Eexc_beg = sqrt(Ptot*Ptot) - mf;
  cout << "Full excitation energy range covered in " << TotalLength << "cm:\n Eexc(beg) = " 
       << Eexc_beg << " MeV";
  // Now get the energy of excitation at the end of the segments.
  Kb_min = Beam->GetFinalEnergy(0, Kb, TotalLength, 0.001);
  pb = sqrt(2*mb*Kb_min*(1 + Kb_min/(2*mb)));
  Eb = sqrt(mb*mb + pb*pb);
  Pb.SetCoords(Eb, 0, 0, pb);
  Ptot = Pb + Pt;
  EexcMin = Eexc_end = sqrt(Ptot*Ptot) - mf;
  cout << "   Eexc(end) = " << Eexc_end << " MeV" << endl;
  cout << "Exc. energy range in each segment:" << endl;
  for (int i=0; i<AnodeRows; i++) {
    Kb = Beam->GetFinalEnergy(0, Kb, AnodeDZ[i][0], 0.001);
    // Linear momentum and total energy of the beam particle in the lab with the current
    // value of the kinetic energy.
    pb = sqrt(2*mb*Kb*(1 + Kb/(2*mb)));
    Eb = sqrt(mb*mb + pb*pb);
    // Four-momentum of the beam.
    Pb.SetCoords(Eb, 0, 0, pb);
    // Total four momentum
    Ptot = Pb + Pt;
    Eexc_end = sqrt(Ptot*Ptot) - mf;
    SegEexcRange[i] = Eexc_beg - Eexc_end;
    cout << i << "\t" << AnodeDZ[i][0] << " cm \t" << SegEexcRange[i] << " MeV" << endl;
    // The excitation energy at the end of this segment is the excitation energy at the beginnig
    // of the next segment.
    if (i+1<AnodeRows) 
      Eexc_beg = Eexc_end;
  }
#endif
  return;
}

///////////////////////////////////////////////////////////////////////////////////
// Simple function that returns 0 if the used RAM is more than the memory limit
// (or if the 'System' pointer hasn't been set) and 1 if the used RAM is 
// still less than the memory limit.
///////////////////////////////////////////////////////////////////////////////////
int MUSIC_Simulator::CheckMemoryUsage(int Print)
{
  int status = 1;
  MemInfo_t* Memory;
  if (gSystem!=0) {
    Memory = new MemInfo_t();
    gSystem->GetMemInfo(Memory);
    float MemoryLimit = 0.95*Memory->fMemTotal;
    if (Print) {
      cout << "> Total memory used: "  << Memory->fMemUsed << " MB = "
	   << 100.0*Memory->fMemUsed/Memory->fMemTotal 
	   << " % of max memory (" << Memory->fMemTotal << " MB)" << endl;
    }
    if ((float)Memory->fMemUsed > MemoryLimit) {
      cout << "Memory limit exceded! Limit at " << MemoryLimit << " MB" << endl;
      status = 0;
    }
    delete Memory;
  }
  else
    cout << "> Warning: CheckMemoryUsage(), gSystem pointer not set. "
	 << "Use SetROOTSystemPointer() method."
	 << endl;
  return status;
}



///////////////////////////////////////////////////////////////////////////////////
// Private method, to be used in an event loop.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::ComputeDetectorResponse(int evt, int reacStp, int UpdateVis)
{
#if 1
  if (PrintLevel>0)
    Log << "Compute detector response evt " << evt << endl;
  double Ethresh = 0.02;  // Right now the threshold for the multiplicity is hard-coded here.
  double DeltaE = 0;
  TraceMult->Reset();
  
  // Reset the SimTree leaves (already reset in Simulate())
  if (SimTree!=0) {
     this->reacStp = reacStp;
  }
  
  // Loop over the anode's rows and columns
  for (int row=0; row<AnodeRows; row++) {
    DeltaE = 0;
    int mult = 0;
    for (int col=0; col<AnodeCols+1; col++) {
      DeltaE = DeltaEB[row][col];
      for (int er=0; er<CurEva; er++) {
	DeltaE += DeltaE_EvaP[er][row][col];
	DeltaE += DeltaE_EvaR[er][row][col];
      }

      if (ctf.Eres>0.0) {
	// Adding randomness to the energy loss to mimic experimental jitter
	Gaussian->SetParameters(1.0, 0.0, ctf.Eres);
	DeltaE += Gaussian->GetRandom();
      }	

      if (DeltaE>Ethresh && col<AnodeCols)
	mult++;

      // Fill tree leaves
      if (SimTree!=0) {
	int stpid = AnodeStpID[row][col];
	if (stpid>=0 && stpid<=17) {
	  cathode += DeltaE;
	  if (stpid==0 && col==0)
	    strip0 += DeltaE;
	  else if (stpid==17 && col==0)
	    strip17 += DeltaE;
	  else if (stpid-1<ExpAnodeStps) {
	    seg[stpid-1] = stpid;
	    if (col==0) 
	      de_r[stpid-1] += DeltaE;
	    else if (col==1)
	      de_l[stpid-1] += DeltaE;
	  }
	}
      }
      // Changed randomness to PropagateParticle method
      // if (EneSigma!=0 && Gaussian!=0 && DeltaE>0) {
      // 	Gaussian->SetRange(0.0, 2*DeltaE);
      // 	Gaussian->SetParameters(1.0, DeltaE, EneSigma);
      // 	DeltaE = Gaussian->GetRandom();
      // }
      Trace[col]->SetPoint(row, row, DeltaE);
      // if (PrintLevel>0)
      // 	cout << row << " " << col << ": " << DeltaE << " = " << DeltaEB[row][col] << "+" 
      // 	     << DeltaEL[row][col] << "+" << DeltaEH[row][col] << "+" << DeltaED1[row][col] 
      // 	     << "+" << DeltaED2[row][col] << "+E.R."<< endl;
    }
    TraceMult->Fill(row, mult);
  }
  

  if (Trace!=0) {
    for (int col=0; col<AnodeCols+1; col++) {
      if (col==AnodeCols) 
	Trace[col]->Write(Form("Trace%d",evt), TObject::kOverwrite);
      else
	Trace[col]->Write(Form("Trace%dc%d",evt,col), TObject::kOverwrite);
    }
  }
  
#endif
  return;
}

///////////////////////////////////////////////////////////////////////////////////
// Create energy loss trace for this event (detector response). The last column has
// the energy loss deposited in the whole strip.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::CreateTracesAndTrajectories(int NEvents)
{
#if 1
  // Look for the first strip in which all the column colors are
  // specified and assign those colors to the Chroma elements.
  short* Chroma = new short[AnodeCols];
  for (int col=0; col<AnodeCols; col++)
    Chroma[col] = 7 - col;
  for (int stp=0; stp<AnodeRows; stp++) {
    int NotWhiteColumns = 0;
    for (int col=0; col<AnodeCols; col++)
      if (AnodeColor[stp][col]!=kWhite)
	NotWhiteColumns++;
    if (NotWhiteColumns==AnodeCols) {
      for (int col=0; col<AnodeCols; col++)
	Chroma[col] = AnodeColor[stp][col];
      break;
    }
  }

  // Initialize the 'detector' traces for each column
  Trace = new TGraph*[AnodeCols+1];
  for (int col=0; col<AnodeCols+1; col++) {
    Trace[col] = new TGraph();
    if (col==AnodeCols) {
      Trace[col]->SetName("full trace");
      Trace[col]->SetLineColor(kBlack);
      Trace[col]->SetLineWidth(2);
    }
    else {
      Trace[col]->SetName(Form("trace col %d", col));
      Trace[col]->SetLineColor(Chroma[col]);
      Trace[col]->SetLineStyle(2);
      Trace[col]->SetLineWidth(2);
    }
  }

  TraceMult = new TH1F("TraceMult","Mult.",AnodeRows,-0.5,AnodeRows-0.5);
  TraceMult->GetXaxis()->SetTitle("Strip index (in AnodeGeometry file)");
  TraceMult->GetXaxis()->CenterTitle();
 
  // Initialize the traces corresponding to the unreacted beam (UB)
  TraceUB = new TGraph*[AnodeCols+1];
  for (int col=0; col<AnodeCols+1; col++) {
    TraceUB[col] = new TGraph();
    if (col==AnodeCols) {
      TraceUB[col]->SetName("Beam trace");
      TraceUB[col]->SetLineColor(kGray);
      TraceUB[col]->SetLineWidth(3);
    }
    else {
      TraceUB[col]->SetName(Form("Beam trace col %d", col));
      TraceUB[col]->SetLineColor(kGray);
      TraceUB[col]->SetLineStyle(2);
      TraceUB[col]->SetLineWidth(2);
    }
  }

  // Initialize the traces corresponding to the beam particle
  TraceB = new TGraph*[AnodeCols+1];
  for (int col=0; col<AnodeCols+1; col++) {
    TraceB[col] = new TGraph();
    if (col==AnodeCols) {
      TraceB[col]->SetName("Beam trace");
      TraceB[col]->SetLineColor(kGray+2);
      TraceB[col]->SetLineWidth(3);
    }
    else {
      TraceB[col]->SetName(Form("Beam trace col %d", col));
      TraceB[col]->SetLineColor(kGray);
      TraceB[col]->SetLineStyle(2);
      TraceB[col]->SetLineWidth(2);
    }
  }
  
  // Initialize the traces corresponding to the evaporation residues
  // (heavy particles)
  TraceER = new TGraph**[CurEva];
  for (int er=0; er<CurEva; er++) {
    TraceER[er] = new TGraph*[AnodeCols+1];
    for (int col=0; col<AnodeCols+1; col++) {
      TraceER[er][col] = new TGraph();
      if (col==AnodeCols) {
	TraceER[er][col]->SetName(Form("%s trace",EvaR[er]->Name.c_str()));
	TraceER[er][col]->SetLineColor(EvaR[er]->GetColor());
	TraceER[er][col]->SetLineWidth(3);
      }
      else {
	TraceER[er][col]->SetName(Form("%s trace col %d", EvaR[er]->Name.c_str(), col));
	TraceER[er][col]->SetLineColor(EvaR[er]->GetColor()-er-1);
	TraceER[er][col]->SetLineStyle(2);
	TraceER[er][col]->SetLineWidth(2);
      }
    }
  }

  // Initialize the traces corresponding to the evaporated particles
  // (light particles, p, n, 4He)
  TraceEP = new TGraph**[CurEva];
  for (int er=0; er<CurEva; er++) {
    TraceEP[er] = new TGraph*[AnodeCols+1];
    for (int col=0; col<AnodeCols+1; col++) {
      TraceEP[er][col] = new TGraph();
      if (col==AnodeCols) {
	TraceEP[er][col]->SetName(Form("%s trace",EvaP[er]->Name.c_str()));
	TraceEP[er][col]->SetLineColor(EvaP[er]->GetColor());
	TraceEP[er][col]->SetLineWidth(3);
      }
      else {
	TraceEP[er][col]->SetName(Form("%s trace col %d", EvaP[er]->Name.c_str(), col));
	TraceEP[er][col]->SetLineColor(EvaP[er]->GetColor()-er-1);
	TraceEP[er][col]->SetLineStyle(2);
	TraceEP[er][col]->SetLineWidth(2);
      }
    }
  }

#endif
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Show a 3D image of MUSIC with an given transparency level into the specified
// TEveManager object pointer.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::DrawMUSIC(TEveManager* gEve, short Transparency /*From 0 to 100*/)
{
#if 1
  if (VolAnode!=0) {
    if (Transparency<0 || Transparency>100) {
      cout << "Warning: Transparency level must be from 0 to 100." << endl;
      Transparency = 0;
    }
    
    // From: http://root.cern.ch/root/html/TGeoManager.html#TGeoManager:CloseGeometry
    // Closing geometry implies checking the geometry validity, fixing shapes
    // with negative parameters (run-time shapes) building the cache manager,
    // voxelizing all volumes, counting the total number of physical nodes and
    // registring the manager class to the browser.
    Geo->CloseGeometry();
    
    TopNode = new TEveGeoTopNode(Geo, Geo->GetTopNode());
    gEve->AddGlobalElement(TopNode);
    // Axes
    TEveArrow* Xaxis = new TEveArrow(20,0,0,-10,0,0);
    Xaxis->SetName("x axis");
    Xaxis->SetMainColor(kGray); 
    Xaxis->SetMainTransparency(65);
    Xaxis->SetTubeR(0.01);
    gEve->AddElement(Xaxis);
    TEveArrow* Yaxis = new TEveArrow(0,20,0,0,-10,0);
    Yaxis->SetName("y axis");
    Yaxis->SetMainColor(kYellow);
    Yaxis->SetMainTransparency(65);
    Yaxis->SetTubeR(0.01);
    gEve->AddElement(Yaxis);
    TEveArrow* Zaxis = new TEveArrow(0,0,1,0,0,0);
    Zaxis->SetName("z axis");
    Zaxis->SetMainColor(kWhite);
    Zaxis->SetTubeR(0.1);
    Zaxis->SetMainTransparency(65);
    gEve->AddElement(Zaxis);

    gEve->Redraw3D(kTRUE);
  }
#endif
  return;
}


///////////////////////////////////////////////////////////////////////////////////
//   _____                           _         _______                 
//  / ____|                         | |       |__   __|                
// | |  __  ___ _ __   ___ _ __ __ _| |_ ___     | |_ __ __ _  ___ ___ 
// | | |_ |/ _ \ '_ \ / _ \ '__/ _` | __/ _ \    | | '__/ _` |/ __/ _ \
// | |__| |  __/ | | |  __/ | | (_| | ||  __/    | | | | (_| | (_|  __/
//  \_____|\___|_| |_|\___|_|  \__,_|\__\___|    |_|_|  \__,_|\___\___|
//  _____        _        _                    
// |  __ \      | |      | |                   
// | |  | | __ _| |_ __ _| |__   __ _ ___  ___ 
// | |  | |/ _` | __/ _` | '_ \ / _` / __|/ _ \
// | |__| | (_| | || (_| | |_) | (_| \__ \  __/
// |_____/ \__,_|\__\__,_|_.__/ \__,_|___/\___|
//
// Angles in degrees
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::GenerateTraceDatabase(string FileName, 
					    double ThCMMin, double ThCMMax, int ThSteps,
					    double PhiCMMin, double PhiCMMax, int PhiSteps,
					    double MaxTime, double UserStep, int UpdateVis,
					    int Wait)
{
  double ti,xi,yi,zi, tf,xf,yf,zf;
#if 1
  // Verify that the anode geometry has been set
  if (VolAnode==0) {
    cout << "Anode geometry not specified. Use SetAnode method." << endl;
    return;
  }
  
  // For progress monitor
  TStopwatch StpWatch;
  long EvtsProcessed = 0;
  long double Frac[6] = {0.01, 0.25, 0.5, 0.75, 0.9, 1.0};
  int FIndex = 0;
  
  // ROOT file where the traces will be saved
  TFile* TDB = new TFile(FileName.c_str(), "recreate");

  // Tree similar to the one used for experimental data
  SimTree = InitTree(TDB, "recreate");

  // Angles in radians
  double theta = 0;
  double phi = 0;
  
  NEvents = PhiSteps*ThSteps*AnodeRows;
  NTraces = 0;
    
  // Create new traces and trajectories (objectrs) for visualizing the
  // detector response
  CreateTracesAndTrajectories(NEvents);
  
  cout << "Generating " << NEvents << " MUSIC traces ..." << endl;
  
  SetInitialKinematics(ctf.Kb);   

  // Get the average beam energy loss and print the exc energy of the
  // compound nucleus sampled by each strip.
  Particle* BeamInit = new Particle("beam init");
  BeamInit->Copy(Beam);
  Particle* BeamCopy = new Particle("beam copy");
  BeamCopy->Copy(Beam);
  // if (PrintLevel>0)
  //   BeamCopy->Print(Log);
  //  DeltaEB_ave = PropagateParticle(BeamCopy, ctf.Kb, MaxTime, UserStep); 
  PropagateParticle(BeamCopy, 0, MaxTime, UserStep, DeltaEB_ave);
  for (int stp=0; stp<AnodeRows; stp++)
    for (int col=0; col<AnodeCols+1; col++) 
      TraceUB[col]->SetPoint(stp, stp, DeltaEB_ave[stp][col]);


  //  PrintEnergetics(ctf.Kb, DeltaEB_ave);

  //-------------------------------------------------------------------------------
  // Some kinematic variables
  double Kb_min, Kb_max, MinZ, MaxZ, MinT, MaxT;
  int evt = 0; 
  // Converting degrees in radians
  double theta_min = ThCMMin*pi/180;
  double theta_max = ThCMMax*pi/180;
  double phi_min = PhiCMMin*pi/180;
  double phi_max = PhiCMMax*pi/180;

  // Loop over all strips with a non-negative ID (defined in the
  // argument of the SetAnode method).
  //for (int stp_base=11; stp_base<AnodeRows; stp_base++) {
  for (int stp_base=0; stp_base<AnodeRows; stp_base++) {
    if (AnodeStpID[stp_base][0]>=0) {
      
      // Get the beam energy limits in the selected strip (assuming
      // the beam direction is parallel to the z-axis).
      MinZ = 0;
      MaxZ = AnodeDZ[0][0];
      for (int stp=1; stp<AnodeRows; stp++) {
	MinZ += AnodeDZ[stp-1][0];
	MaxZ += AnodeDZ[stp][0];
	if (AnodeStpID[stp][0]==AnodeStpID[stp_base][0])
	  break;
      }
      Kb_max = Beam->GetFinalEnergy(0, ctf.Kb, MinZ, 1E-3/*step size in cm*/);
      MinT = Beam->GetTimeOfFlight(0);
      Kb_min = Beam->GetFinalEnergy(0, ctf.Kb, MaxZ, 1E-3/*step size in cm*/);
      MaxT = Beam->GetTimeOfFlight(0);      
      if (PrintLevel>0) {
	Log << "|---- Kinematic constraints for strip ";
	Log.width(3); Log << AnodeStpID[stp_base][0];
	Log << " ---------\n"
	     << "|     |   In    |   Out   |  Units  |\n"
	     << "| zr  |";
	Log.width(9);    Log << MinZ;   Log << "|";
	Log.width(9);    Log << MaxZ;   Log << "|";
	Log.width(9);    Log << "cm";   Log << "|\n";
	Log << "| tof |";
	Log.width(9);    Log << MinT;   Log << "|";
	Log.width(9);    Log << MaxT;   Log << "|";
	Log.width(9);    Log << "ns";   Log << "|\n";
	Log << "| Kb  |";
	Log.width(9);    Log << Kb_max; Log << "|";
	Log.width(9);    Log << Kb_min; Log << "|";
	Log.width(9);    Log << "MeV";   Log << "|\n";
	Log << "|--------------------------------------------------" << endl; 
      }
      for (int ths=0; ths<ThSteps; ths++) {
	theta = ths*(theta_max - theta_min)/(ThSteps - 1) + theta_min;
	for (int phs=0; phs<PhiSteps; phs++) {
	  phi = phs*(phi_max - phi_min)/(PhiSteps - 1) + phi_min;
	  
	  if (PrintLevel>0) {
	    Log << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	    Log << "!!!       EVENT " << evt << endl;
	    Log << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" << endl;
	  }
	  
	  TraceCan->cd(1);
	  HCT->Draw();
	  
	  // Reset the detector response
	  for (int stp=0; stp<AnodeRows; stp++) 
	    for (int col=0; col<AnodeCols+1; col++) {
	      DeltaEB[stp][col] = 0;
	      DeltaEL[stp][col] = 0;
	      DeltaEH[stp][col] = 0;
	      for (int er=0; er<CurEva; er++) {
		DeltaE_EvaP[er][stp][col] = 0;
		DeltaE_EvaR[er][stp][col] = 0;
	      }
	    }
	  // 1. Set beam inital conditions (beam energy, position)
	  SetInitialKinematics(ctf.Kb);   
	  
	  // 2. Within the selected strip randomly select the position
	  // at which the beam particle interacts with the target and
	  // calculate the kinetic energy at the reaction point
	  this->zr = Rdm->Uniform(MinZ, MaxZ);
	  double Kbr = Beam->GetFinalEnergy(0, ctf.Kb, this->zr, 1E-3/*cm*/);
	  double TOF = Beam->GetTimeOfFlight(0);
	  if (PrintLevel>0)
	    Log << "Kbr = " << Kbr << "  zr = " << this->zr << "  tof = " << TOF << endl;
	  
	  // 3. Set the kinematics of all particles at the reaction point
	  int ReacAllowed = SetReactionKinematics(Kbr, this->zr, TOF, theta, phi);
	  // Check conservation of 4-momentum
	  if (PrintLevel>0) {
	    Log << "Conservation of 4-momentum at reaction point (zr)" 
		<< endl;
	    FourVector Pi("initial 4-momemtum (lab)",0,0,0,0);
	    Pi += Beam->GetP() + Target->GetP();
	    FourVector Pf("final 4-momentum (lab)",0,0,0,0);
	    for (int er=0; er<CurEva; er++) {
	      if (!EvaR[er]->DoNotPropagate)
		Pf += EvaR[er]->GetP();
	      if (!EvaP[er]->DoNotPropagate)
		Pf += EvaP[er]->GetP();
	    }
	    Pi.Print(Log);
	    Pf.Print(Log);
	  }
#if 0
	  if (ReacAllowed==0) {
	    cout << "Warninig: reaction energetically not allowed for event " << evt 
		 << " (Kbr= " << Kbr << " MeV)." << endl;
	    continue;
	  }
	  
	  // 4. Propagate the beam particle (backwards in time) from
	  // the reaction point to the entrance of MUSIC
	  //	  DeltaEB = PropagateParticle(Beam, evt, MaxTime, -UserStep);
	  PropagateParticle(Beam, evt, MaxTime, -UserStep, DeltaEB);
	  
	  // 5. Propagate heavy particle (or decay daughters) and calculate energy 
	  // loss in the anode elements
	  if (DeDau1 && DeDau2) {
	    //	    DeltaED1 = PropagateParticle(DeDau1, evt, MaxTime, UserStep);
	    // DeltaED2 = PropagateParticle(DeDau2, evt, MaxTime, UserStep);
	    PropagateParticle(DeDau1, evt, MaxTime, UserStep, DeltaED1);
	    PropagateParticle(DeDau2, evt, MaxTime, UserStep, DeltaED2);
	  }
	  else {
	    //	    DeltaEH = PropagateParticle(Heavy, evt, MaxTime, UserStep);
	    PropagateParticle(Heavy, evt, MaxTime, UserStep, DeltaEH);
	  }
	  // 6. Propagate light particle and calculate energy loss in the
	  // anode elements
	  //	  DeltaEL = PropagateParticle(Light, evt, MaxTime, UserStep);
	  PropagateParticle(Light, evt, MaxTime, UserStep, DeltaEL);
#endif	  

	  if (ReacAllowed) {
	    // 4. Propagate the beam particle (backwards in time) from the
	    // reaction point to the entrance of MUSIC
	    Beam->GetX(tf,xf,yf,zf);
	    //	  DeltaEB = PropagateParticle(Beam, evt, MaxTime, -UserStep);
	    PropagateParticle(Beam, evt, MaxTime, -UserStep, DeltaEB);
	    Beam->GetX(ti,xi,yi,zi);                      // <- This is not a mistake
	    TrackBeam->SetOrigin(xi,yi,zi);
	    TrackBeam->SetVector(xf-xi,yf-yi,zf-zi);
	  
	    // 5-6. Propagate outgoing particles (evaporation residues)
	    for (int er=0; er<CurEva; er++) {
	      // evaporated (light) particle (p,n,4He)
	      EvaP[er]->GetX(ti,xi,yi,zi);
	      //	  DeltaE_EvaP[er] = PropagateParticle(EvaP[er], evt, MaxTime, UserStep);
	      PropagateParticle(EvaP[er], evt, MaxTime, UserStep, DeltaE_EvaP[er]);
	      EvaP[er]->GetX(tf,xf,yf,zf);
	      TrackEvaP[er]->SetOrigin(xi,yi,zi);
	      TrackEvaP[er]->SetVector(xf-xi,yf-yi,zf-zi);
	      // evaporation residue (heavy particle)
	      EvaR[er]->GetX(ti,xi,yi,zi);
	      //	  DeltaE_EvaR[er] = PropagateParticle(EvaR[er], evt, MaxTime, UserStep);
	      PropagateParticle(EvaR[er], evt, MaxTime, UserStep, DeltaE_EvaR[er]);
	      EvaR[er]->GetX(tf,xf,yf,zf);
	      if (!EvaR[er]->DoNotPropagate) {
		xfe = xf;
		yfe = yf;
		zfe = zf;
	      }	      
	      TrackEvaR[er]->SetOrigin(xi,yi,zi);
	      TrackEvaR[er]->SetVector(xf-xi,yf-yi,zf-zi);
	    }
	  }
	  else {
	    Beam->Copy(BeamInit);
	    PropagateParticle(Beam, evt, MaxTime, UserStep, DeltaEB);
	    cout << "Warninig: reaction energetically not allowed for event " << evt 
		 << " (Kbr= " << Kbr << " MeV)." << endl;
	  }
	  
	  // 7. Compute detector response (i.e. DE for beam + light + heavy)
	  // Clone the particle trajectories
	  ComputeDetectorResponse(evt, stp_base, UpdateVis);
	  SimTree->Fill();
	  
	  // 8. Display trace and particle trajecories   
	  if (UpdateVis) 
	    UpdateVisuals(evt, Kbr, this->zr, TOF, Wait);
	  
	  // Simple progress monitor
	  EvtsProcessed++;
	  if (NEvents>99) {
	    if ((long double)(EvtsProcessed)>=Frac[FIndex]*NEvents) {
	      cout << "\t" << Frac[FIndex]*100 << "% processed (" 
		   << StpWatch.RealTime() << " s)" << endl;
	      StpWatch.Start(kFALSE);
	      FIndex++;
	    }
	  }

	  NTraces++;
	  evt++;
	}
      }
    }
  }
  cout << "Saving traces ..." << endl;
  StpWatch.Print();

  TDB->cd();
  SimTree->Write("", TObject::kSingleKey);
  TDB->Close();
  StpWatch.Stop();
  StpWatch.Print();
#endif
  return;
}




////////////////////////////////////////////////////////////////////////////////////
// Method to load the parameters from the control file
// Control file example below horizontal line (---), first two lines are skipped
// and can be used for comments. Third column can be used for comments:
// ------------------------------------------------------
// Control file for DecayAnalyzer class, run005         | File description
// Parameter      Value     Comment                     | Column description
// colorScheme    1         Int 1,2,3                   | first data line
// ...                                                  | 
////////////////////////////////////////////////////////////////////////////////////
int MUSIC_Simulator::loadCtrlFile(char* fileName)
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
      // Detector parameters
      if (ParName=="AnodeGeom")
	ctf.AnodeGeom = ParVal;
      else if (ParName=="pressure")
	ctf.pressure = atoi(ParVal.c_str());
      else if (ParName=="ELossBins")
	ctf.ELossBins = atoi(ParVal.c_str());
      else if (ParName=="MaxELoss")
	ctf.MaxELoss = atof(ParVal.c_str());
      else if (ParName=="strip")
	ctf.strip = atoi(ParVal.c_str());
      else if (ParName=="Eres")
	ctf.Eres = atof(ParVal.c_str());

      // Beam parameters
      else if (ParName=="beam")
	ctf.beamName = ParVal;
      else if (ParName=="Ebeam" || ParName=="Kb")
	ctf.Kb = atof(ParVal.c_str());
      else if (ParName=="EbeamFWHM" || ParName=="KbFWHM")
	ctf.KbFWHM = atof(ParVal.c_str());
      else if (ParName=="SRIMbeam")
	ctf.SRIMbeam = ParVal;
      else if (ParName=="dEdxScaleBeam")
	ctf.dEdxScaleBeam = atof(ParVal.c_str());
  
      // Target parameters
      else if (ParName=="target")
	ctf.target = ParVal;
      
      // Compound parameters
      else if (ParName=="compound")
	ctf.compound = ParVal;

      // Number of evaporated particles (e.g. 1H, n)
      else if (ParName=="NumEvapPart")
      	ctf.NumEvapPart = atoi(ParVal.c_str());

      // Particle 0
      else if (ParName=="evap0Name")
	ctf.evap[0] = ParVal;
      else if (ParName=="evap0Color")
	ctf.colorEvap[0] = atoi(ParVal.c_str());
      else if (ParName=="SRIMevap0")
	ctf.SRIMevap[0] = ParVal;
      else if (ParName=="dEdxScaleEvap0")
	ctf.dEdxScaleEvap[0] = atof(ParVal.c_str());
 
      // Evaporation residue 0
      else if (ParName=="res0Name")
	ctf.res[0] = ParVal;
      else if (ParName=="res0Color")
	ctf.colorRes[0] = atoi(ParVal.c_str());
      else if (ParName=="SRIMres0")
	ctf.SRIMres[0] = ParVal;
      else if (ParName=="dEdxScaleRes0")
	ctf.dEdxScaleRes[0] = atof(ParVal.c_str());
      
      // Particle 1
      else if (ParName=="evap1Name")
	ctf.evap[1] = ParVal;
      else if (ParName=="evap1Color")
	ctf.colorEvap[1] = atoi(ParVal.c_str());
      else if (ParName=="SRIMevap1")
	ctf.SRIMevap[1] = ParVal;
       else if (ParName=="dEdxScaleEvap1")
	ctf.dEdxScaleEvap[0] = atof(ParVal.c_str());
     
      // Evaporation residue 1
      else if (ParName=="res1Name")
	ctf.res[1] = ParVal;
      else if (ParName=="res1Color")
	ctf.colorRes[1] = atoi(ParVal.c_str());
      else if (ParName=="SRIMres1")
	ctf.SRIMres[1] = ParVal;
      else if (ParName=="dEdxScaleRes1")
	ctf.dEdxScaleRes[1] = atof(ParVal.c_str());

      // Particle 2
      else if (ParName=="evap2Name")
	ctf.evap[2] = ParVal;
      else if (ParName=="evap2Color")
	ctf.colorEvap[2] = atoi(ParVal.c_str());
      else if (ParName=="SRIMevap2")
	ctf.SRIMevap[2] = ParVal;
       else if (ParName=="dEdxScaleEvap2")
	ctf.dEdxScaleEvap[0] = atof(ParVal.c_str());
     
      // Evaporation residue 2
      else if (ParName=="res2Name")
	ctf.res[2] = ParVal;
      else if (ParName=="res2Color")
	ctf.colorRes[2] = atoi(ParVal.c_str());
      else if (ParName=="SRIMres2")
	ctf.SRIMres[2] = ParVal;
      else if (ParName=="dEdxScaleRes2")
	ctf.dEdxScaleRes[2] = atof(ParVal.c_str());

      // DSG - need to generalize the residue/particle lines above and
      // make it go to MaxNumEvepPart. Maybe with sprintf().

      // Simulation parameters
      else if (ParName=="NEvents")
	ctf.NEvents = atoi(ParVal.c_str()); 
      else if (ParName=="Wait")
	ctf.Wait = atoi(ParVal.c_str()); 
      else if (ParName=="Update")
	ctf.Update = atoi(ParVal.c_str());
      else if (ParName=="MaxTime")
	ctf.MaxTime = atof(ParVal.c_str());
      else if (ParName=="SimStep")
	ctf.SimStep = atof(ParVal.c_str());
      else if (ParName=="Method")
	ctf.Method = atoi(ParVal.c_str());
      else if (ParName=="FileName")
      	ctf.FileName = ParVal;
      else if (ParName=="FileOpt")
      	ctf.FileOpt = ParVal;
      else if (ParName=="PrintOpt")
      	ctf.PrintOpt = atoi(ParVal.c_str());      
      else
	cout << "musicsim warning: control file parameter \'" << ParName << "\' not recognized."
	     << endl;

      
#if 0
      else if (ParName=="str")
      	ctf.var = ParVal;
      else if (ParName=="int")
      	ctf.var = atoi(ParVal.c_str());
      else if (ParName=="float")
      	ctf.var = atof(ParVal.c_str());
#endif
    }
    Status = 1;
  } 
  return Status;
}

///////////////////////////////////////////////////////////////////////////////////
// Basic initialization of control file parameters.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::InitCTF()
{
    ctf.pressure = 100; // Torr
    ctf.AnodeGeom = "";
    ctf.ELossBins = 300;
    ctf.MaxELoss = 10.0;
    ctf.beamName = "unassigned beam";
    ctf.SRIMbeam = "";
    ctf.dEdxScaleBeam = 1.0;
    ctf.target = "unassigned target";
    ctf.compound = "unassigned compound";
    ctf.NumEvapPart = ctf.MaxNumEvapPart;
    for (int i=0; i<ctf.NumEvapPart; i++) {
      ctf.res[i] = "unassigned res";
      ctf.SRIMres[i] = "";
      ctf.dEdxScaleRes[i] = 1.0;
      ctf.colorRes[i] = 416;
      ctf.evap[i] = "unassigned evap";
      ctf.SRIMevap[i] = "";
      ctf.dEdxScaleEvap[i] = 1.0;
      ctf.colorEvap[i] = 616;
    }
    ctf.Kb = 100;      // MeV - Energy of the beam after the Ti window and degrader (if any)
    ctf.KbFWHM = 0.0;  // MeV - Beam energy spread (full-width half maximum)
    ctf.strip = 5;     // Strip where reaction takes place
    ctf.Eres = 0.0;    // MeV - Strip energy resolution (larger values increase signal randomness)
    ctf.NEvents = 10;  // Number of simulated events (recommendation: keep it <1000)
    ctf.Wait = 1;      // 1 - canvas waits for user's double click, 0 - no wait
    ctf.Update = 1;    // 1 - update visuals for every event, 0 - don't
    ctf.MaxTime = 2000.0; // ns - max time for an event
    ctf.SimStep = 0.001;  // cm - simulation steps size
    ctf.Method = 0;    // Select the simulation method: 0 - Simulate, 1 - GenerateTraceDatabase
    ctf.FileName = "";
    ctf.FileOpt = "";
    ctf.PrintOpt = 0;
    return;
}
    
///////////////////////////////////////////////////////////////////////////////////
// Initialize the TTree (SimTree) similar to the one used for experimental data.
///////////////////////////////////////////////////////////////////////////////////
TTree* MUSIC_Simulator::InitTree(TFile* ROOTfile, string FileOpt)
{
  TTree* tree;
  // Tree similar to the one used for experimental data
  if (ROOTfile && (FileOpt=="update" || FileOpt=="UPDATE")) {
    tree = (TTree*)ROOTfile->Get("simt");
    tree->SetBranchAddress("de_l",  de_l);
    tree->SetBranchAddress("de_r",  de_r);
    tree->SetBranchAddress("seg",   seg);
    tree->SetBranchAddress("stp0",  &strip0);
    tree->SetBranchAddress("stp17", &strip17);
    tree->SetBranchAddress("cath",  &cathode);
    // The following branches are for physical quantities that at the
    // moment can only be obtained from the simulation
    tree->SetBranchAddress("reacStp",   &reacStp);
    tree->SetBranchAddress("Kb", &Kb);
    tree->SetBranchAddress("Kl", Kl);
    tree->SetBranchAddress("Kh", Kh);
    tree->SetBranchAddress("theta_CM", theta_CM);
    tree->SetBranchAddress("theta_l",  theta_l);
    tree->SetBranchAddress("theta_h",  theta_h);
    tree->SetBranchAddress("phi_l",    phi_l);
    tree->SetBranchAddress("phi_h",    phi_h);
    // Final coordinates of the light evaporated particles in the reaction (arrays)
    tree->SetBranchAddress("xfl",      xfl);
    tree->SetBranchAddress("yfl",      yfl);
    tree->SetBranchAddress("zfl",      zfl);
    // Reaction point coordinates
    tree->SetBranchAddress("xr", &xr);
    tree->SetBranchAddress("yr", &yr);
    tree->SetBranchAddress("zr", &zr);
    // Final coordinates of the lightest evaporation residue in the reaction
    tree->SetBranchAddress("xfe", &xfe);
    tree->SetBranchAddress("yfe", &yfe);
    tree->SetBranchAddress("zfe", &zfe);
  }
  else {
    tree = new TTree("simt","Simulated MUSIC data");
    tree->Branch("de_l",  de_l,     Form("de_l[%d]/F",ExpAnodeStps));
    tree->Branch("de_r",  de_r,     Form("de_r[%d]/F",ExpAnodeStps));
    tree->Branch("seg",   seg,      Form("seg[%d]/I",ExpAnodeStps));
    tree->Branch("stp0",  &strip0,  "stp0/F");
    tree->Branch("stp17", &strip17, "stp17/F");
    tree->Branch("cath",  &cathode, "cath/F");
    // The following branches are for physical quantities that at the
    // moment can only be obtained from the simulation
    tree->Branch("reacStp",   &reacStp,   "reacStp/I");
    tree->Branch("Kb", &Kb, "Kb/F");
    tree->Branch("Kl", Kl,  Form("Kl[%d]/F",MaxEva));
    tree->Branch("Kh", Kh,  Form("Kh[%d]/F",MaxEva));
    tree->Branch("theta_CM", theta_CM, Form("theta_CM[%d]/F",MaxEva));
    tree->Branch("theta_l",  theta_l,  Form("theta_l[%d]/F",MaxEva));
    tree->Branch("theta_h",  theta_h,  Form("theta_h[%d]/F",MaxEva));
    tree->Branch("phi_l",    phi_l,    Form("phi_l[%d]/F",MaxEva));
    tree->Branch("phi_h",    phi_h,    Form("phi_h[%d]/F",MaxEva));
    // Final coordinates of the light evaporated particles in the reaction (arrays)
    tree->Branch("xfl",      xfl,      Form("xfl[%d]/F",MaxEva));
    tree->Branch("yfl",      yfl,      Form("yfl[%d]/F",MaxEva));
    tree->Branch("zfl",      zfl,      Form("zfl[%d]/F",MaxEva));
    // Reaction point coordinates
    tree->Branch("xr", &xr, "xr/F");
    tree->Branch("yr", &yr, "yr/F");
    tree->Branch("zr", &zr, "zr/F");
    // Final coordinates of the lightest evaporation residue in the reaction
    tree->Branch("xfe", &xfe, "xfe/F");
    tree->Branch("yfe", &yfe, "yfe/F");
    tree->Branch("zfe", &zfe, "zfe/F");
  }
  ResetBranches();
  if (PrintLevel>0)
    tree->Print();
  return tree;
}



///////////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::PrintEnergetics(double Kb, double** DeltaEB)
{
#if 1
  double mb = Beam->Mass;
  double mt = Target->Mass;
  double mc = Compound->Mass;
  double pb, Eb;
  FourVector Ptot;

  cout.width(5);
  cout << "Stp";
  cout.width(15);
  cout << "ExMax[MeV]";
  cout.width(15);
  cout << "ExMin[MeV]";
  cout.width(15);
  cout << "Exmax[MeV]";
  cout.width(15);
  cout << "Exmin[MeV]";
  cout.width(15);
  cout << "ECMmax[MeV]";
  cout.width(15);
  cout << "ECMmin[MeV]" << endl;

  for (int stp=0; stp<AnodeRows+1; stp++) {
    // Exc. energy of compound particle.
    cout.width(5);
    cout << stp;
    // Upper limit
    // Linear momentum and total energy of the beam particle in the lab.
    pb = sqrt(2*mb*Kb*(1 + Kb/2/mb));
    Eb = sqrt(mb*mb + pb*pb);
    // Total four momentum in the lab.
    Ptot.SetCoords(Eb+mt, 0, 0, pb);
    cout.width(15);
    cout.precision(4);
    cout << sqrt(Ptot*Ptot) - mc; // excitation energy in compound nucleus
    // Initial momentum in the CM reference frame (Ptot*Ptot is an invariant quantity).
    double pCM_max = sqrt((Ptot*Ptot-pow(mt+mb,2))*(Ptot*Ptot-pow(mb-mt,2))/(4*(Ptot*Ptot)));    
    // Lower limit
    Kb -= DeltaEB[stp][AnodeCols];
    if (Kb>0) {
      // Linear momentum and total energy of the beam particle in the lab.
      pb = sqrt(2*mb*Kb*(1 + Kb/2/mb));
      Eb = sqrt(mb*mb + pb*pb);
      // Total four momentum in the lab.
      Ptot.SetCoords(Eb+mt, 0, 0, pb);
      cout.width(15);
      cout.precision(4);
      cout << sqrt(Ptot*Ptot) - mc;
      // Non-rel center-of-mass energy
      double Ecm = mt*Kb/(mb+mt);
      // Energy of target particle (in normal kinematics)
      double Kt = Ecm*(mt+mb)/mb;
      cout.width(15);
      cout.precision(4);
      cout << Kt;
      cout.width(15);
      cout.precision(4);
      cout << sqrt(Ptot*Ptot) - mc;
      // Initial momentum in the CM reference frame (Ptot*Ptot is an invariant quantity).
      double pCM_min = sqrt((Ptot*Ptot-pow(mt+mb,2))*(Ptot*Ptot-pow(mb-mt,2))/(4*(Ptot*Ptot)));
      cout.width(15);
      cout.precision(4);
      cout << sqrt(mt*mt+pCM_max*pCM_max) - mt;
      cout.width(15);
      cout.precision(4);
      cout << sqrt(mt*mt+pCM_min*pCM_min) - mt;
      cout.width(15);
      cout.precision(4);
      //    cout << sqrt(mb*mb+pCM_max*pCM_max) - mb;
      cout << Kb*mt/(mb+mt) << endl;
      // cout.width(15);
      // cout.precision(4);
      // cout << sqrt(mb*mb+pCM_min*pCM_min) - mb << endl;
    }
    else {
      cout.width(15);
      cout << "0";
      cout.width(15);
      cout << "0";
      cout.width(15);
      cout << "0" << endl;
    }
  }
#endif
    return;
  }


///////////////////////////////////////////////////////////////////////////////////
// Private method, to be used in an event loop.
// This method transports the particle from its initial position until the maximum
// transport time has been reached or if the particle stops or if it leaves the 
// detector volume. As the particle is propagated through the gaseous medium, 
// energy loses in each anode segment are saved in a bidimensional array (DE) which 
// is retured when this method ends. The last column element (AnodeCols) is
// reserved for the total energy loss in all other columns for that strip.
///////////////////////////////////////////////////////////////////////////////////
int MUSIC_Simulator::PropagateParticle(Particle* PO, int Event, double MaxTime, double UserStep, double** DE)
{
  for (int stp = 0; stp<AnodeRows; stp++) 
    for (int col = 0; col<AnodeCols+1; col++) 
      DE[stp][col] = 0;

  if (PrintLevel>0) {
    Log << "\n*******************************************************************"
	<< "\nmusicsim::PropagateParticle ***************************************\n"
	<< PO->Name << endl;
  }
  if (PO->DoNotPropagate) {
    if (PrintLevel>0) 
      Log << "Not propagating!" << endl;
    return 0;
  }
  
#if 1
  TGeoVolume* Vol = 0;
  bool Skip = 0;
  int step = 0;
  // Get the initial conditions from the Particle object.
  double m = PO->Mass;
  double tf=0, xf=0, yf=0, zf=0;
  double ti, xi, yi, zi;
  PO->GetX(ti, xi, yi, zi);
  double Ene, px, py, pz;
  PO->GetP(Ene, px, py, pz);
  double p_mag = sqrt(px*px + py*py + pz*pz);
  double Ki = PO->GetKE();
  
  PO->ResetTrace();
  double tt, xt, yt, zt, Kt; // initial coordinates for the particle trajectory
  tt = ti;
  xt = xi;
  yt = yi;
  zt = zi;
  Kt = Ki;
  PO->SetTracePoint((float)tt, (float)xt, (float)yt, (float)zt, (float)Kt);
  //  PO->AllTraj[Event]->SetName(Form("%s evt %d", PO->Name.c_str(), Event));
  PO->Trajectory->RemoveElements();
  PO->Trajectory->SetName(Form("%s evt %d", PO->Name.c_str(), Event));
  
  if (PrintLevel>0)
    Log << "MaxTime=" << MaxTime << " ns\n" << "Initial time ti=" << ti << " ns"<< endl;
  
  // Get azimuthal and polar angles (direction of the particle)
  double phi = PO->GetPhi();
  double theta = PO->GetTheta();
  
  // The particle trajectory is calculated numerically with very fine
  // time steps.
  while (ti<MaxTime) {
    // For the current point with given components of the momentum,
    // compute the velocity in cm/ns. Although the formula in SI is
    // v=c^2p/E, in my units ([p]=MeV/c, [E]=MeV and c=30cm/ns), this
    // becomes v=cp/E. The formula is correct!
    px = p_mag*cos(phi)*sin(theta);
    py = p_mag*sin(phi)*sin(theta);
    pz = p_mag*cos(theta);
    double vx = c*px/Ene;
    double vy = c*py/Ene;
    double vz = c*pz/Ene;
    double vel = sqrt(vx*vx+vy*vy+vz*vz);
    
    // Use short steps for dense media
    double Dt = 0;
    if (vel>0) 
      Dt = UserStep/vel;   // from p/c = m*dx/dt -> dt = c*m*dx/p
    
    tf = ti + Dt;
    xf = xi + vx*Dt;
    yf = yi + vy*Dt;
    zf = zi + vz*Dt;
    
    if (PrintLevel>0)
      if (PrintLevel>1 || step==0) {
	Log << "step " << step << " (";
	Log.precision(7);
	Log << tf << "ns, " << xf << "cm, " << yf << "cm, " << zf << "cm)  Dt=" << Dt << "ns" << endl;
      } 
    
    // Exit the while loop if the particle leaves the anode volume.
    if (zf>AnodeDepth || zf<0 || xf>AnodeLength/2 || xf<-AnodeLength/2 || 
	yf>AnodeHeight/2 || yf<-AnodeHeight/2) {
      if (PrintLevel>0)
	Log << "Particle reached end of active volume." << endl;
      break;
    }
    
    // Find the anode segment in which the particle is moving.
    int stp = -1;
    int col = -1;
    Vol = Geo->FindNode(xf, yf, zf)->GetVolume(); 
    for (int s=0; s<AnodeRows; s++) 
      for (int c=0; c<AnodeCols; c++) 
	if (Vol==VolAnode[s][c]) {
	  stp = s;
	  col = c;
	  break;
	} 
    // Exit the while loop if the anode segment was not found.
    if (stp==-1 || col==-1) {
      if (PrintLevel>0)
	Log << "WARNING: anode segment not found." << endl;
      break;
    }

    // Get the energy loss for the initial energy (Ki) in the gas
    // medium (0) over a small (differential) path length (dist).
    double dist = sqrt(pow(xf-xi,2) + pow(yf-yi,2) + pow(zf-zi,2));
    double Kf = 0;
    if (UserStep>0) 
      Kf = PO->GetFinalEnergy(0, Ki, dist, dist/10);
    else
      Kf = PO->GetInitialEnergy(0, Ki, dist, dist/10);

    // Exit the while loop if the particle has stopped.
    if (Kf<0.001) {
      if (PrintLevel>0) {
	Log << "Particle stops inside active volume (Kf<1 keV)." << endl;
	if (Kf<0)	
	  Log << "WARNING: Less than ZERO K.E.! " << Kf << endl;
      }    
      break;
    }
    
    DE[stp][col] += fabs(Ki - Kf);
    // DE[stp][AnodeCols] += fabs(Ki - Kf); // AnodeCols component determined at the end
    
    if (PrintLevel>0) 
      if (PrintLevel>1 || step==0) 
	Log << "d=" << dist << " cm" << " \tKf=" << Kf << " \tDE=" << DE[stp][col] << endl; 
    
    // Get the momentum magnitude using the new (lower) kinetic energy
    // (I'm using a relativistic formula, although it is unlikely for
    // us to use relativistic energies).
    p_mag = sqrt(2*m*Kf*(1+Kf/2/m));
    // Reduce (or increase) the total energy by amount of energy
    // deposited in the medium.
    if (UserStep>0)
      Ene -= fabs(Ki - Kf);
    else
      Ene += fabs(Ki - Kf);
    
    // Move the coordinates of the initial point to the ones of the
    // final point.
    ti = tf;
    xi = xf;
    yi = yf;
    zi = zf;
    Ki = Kf;
    step++;
    // Add a new point to the trace and trajectory if enough time has
    // passed.
    if (tf-tt>UserStep) {   
      // if (PO->SaveTrajectory) {
      // 	if (!Skip) {
      // 	  //	  PO->AllTraj[Event]->AddLine(xt,yt,zt, xf,yf,zf);
      // 	  PO->Trajectory->AddLine(xt,yt,zt, xf,yf,zf);
      // 	  Skip = 1;
      // 	}
      // 	else 
      // 	  Skip = 0;
      // }
      tt = tf;
      xt = xf;
      yt = yf;
      zt = zf;
      Kt = Kf;
      PO->SetTracePoint((float)tt, (float)xt, (float)yt, (float)zt, (float)Kt);
    }
  } // end 


  for (int stp = 0; stp<AnodeRows; stp++) 
    for (int col = 0; col<AnodeCols; col++)  {
      if (DE[stp][col]>0) {
	// if (ctf.Eres>0.0) {
	//   // Adding randomness to the energy loss to mimic eperimental jitter
	//   Gaussian->SetParameters(1.0, 0.0, ctf.Eres);
	//   DE[stp][col] += Gaussian->GetRandom();
	// }	
	DE[stp][AnodeCols] += DE[stp][col];
      }
    }
  

  PO->SetX(tf,xf,yf,zf);
  PO->SetP(Ene, p_mag*cos(phi)*sin(theta), p_mag*sin(phi)*sin(theta), p_mag*cos(theta));
#endif
  return 1;
}


/////////////////////////////////////////////////////////////////////////////////////////////////
// Reset the values of the branches in the TTree
/////////////////////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::ResetBranches()
{

  for (int stp=0; stp<ExpAnodeStps; stp++) {
    de_r[stp] = 0;
    de_l[stp] = 0;
    seg[stp] = -1;
  }
  strip0 = 0;
  strip17 = 0;
  cathode = 0;
  reacStp = -1;
  Kb = 0;
  for (int er=0; er<MaxEva; er++) {
    Kl[er] = 0;
    Kh[er] = 0;
    phi_CM[er] = -1;
    theta_CM[er] = -1;
    phi_l[er] = -1;
    theta_l[er] = -1;
    phi_h[er] = -1;
    theta_h[er] = -1;
    xfl[er] = 0;
    yfl[er] = 0;
    zfl[er] = -1000;
  }
  xr = yr = 0;
  zr = -1000;
  xfe = yfe = 0;
  zfe = -1000;
}



/////////////////////////////////////////////////////////////////////////////////////////////////
// This method intents to substitute the use of ROOT scripts.
/////////////////////////////////////////////////////////////////////////////////////////////////
int MUSIC_Simulator::run()
{
  //SetROOTSystemPointer(gSystem);
  // TEveManager for drawing 3D particle trajectories
  Eve = new TEveManager(960, 1018, kTRUE, "V");
  Eve->GetDefaultGLViewer()->SetClearColor(kWhite);
  // 3D particle tracks
  TrackBeam = new TEveArrow();
  TrackEvaP = new TEveArrow*[MaxEva];
  TrackEvaR = new TEveArrow*[MaxEva];
  for (int er=0; er<MaxEva; er++) {
    TrackEvaP[er] = new TEveArrow();
    TrackEvaR[er] = new TEveArrow();
  }
  
  // Canvas and legends for traces
  TraceCan = new TCanvas("TraceCan","Traces", 0, 0, 960, 1018);
  TraceCan->Divide(2,2);
  TraceCan->cd(1)->SetGrid();
  TraceCan->cd(2)->SetGrid();
  TraceCan->cd(3)->SetGrid();
  TraceCan->cd(4)->SetGrid();
  LegCol = new TLegend(0.692,0.616,0.826,0.861);
  LegPart = new TLegend(0.692,0.616,0.826,0.861);
  LabelKine = new TPaveText(0.152,0.679,0.437,0.875,"NDC");
  
  
  SetPrintLevel(ctf.PrintOpt);
  SetStripEnergyResolution(ctf.Eres);
  // Geometry
  SetAnode(ctf.AnodeGeom, 90, ctf.ELossBins, ctf.MaxELoss);

  // Beam
  SetBeamParticle(ctf.beamName, kBlack, ctf.SRIMbeam, ctf.dEdxScaleBeam);
  // Target
  SetTargetParticle(ctf.target);
  // Compound particle
  SetCompoundParticle(ctf.compound);
  // Evaporation residues and particles
  for (int i=0; i<ctf.NumEvapPart; i++)
    SetEvapResAndPart(ctf.res[i],
		      ctf.SRIMres[i],
		      ctf.colorRes[i],
		      ctf.evap[i],
		      ctf.SRIMevap[i],
		      ctf.colorEvap[i],
		      ctf.dEdxScaleRes[i],
		      ctf.dEdxScaleEvap[i]);
  
  if (ctf.Method==0) {
    // Simulate events for one strip or generate trace data base (see below)
    Simulate(ctf.strip,
	     ctf.NEvents,
	     ctf.MaxTime,
	     ctf.SimStep,
	     ctf.Update,
	     ctf.Wait,
	     ctf.FileName,
	     ctf.FileOpt);
  }
  else if (ctf.Method==1) {
    cout << "musicsim warning: GenerateTraceDataBase method not ready." << endl;
    // GenerateTraceDatabase("TraceDB.root", 
    // 			  ThCMMin, ThCMMax, ThSteps, 
    // 			  PhiCMMin, PhiCMMax, PhiSteps,
    // 			  MaxTime, SimStep, Update, Wait);
  }
  
  return 0;
}




///////////////////////////////////////////////////////////////////////////////////
// Establish the dimensions of the MUSIC anode segments by reading a text file with
// the following structre
//------------------------------------------------------------------------------|
// Title                                                                        |
// Number of strips: #                                                          |
// Number of columns: #                                                         |
// Stp  Col  StpID  Name    Dx   Dy   Dz  Color	 Comment                        |
// #    #    #	    string  #f   #f   #f  #	 string (spaces allowed)        |
// ... (until all segments have been specified)                                 |
//------------------------------------------------------------------------------|
// , where '#' is meant to be replaced by an integer and '#f' by a floating point 
// number.
///////////////////////////////////////////////////////////////////////////////////
 void MUSIC_Simulator::SetAnode(string AnodeGeomFile, short Trans, int ELossBins, float MaxELoss)
{
#if 1
  ifstream GeomFile;
  string line;
  int line_counter = 0;

  // It's required for column description in the AnodeGeomFile to
  // match the string below.
  string ColDescription = "Stp	Col	StpID	Name	Dx	Dy	Dz	Color	Comment";

  AnodeDepth = 0;
  AnodeLength = 0;
  AnodeHeight = 0;
  AnodeRows = 0;
  AnodeCols = 0;

  if (PrintLevel>0) {
    Log << "\n*******************************************************************"
	<< "\nmusicsim::SetAnode ************************************************" << endl;
  }
    
  
  GeomFile.open(AnodeGeomFile.c_str());
  if (!GeomFile.is_open()) 
    cout << "ERROR: Anode geometry file \"" << AnodeGeomFile << "\" couldn't be opened." << endl;
  else {
    cout << "Loading anode geometry loaded from \"" << AnodeGeomFile << "\"." << endl;

    // Get the number of strips and columns for the anode segments
    do {
      getline(GeomFile, line);
      if (line.find("strips:")!=string::npos) 
	AnodeRows = atoi(line.substr(line.find(':')+1).c_str());
      if (line.find("columns:")!=string::npos)
	AnodeCols = atoi(line.substr(line.find(':')+1).c_str());     
    } while (!GeomFile.eof());
    GeomFile.close();
    
    cout << "Anode strips: " << AnodeRows << endl;
    cout << "Anode columns: " << AnodeCols << endl;    
    if (AnodeRows>0 && AnodeCols>0) {
      // Initialize the anode segment dimensions
      AnodeDX = new double*[AnodeRows];
      AnodeDY = new double*[AnodeRows];
      AnodeDZ = new double*[AnodeRows];
      AnodeColor = new short*[AnodeRows];
      AnodeSegName = new string*[AnodeRows];
      AnodeStpID = new int*[AnodeRows];
      for (int stp=0; stp<AnodeRows; stp++) {
	AnodeDX[stp] = new double[AnodeCols];
	AnodeDY[stp] = new double[AnodeCols];
	AnodeDZ[stp] = new double[AnodeCols];
	AnodeColor[stp] = new short[AnodeCols];
	AnodeSegName[stp] = new string[AnodeCols];
	AnodeStpID[stp] = new int[AnodeCols];
	for (int col=0; col<AnodeCols; col++) {
	  AnodeDX[stp][col] = 0;
	  AnodeDY[stp][col] = 0;
	  AnodeDZ[stp][col] = 0;
	  AnodeColor[stp][col] = kWhite; // kWhite = (const enum EColor) 0
	  AnodeSegName[stp][col] = "";
	  AnodeStpID[stp][col] = -1;
	}
      }
 
      // Get the name and dimensions for each anode segments(strip and
      // column)
      GeomFile.open(AnodeGeomFile.c_str());
      do {
	getline(GeomFile, line);
	if (line.find(ColDescription)<line.npos)
	  break;
      } while (!GeomFile.eof());
      // First we need to count how many (not empty) lines are listed
      do {
	getline(GeomFile, line);
	if (!line.empty())
	  line_counter++;
      } while (!GeomFile.eof());
      GeomFile.close();
      if (PrintLevel>1)
	Log << "Total lines: " << line_counter << endl;
      // Now that the number of lines has been established reopen the
      // text file to load the parameters
      GeomFile.open(AnodeGeomFile.c_str());
      do {
	getline(GeomFile, line);
	if (line.find(ColDescription)<line.npos)
	  break;
      } while (!GeomFile.eof());
      // Loading parameters and printing them to confirm they have been
      // read correctly
      for (int nl=0; nl<line_counter; nl++) {
	int stp = 0;
	int col = 0;
	int id = -1;    
	string name;
	double dx, dy, dz;
	short color;
	GeomFile >> stp >> col >> id >> name >> dx >> dy >> dz >> color;
	getline(GeomFile, line); // The last column is for comments
	AnodeStpID[stp][col] = id;
	AnodeSegName[stp][col] = name;
	AnodeDX[stp][col] = dx;
	AnodeDY[stp][col] = dy;
	AnodeDZ[stp][col] = dz;
	AnodeColor[stp][col] = color;
	if (PrintLevel>0)
	  Log << nl << "\t" << AnodeStpID[stp][col] << "\t" << AnodeSegName[stp][col] << "\t" 
	      << AnodeDX[stp][col] << "\t" << AnodeDY[stp][col] << "\t" << AnodeDZ[stp][col] 
	      << "\t" << AnodeColor[stp][col] << endl;
      }
      GeomFile.close();
    
      // The total anode depth (distance along the z axis) is the sum of
      // all segment depths for the first column
      for (int stp=0; stp<AnodeRows; stp++)
	AnodeDepth += AnodeDZ[stp][0];

      // The total anode length (distance along the x axis) is the sum
      // of all segment lengths for the first strip
      for (int col=0; col<AnodeCols; col++)
	AnodeLength += AnodeDX[0][col];

      // The total anode height (distance along the y axis) is the hight
      // of the segment in the first strip and first column
      AnodeHeight = AnodeDY[0][0];

      // Define anode volumes. The Vacuum medium used below is just
      // because a medium is needed when defining a Volume, it has
      // nothing to do with the energy loss. The actual stopping power
      // tables are contained in the Particle objects.
      double z0 = 0;
      VolAnode = new TGeoVolume**[AnodeRows];
      for (int stp=0; stp<AnodeRows; stp++) {
	VolAnode[stp] = new TGeoVolume*[AnodeCols];
	z0 += AnodeDZ[stp][0]/2;
	double x0 = -AnodeLength/2;
	for (int col=0; col<AnodeCols; col++) {
	  if (AnodeDX[stp][col]>0) {
	    VolAnode[stp][col] = Geo->MakeBox(Form("VolAnode%d%d",stp,col),
					      Vacuum /*just because a medium is need*/,
					      AnodeDX[stp][col]/2, AnodeDY[stp][col]/2,
					      AnodeDZ[stp][col]/2);
	    VolAnode[stp][col]->SetLineColor(AnodeColor[stp][col]);
	    VolAnode[stp][col]->SetTransparency(Trans);
	    x0 += AnodeDX[stp][col]/2;
	    VolTop->AddNode(VolAnode[stp][col], 1, new TGeoTranslation(x0,0,z0));
	    x0 += AnodeDX[stp][col]/2;
	  }
	  else 
	    VolAnode[stp][col] = 0;
	}
	z0 += AnodeDZ[stp][0]/2;
      }

      // Arrays where the energy loss in each strip will be saved. The
      // last column (AnodeCols) is reserved for the summed energy
      // loss in the other columns for that strip.
      DeltaEB_ave = new double*[AnodeRows];// average beam energy loss
      DeltaEB = new double*[AnodeRows];    // beam
      DeltaEL = new double*[AnodeRows];    // light
      DeltaEH = new double*[AnodeRows];    // heavy
      DeltaED1 = new double*[AnodeRows];   // decay daughter1
      DeltaED2 = new double*[AnodeRows];   // decay daughter2
      for (int stp=0; stp<AnodeRows; stp++) {
	DeltaEB_ave[stp] = new double[AnodeCols+1];
	DeltaEB[stp] = new double[AnodeCols+1];
	DeltaEL[stp] = new double[AnodeCols+1];
	DeltaEH[stp] = new double[AnodeCols+1];
	DeltaED1[stp] = new double[AnodeCols+1];
	DeltaED2[stp] = new double[AnodeCols+1];
	for (int col=0; col<AnodeCols+1; col++) {
	  DeltaEB_ave[stp][col] = 0;
	  DeltaEB[stp][col] = 0;
	  DeltaEL[stp][col] = 0;
	  DeltaEH[stp][col] = 0;
	  DeltaED1[stp][col] = 0;
	  DeltaED2[stp][col] = 0;
	}
      }
      // New stuff for evaporation residues
      DeltaE_EvaR = new double**[MaxEva];
      DeltaE_EvaP = new double**[MaxEva];
      for (int er=0; er<MaxEva; er++) {	
	DeltaE_EvaR[er] = new double*[AnodeRows];
	DeltaE_EvaP[er] = new double*[AnodeRows];
	for (int stp=0; stp<AnodeRows; stp++) {	  
	  DeltaE_EvaR[er][stp] = new double[AnodeCols+1];
	  DeltaE_EvaP[er][stp] = new double[AnodeCols+1];
	  for (int col=0; col<AnodeCols+1; col++) {
	    DeltaE_EvaR[er][stp][col] = 0;
	    DeltaE_EvaP[er][stp][col] = 0;
	  }
	}
      }

    } // end if (AnodeRows>0 && AnodeCols>0)
    
    HCTB = new TH2F("HCTB","Beam", AnodeRows,-0.5, AnodeRows-0.5, ELossBins,0,MaxELoss);
    HCTB->GetXaxis()->SetTitle("Strip index (in AnodeGeometry file)");
    HCTB->GetXaxis()->CenterTitle();
    HCTB->GetYaxis()->SetTitle("Energy loss [MeV]");
    HCTB->GetYaxis()->CenterTitle(); 
    
    HCT = new TH2F("HCT","Column traces", AnodeRows,-0.5, AnodeRows-0.5, ELossBins,0,MaxELoss);
    HCT->GetXaxis()->SetTitle("Strip index (in AnodeGeometry file)");
    HCT->GetXaxis()->CenterTitle();
    HCT->GetYaxis()->SetTitle("Energy loss [MeV]");
    HCT->GetYaxis()->CenterTitle(); 
    
    HPT = new TH2F("HPT","Particle traces", AnodeRows,-0.5, AnodeRows-0.5, ELossBins,0,MaxELoss);
    HPT->GetXaxis()->SetTitle("Strip index (in AnodeGeometry file)");
    HPT->GetXaxis()->CenterTitle();
    HPT->GetYaxis()->SetTitle("Energy loss [MeV]");
    HPT->GetYaxis()->CenterTitle(); 
    
  }
  
  DrawMUSIC(Eve, 85);

  cout << "Anode dimensions: " << AnodeLength << "x" << AnodeHeight << "x" 
       << AnodeDepth << "cm^3" << endl;
#endif
  return;
}

///////////////////////////////////////////////////////////////////////////////////
// Create the beam particle object.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetBeamParticle(string Name, int Color, string ELossFile, float dEdxScale)
{
#if 1
  // Use the nuclide finder object to determine the mass and atomic
  // number of this particle. Mass must be in MeV/c^2) and Z (in e).
  double m = NuF->GetMass(Name, "MeV/c^2");
  int Z = NuF->GetZ(Name);
  Beam = new Particle(Name, m, Z, false /*SaveTrajectories off*/);
  Beam->SetTrajectoryAtt((short)Color);
  // Currently, this simulation is restricted to one medium (gas) in MUSIC.
  Beam->SetMedium(ELossFile, dEdxScale);
  TrackBeam->SetName(Name.c_str());
  TrackBeam->SetMainColor(Color);
  TrackBeam->SetPickable(kTRUE);
  Eve->AddElement(TrackBeam);
  //  Kb_after_window = K;
#endif
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Create the fused (a.k.a. compound) particle object.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetCompoundParticle(string Name)
{
#if 1
  // Use the nuclide finder object to determine the mass and atomic
  // number of this particle. Mass must be in MeV/c^2) and Z (in e).
  double m = NuF->GetMass(Name, "MeV/c^2");
  int Z = NuF->GetZ(Name);
  Compound = new Particle(Name, m, Z, false /*SaveTrajectories off*/);
  // If both beam and target particles have been set, print the
  // maximum compund excitation energy, ExMax
  if (Beam && Target) {
    double mb = Beam->Mass;
    double Kb = ctf.Kb;
    double pb = sqrt(2*mb*Kb*(1 + Kb/(2*mb)));
    double Eb = sqrt(mb*mb + pb*pb);
    FourVector Pb("Pb", Eb, 0, 0, pb);
    FourVector Pt("Pt", Target->Mass, 0, 0, 0);
    FourVector Ptot("Total four-mom. in the lab");
    Ptot = Pb + Pt;
    double ExMax = sqrt(Ptot*Ptot) - Compound->Mass;
    cout << "Maximum excitation energy of " << Name << " (compound) = " << ExMax << " MeV" << endl;
  }
#endif
  return;
}

///////////////////////////////////////////////////////////////////////////////////
// Create the decay daughter1 of the heavy particle.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetDecayDaughter1(string Name, int Color, string ELossFile)
{
#if 1
  // Use the nuclide finder object to determine the mass and atomic
  // number of this particle. Mass must be in MeV/c^2) and Z (in e).
  double m = NuF->GetMass(Name, "MeV/c^2");
  int Z = NuF->GetZ(Name);
  DeDau1 = new Particle(Name, m, Z, false /*SaveTrajectories off*/);
  DeDau1->SetTrajectoryAtt((short)Color);
  // Currently, this simulation is restricted to one medium (gas) in MUSIC.
  DeDau1->SetMedium(ELossFile);
#endif
  return;
}

///////////////////////////////////////////////////////////////////////////////////
// Create the decay daughter2 of the heavy particle.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetDecayDaughter2(string Name, int Color, string ELossFile)
{
  // Use the nuclide finder object to determine the mass and atomic
  // number of this particle. Mass must be in MeV/c^2) and Z (in e).
  double m = NuF->GetMass(Name, "MeV/c^2");
  int Z = NuF->GetZ(Name);
  DeDau2 = new Particle(Name, m, Z, false /*SaveTrajectories off*/);
  DeDau2->SetTrajectoryAtt((short)Color);
  // Currently, this simulation is restricted to one medium (gas) in MUSIC.
  DeDau2->SetMedium(ELossFile);
  return;
}

///////////////////////////////////////////////////////////////////////////////////
// Create evaporation residue and particle objects.
// Use the nuclide finder object to determine the mass and atomic
// number of this particle. Mass must be in MeV/c^2 and Z in e.
// Currently particles restricted to one medium (gas).
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetEvapResAndPart(string ResName, string ResELossFile, int ResColor, 
					string ParName,	string ParELossFile, int ParColor,
					float dEdxScaleRes, float dEdxScalePar)
{
  if (CurEva>=MaxEva) {
    cout << "Warning: No more than " << MaxEva << " evaporation particles allowed." << endl;  
    return;
  }
  
  // Evaporated particle (p,n,4He)
  double mp = NuF->GetMass(ParName, "MeV/c^2");
  int Zp = NuF->GetZ(ParName);
  ParName += std::to_string(CurEva);
  EvaP[CurEva] = new Particle(ParName, mp, Zp, false /*SaveTrajectories off*/);
  EvaP[CurEva]->SetTrajectoryAtt((short)ParColor);
  EvaP[CurEva]->SetMedium(ParELossFile, dEdxScalePar);
  EvaP[CurEva]->Print();
  TrackEvaP[CurEva]->SetName(ParName.c_str());
  TrackEvaP[CurEva]->SetMainColor(ParColor);
  TrackEvaP[CurEva]->SetPickable(kTRUE);
  Eve->AddElement(TrackEvaP[CurEva]);

  // Evaporation residue (heavy particle)
  double mr = NuF->GetMass(ResName, "MeV/c^2");
  int Zr = NuF->GetZ(ResName);
  EvaR[CurEva] = new Particle(ResName, mr, Zr, false /*SaveTrajectories off*/);
  EvaR[CurEva]->SetTrajectoryAtt((short)ResColor);
  EvaR[CurEva]->SetMedium(ResELossFile, dEdxScaleRes);
  EvaR[CurEva]->Print();  
  TrackEvaR[CurEva]->SetName(ResName.c_str());
  TrackEvaR[CurEva]->SetMainColor(ResColor);
  TrackEvaR[CurEva]->SetPickable(kTRUE);
  Eve->AddElement(TrackEvaR[CurEva]);

  CurEva++;
  return;
}
 
 
///////////////////////////////////////////////////////////////////////////////////
// Create the heavy particle (evaporation residue) object.
///////////////////////////////////////////////////////////////////////////////////
 void MUSIC_Simulator::SetHeavyParticle(string Name, int Color, string ELossFile, int NEexc,
    double* Eexc)
{
#if 1
  // Use the nuclide finder object to determine the mass and atomic
  // number of this particle. Mass must be in MeV/c^2) and Z (in e).
  double m = NuF->GetMass(Name, "MeV/c^2");
  int Z = NuF->GetZ(Name);
  Heavy = new Particle(Name, m, Z, false /*SaveTrajectories off*/);
  Heavy->SetTrajectoryAtt((short)Color);
  // Currently, this simulation is restricted to one medium (gas) in MUSIC.
  Heavy->SetMedium(ELossFile);
  Heavy->SetExcEnergies(NEexc, Eexc);
#endif
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Private method, to be used in an event loop.  Establishes the kinematics of the
// particles right after the window.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetInitialKinematics(double Kbi/*MeV*/)
{
#if 1
  // Mass of the beam particle.
  double mb = Beam->Mass;

  // Linear momentum and total energy of the beam particle in the lab.
  double pb = sqrt(2*mb*Kbi*(1 + Kbi/2/mb));
  double theta_b = 0;
  double phi_b = 0;
  double Eb = sqrt(mb*mb + pb*pb);

  // Four-momentum of the beam.
  Beam->SetP(Eb, pb*sin(theta_b)*cos(phi_b), pb*sin(theta_b)*sin(phi_b), pb*cos(theta_b));
  Beam->SetX(0, 0, 0, 0);
  if (PrintLevel>0) {
    Log << "\n*******************************************************************" << endl;
    Log << "musicsim::SetInitialKinematics ************************************" << endl;    
    Beam->Print(Log);
  }
#endif
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Create the light particle object. It is assumed that the light particle 
// (e.g. p, n, alpha) is in its ground state.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetLightParticle(string Name, int Color, string ELossFile)
{
#if 1
  // Use the nuclide finder object to determine the mass and atomic
  // number of this particle. Mass must be in MeV/c^2) and Z (in e).
  double m = NuF->GetMass(Name, "MeV/c^2");
  int Z = NuF->GetZ(Name);
  Light = new Particle(Name, m, Z, false /*SaveTrajectories off*/);
  Light->SetTrajectoryAtt((short)Color);
  // Currently, this simulation is restricted to one medium (gas) in MUSIC.
  Light->SetMedium(ELossFile);
#endif
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Select the amount of information displayed by the Simulator.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetPrintLevel(int Level/*0-2*/)
{
#if 1
  if (Level<0 || Level>2) {
    cout << "Warning: Invalid information printing level.\n"
	 << "Valid options are:\n"
	 << "\t0 - minimum printing\n"
	 << "\t1 - info per event (e.g. propagator initial and final conditions)"
	 << "\t2 - print all" << endl;
    PrintLevel = 0;
  }
  else {
    cout << "See musicsim.log file for detailed information" << endl;
    PrintLevel = Level;
    Log.open("musicsim.log");
    Log << "================================================================================" << endl;
    Log << "|--- MUSIC simulator log file -------------------------------------------------|" << endl;
  }
#endif
  return;
}


///////////////////////////////////////////////////////////////////////////////////
//   _____      _     _____                 _   _                                //
//  / ____|    | |   |  __ \               | | (_)                               //
// | (___   ___| |_  | |__) |___  __ _  ___| |_ _  ___  _ __                     //
//  \___ \ / _ \ __| |  _  // _ \/ _` |/ __| __| |/ _ \| '_ \                    //
//  ____) |  __/ |_  | | \ \  __/ (_| | (__| |_| | (_) | | | |                   //
// |_____/ \___|\__| |_|  \_\___|\__,_|\___|\__|_|\___/|_| |_|                   //
//  _  ___                            _   _                                      //
// | |/ (_)                          | | (_)                                     //
// | ' / _ _ __   ___ _ __ ___   __ _| |_ _  ___ ___                             //
// |  < | | '_ \ / _ \ '_ ` _ \ / _` | __| |/ __/ __|                            //
// | . \| | | | |  __/ | | | | | (_| | |_| | (__\__ \                            //
// |_|\_\_|_| |_|\___|_| |_| |_|\__,_|\__|_|\___|___/                            //
//                                                                               //
// Private method, to be used in an event loop. Establish the kinematics of the  //
// particles at the reaction point. If theta_CM and phi_CM are both -1 (default  //
// values) then their values are assigned randomly.                              //
///////////////////////////////////////////////////////////////////////////////////
int MUSIC_Simulator::SetReactionKinematics(double Kbr/*MeV*/, double zr/*cm*/, double tof/*ns*/,
					   double theta_CM, double phi_CM)
{
  int ReactionAllowed = 1;
#if 1
  if (PrintLevel>0) {
    Log << "\n*******************************************************************" << endl;
    Log << "musicsim::SetReactionKinematics ***********************************" << endl;
  }
  double mb = Beam->Mass;
  double mt = Target->Mass;
  double mc = Compound->Mass;

  // Linear momentum and total energy of the beam particle in the lab.
  double pb = sqrt(2*mb*Kbr*(1 + Kbr/2/mb));
  double theta_b = 0;
  double phi_b = 0;
  double Eb = sqrt(mb*mb + pb*pb);

  // Center of mass beta (v/c)
  double BetaX = pb*sin(theta_b)*cos(phi_b)/(Eb+mt);
  double BetaY = pb*sin(theta_b)*sin(phi_b)/(Eb+mt);
  double BetaZ = pb*cos(theta_b)/(Eb+mt);
  if (PrintLevel>0) {
    Log << "Center-of-mass velocity (v/c):" << endl;
    Log << "BetaX=" << BetaX << "  BetaY=" << BetaY << "  BetaZ=" << BetaZ << endl; 
    Log << "--- Beam and target particles --------------------------------------" << endl;
  }

  // Four-momentum of the beam.
  Beam->SetP(Eb, pb*sin(theta_b)*cos(phi_b), pb*sin(theta_b)*sin(phi_b), pb*cos(theta_b));
  Beam->SetX(tof, 0, 0, zr);
  if (PrintLevel>0)
    Beam->Print(Log);
  
  // Four-momentum of the target.
  Target->SetP(mt, 0, 0, 0);
  Target->SetX(tof, 0, 0, zr);
  if (PrintLevel>0)
    Target->Print(Log);
  
  // Total four momentum in the lab.
  FourVector Ptot = Beam->GetP() + Target->GetP();
  //   Ptot.SetName("Total four-mom. in the lab");
  
  // Exc. energy of compound particle.
  Compound->SetP(Ptot);
  Compound->SetX(tof, 0, 0, zr);
  Compound->SetExcEnergy(sqrt(Ptot*Ptot) - mc);
  if (PrintLevel>0) {
    Log << "--- Compound particle ----------------------------------------------" << endl;
    //    Log << "Eexc(" << Compound->Name << ") = " << sqrt(Ptot*Ptot) - mc << " MeV" << endl;
    Compound->Print(Log);
  }
  // Assume all the particles will be propagated
  for (int er=0; er<CurEva; er++) {
    EvaP[er]->DoNotPropagate = false;
    EvaR[er]->DoNotPropagate = false;
  }
  
  // If user specifies angles used them for the first reaction
  int user_angles = 1;
  
  // Loop over the evaporation residues (heavy) and particles (light particle always in g.s.)
  string reacstr = Beam->Name + "(" + Target->Name + "," + EvaP[0]->Name + ")" + EvaR[0]->Name; 
  for (int er=0; er<CurEva; er++) {
    double ml = EvaP[er]->Mass;
    double mh = EvaR[er]->Mass;
    double Qvalue = sqrt(Ptot*Ptot) - ml - mh;
    if (PrintLevel>0) {
      Log << "--- Reaction -------------------------------------------------------" << endl;
      if (er>0)
	reacstr = EvaR[er-1]->Name + "->" + EvaP[er]->Name + "+" + EvaR[er]->Name;      
      Log << "reac=" << er << ": " << reacstr << endl;
      Log << "Qvalue=" << Qvalue << " MeV" << endl;
      
    }
    if (Qvalue<0 /*&& er>0*/) {
      // The present reaction is not energetically allowed
      // Propagate the previuos evap res and exit the loop
      //EvaR[er-1]->DoNotPropagate = false;
      if (PrintLevel>0)
	Log << "Negative Qvalue!\nThe following particles will NOT be propagated:" << endl; 
      for (int i=er; i<CurEva; i++) {
	EvaP[i]->DoNotPropagate = true;
	EvaR[i]->DoNotPropagate = true;
	if (PrintLevel>0)
	  Log << i << " " << EvaP[i]->Name << ", " << EvaR[i]->Name << endl;
      }
      break;
    }
    // When the Qvalue is positive, assume that the reactio or decay
    // took place and stop the propagation of the previous evaporation
    // residue
    else {
      if (er>0)
	EvaR[er-1]->DoNotPropagate = true;
    }
    // if (er==CurEva-1)
    //   EvaR[er]->DoNotPropagate = false;
    
    double Ex = Rdm->Uniform(/*0.0*/Qvalue/2, Qvalue);
    //   Ex = 0; // Forcing g.s. of evaporation residue
#if 0
    // Warning: for 17F(alpha,p) only!!
    double ExIndex = Rdm->Uniform(0.0,3);
    if (ExIndex<1.0)
      Ex = 0;
    else if (ExIndex<2.0)
      Ex = 1.63;
    else
      Ex = 4.25;
#endif
    EvaR[er]->SetExcEnergy(Ex);
    
    // If the user did not specify the value of theta and phi (initial
    // theta=phi=-1) or for er>0, randomly select the scattering angle
    // in the center of mass.
    if ((theta_CM==-1 && phi_CM==-1) || er>0) {
      theta_CM = acos(Rdm->Uniform(-1.0,1.0));
      phi_CM = Rdm->Uniform(-pi,pi);
    }
    
    if (PrintLevel>0) {
      Log << "Ex(" << EvaR[er]->Name << ")=" << Ex << " MeV\ntheta_cm=" << theta_CM*180/pi
	  << "\nphi_cm=" << phi_CM*180/pi << endl;
      Log << "--- Outgoing particles (evap res = " << er << ") -------------------------------" 
	  << endl;
    }
    
    // Verify whether more particles can evaporate from the current residue or compound (Ptot)
    // In this case, the current light particle cannot evaporate from the previous reside/compound.
    // if (Ptot*Ptot<pow(ml+mh+Ex,2)) {
    //   EvaP[er]->SetP(ml, 0, 0, 0);
    //   EvaP[er]->SetX(tof, 0, 0, zr);
    //   EvaP[er]->DoNotPropagate = true;
    //   EvaR[er]->SetP(mh+Ex, 0, 0, 0);
    //   EvaR[er]->SetX(tof, 0, 0, zr);
    //   EvaR[er]->DoNotPropagate = true;
    //   break;
    // }
    // // In this case, the reaction IS energetically allowed.
    // else {
    // Final momentum in the CM reference frame (Ptot*Ptot is an invariant quantity).
    double pf_CM = sqrt((Ptot*Ptot-pow(ml+mh+Ex,2))*(Ptot*Ptot-pow(ml-mh-Ex,2))/(4*(Ptot*Ptot)));
      
    // Set the four-momentum components of the light particle in the center of mass.
    double plxCM = -pf_CM*sin(theta_CM)*cos(phi_CM);
    double plyCM = -pf_CM*sin(theta_CM)*sin(phi_CM);
    double plzCM = -pf_CM*cos(theta_CM);
    double ElCM = sqrt(ml*ml + pf_CM*pf_CM);
    EvaP[er]->SetP(ElCM, plxCM, plyCM, plzCM);

    // Set the four-momentum components of the heavy particle in the center of mass.
    double phxCM = pf_CM*sin(theta_CM)*cos(phi_CM);
    double phyCM = pf_CM*sin(theta_CM)*sin(phi_CM);
    double phzCM = pf_CM*cos(theta_CM);
    double EhCM = sqrt((mh+Ex)*(mh+Ex) + pf_CM*pf_CM);
    EvaR[er]->SetP(EhCM, phxCM, phyCM, phzCM);

    if (PrintLevel>0) {
      Log << "(((((((((( Before lorentz boost ))))))))))" << endl;
      EvaR[er]->Print(Log);
      EvaP[er]->Print(Log);
    }
      
    // Do a Lorentz transformation (boost) into the lab reference frame.
    // I've double checked the sign of the boost and is correct (-Beta).
    EvaP[er]->Boost(-BetaX, -BetaY, -BetaZ);
    EvaR[er]->Boost(-BetaX, -BetaY, -BetaZ); 

    // Now the total four-momentum (Ptot) is the 4-mom of the evap
    // res, so that subsequent "reactions" are just decays from the
    // present evaporation residue.
    Ptot = EvaR[er]->GetP();
      
    // Save the angles in degrees
    theta_l[er] = (EvaP[er]->GetTheta())*180/pi;
    phi_l[er] = (EvaP[er]->GetPhi())*180/pi;
    
    // Initial position of the light particle (at the target). This
    // particle will be propagated.
    EvaP[er]->SetX(tof, 0, 0, zr);
    //    EvaP[er]->DoNotPropagate = false;

    
    // Initial position of the evaporation residue (at the
    // target). We don't yet know if this particle will be
    // propagated (this will be known in the next cycle)
    EvaR[er]->SetX(tof, 0, 0, zr);
    // if (er>0)
    //   EvaR[er-1]->DoNotPropagate = true;
    
    // Save the angles in degrees
    theta_h[er] = (EvaR[er]->GetTheta())*180/pi;
    phi_h[er] = (EvaR[er]->GetPhi())*180/pi;


    // Print after lorentz boost
    if (PrintLevel>0) {
      Log << ")))))))))) After lorentz boost ((((((((((" << endl;
      EvaR[er]->Print(Log);
      EvaP[er]->Print(Log);
    }

    // Get the beta (v/c) for the evaporation residue. Will be used in the next 'er'.
    EvaR[er]->GetBeta(BetaX, BetaY, BetaZ);
    // if (PrintLevel>0) {
    //   Log << "Evap residue beta (v/c):" << endl;
    //   Log << "\tBetaX=" << BetaX << "  BetaY=" << BetaY << "  BetaZ=" << BetaZ << endl; 
    // }
  } // end for (er)

  
  // Fill the leaves related to the reaction kinematics
  Kb = Beam->GetKE();
  for (int er=0; er<CurEva; er++) {
    this->theta_CM[er] = theta_CM*180/pi;
    this->phi_CM[er] = phi_CM*180/pi;
    if (!EvaR[er]->DoNotPropagate)
      Kh[er] = EvaR[er]->GetKE();
    else
      Kh[er] = 0;
    if (!EvaP[er]->DoNotPropagate)
      Kl[er] = EvaP[er]->GetKE();
    else
      Kl[er] = 0;    
    theta_l[er] = (EvaP[er]->GetTheta())*180/pi;
    phi_l[er] = (EvaP[er]->GetPhi())*180/pi;
    theta_h[er] = (EvaR[er]->GetTheta())*180/pi;
    phi_h[er] = (EvaR[er]->GetPhi())*180/pi;
  }

  // Update the kinematics label and then draw it
  LabelKine->Clear();
  LabelKine->AddText("Kinematics");
  //  LabelKine->AddLine(0.0,0.76,1.0,0.76);
  LabelKine->AddText(Form("beam: K=%.2f MeV  z_{r}=%.2f cm  tof=%.1f ns", Kbr, zr, tof));
  for (int er=0; er<CurEva; er++) {
    if (EvaP[er] && !EvaP[er]->DoNotPropagate) 
      LabelKine->AddText(Form("%s: K=%.2f MeV  #theta_{lab}=%.1f deg  #phi_{lab}=%.1f deg",
			      EvaP[er]->Name.c_str(), EvaP[er]->GetKE(), EvaP[er]->GetTheta()*180/pi,
			      EvaP[er]->GetPhi()*180/pi));
    if (EvaR[er] && !EvaR[er]->DoNotPropagate) 
      LabelKine->AddText(Form("%s: K=%.2f MeV  #theta_{lab}=%.1f deg  #phi_{lab}=%.1f deg",
			      EvaR[er]->Name.c_str(), EvaR[er]->GetKE(), EvaR[er]->GetTheta()*180/pi,
			      EvaR[er]->GetPhi()*180/pi));
  }

  if (DeDau1 && DeDau2) {
    LabelKine->AddText(Form("%s: K=%.2f MeV  #theta_{lab}=%.1f deg  #phi_{lab}=%.1f deg",
			    DeDau1->Name.c_str(), DeDau1->GetKE(), DeDau1->GetTheta()*180/pi,
			    DeDau1->GetPhi()*180/pi));
    LabelKine->AddText(Form("%s: K=%.2f MeV  #theta_{lab}=%.1f deg  #phi_{lab}=%.1f deg",
			    DeDau2->Name.c_str(), DeDau2->GetKE(), DeDau2->GetTheta()*180/pi,
			    DeDau2->GetPhi()*180/pi));
  }
  else if (Heavy)
    LabelKine->AddText(Form("%s: K=%.2f MeV  #theta_{lab}=%.1f deg  #phi_{lab}=%.1f deg",
			    Heavy->Name.c_str(), Heavy->GetKE(), Heavy->GetTheta()*180/pi,
			    Heavy->GetPhi()*180/pi));
  if (Light)
    LabelKine->AddText(Form("%s: K=%.2f MeV  #theta_{lab}=%.1f deg  #phi_{lab}=%.1f deg",
			    Light->Name.c_str(), Light->GetKE(), Light->GetTheta()*180/pi,
			    Light->GetPhi()*180/pi));
  LabelKine->AddText(Form("#theta_{c.m.}=%.1f deg", this->theta_CM[0]*180/pi));



#endif
  return ReactionAllowed;
}


///////////////////////////////////////////////////////////////////////////////////
// This randomizes the detector response, making it more realistic as it is simiar
// to accunting for the energy loss uncertainty.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetStripEnergyResolution(float Sigma)
{
  EneSigma = Sigma;
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Function that passes the 'gSystem' pointer (generated in each root session to
// the corresponding pointer member of this class. The 'this->gSystem' pointer is 
// mainly used to keep track of the memory used by the program.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetROOTSystemPointer(TSystem* gSystem)
{
  this->gSystem = gSystem;
  cout << "gSystem = " << gSystem << endl;
  CheckMemoryUsage(1);
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Create the target particle object (currently the charge is not used).
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetTargetParticle(string Name)
{ 
#if 1
  // Use the nuclide finder object to determine the mass and atomic
  // number of this particle. Mass must be in MeV/c^2) and Z (in e).
  double m = NuF->GetMass(Name, "MeV/c^2");
  int Z = NuF->GetZ(Name);
  Target = new Particle(Name, m, Z);
#endif
  return;
}




///////////////////////////////////////////////////////////////////////////////////
// Core method of this class, where the simulation takes place. Logic:
// Loop over events
// 1. Set beam inital conditions (beam energy, position)
// 2. Within the selected strip randomly select the position at
//    which the beam particle interacts with the target and calculate
//    the kinetic energy at the reaction point
// 3. Set the kinematics of all particles at the reaction point
// 4. Propagate the beam particle (backwards in time) from the
//    reaction point to the entrance of MUSIC
// 5. Propagate evaporated (light) particle (p,n,4He)
// 6. Propagate evaporation residue (heavy particle)
// 7. Compute detector response (i.e. DE for beam + light + heavy + etc)
// 8. Display trace and particle trajecories
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::Simulate(int StpID, // set to -1 for unreacted beam
			       int NEvents,
			       double MaxTime,
			       double UserStep,
			       int UpdateVis, 
			       int Wait,
			       string FileName,
			       string FileOpt)
{
  double ti,xi,yi,zi, tf,xf,yf,zf;
#if 1
  // Verify that the anode geometry has been set
  if (VolAnode==0) {
    cout << "Anode geometry not specified. Use SetAnode method." << endl;
    return;
  }
  
  this->NEvents = NEvents;

  // For progress monitor
  TStopwatch StpWatch;
  long EvtsProcessed = 0;
  long double Frac[11] = {0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  int FIndex = 0;

  // ROOT file where the traces will be saved
  TFile* ROOTfile = 0;
  
  if (FileName!="")
    ROOTfile = new TFile(FileName.c_str(), FileOpt.c_str());

  // Tree similar to the one used for experimental data
  SimTree = InitTree(ROOTfile, FileOpt);

  TDirectory* trace_dir = 0;
  if (ROOTfile) {
    if (FileOpt=="update" || FileOpt=="UPDATE") {
      trace_dir = (TDirectory*)ROOTfile->Get("traces");
      trace_dir->cd();
    }
    else {
      trace_dir = ROOTfile->mkdir("traces");
      trace_dir->cd();
    }
  }

  // Create new traces and trajectories (objectrs) for visualizing the
  // detector response
  CreateTracesAndTrajectories(NEvents);
  
  cout << "Simulating " << NEvents << " MUSIC traces ... " << endl;
 
  // Get the average beam energy loss and print the exc energy of the
  // compound nucleus sampled by each strip.
  SetInitialKinematics(ctf.Kb); // ideal case with no energy straggling
  Particle* BeamInit = new Particle("beam init");  
  BeamInit->Copy(Beam);
  Particle* BeamCopy = new Particle("beam copy");
  BeamCopy->Copy(Beam);
  if (PrintLevel>0)
    BeamCopy->Print();
  PropagateParticle(BeamCopy, 0, MaxTime, UserStep, DeltaEB_ave);
  for (int stp=0; stp<AnodeRows; stp++)
    for (int col=0; col<AnodeCols+1; col++) 
      TraceUB[col]->SetPoint(stp, stp, DeltaEB_ave[stp][col]);
  
  // PrintEnergetics(ctf.Kb, DeltaEB_ave);

  //-------------------------------------------------------------------------------
  // Some kinematic variables
  double Kb_min, Kb_max, MinZ, MaxZ, MinT, MaxT;
  
  // Get the beam energy limits in the selected strip (assuming the
  // beam direction is parallel to the z-axis).
  MinZ = 0;
  MaxZ = AnodeDZ[0][0];
  for (int stp=1; stp<AnodeRows; stp++) {
    MinZ += AnodeDZ[stp-1][0];
    MaxZ += AnodeDZ[stp][0];
    if (AnodeStpID[stp][0]==StpID)
      break;
  }
  Kb_max = Beam->GetFinalEnergy(0, ctf.Kb, MinZ, 1E-3/*step size in cm*/);
  MinT = Beam->GetTimeOfFlight(0);
  Kb_min = Beam->GetFinalEnergy(0, ctf.Kb, MaxZ, 1E-3/*step size in cm*/);
  MaxT = Beam->GetTimeOfFlight(0);
  if (PrintLevel>0) {
    Log << "|---- Kinematic constraints for strip ";
    Log.width(3); Log << StpID;
    Log << " ---------\n"
	 << "|     |   In    |   Out   |  Units  |\n"
	 << "| zr  |";
    Log.width(9);    Log << MinZ;   Log << "|";
    Log.width(9);    Log << MaxZ;   Log << "|";
    Log.width(9);    Log << "cm";   Log << "|\n";
    Log << "| tof |";
    Log.width(9);    Log << MinT;   Log << "|";
    Log.width(9);    Log << MaxT;   Log << "|";
    Log.width(9);    Log << "ns";   Log << "|\n";
    Log << "| Kb  |";
    Log.width(9);    Log << Kb_max; Log << "|";
    Log.width(9);    Log << Kb_min; Log << "|";
    Log.width(9);    Log << "MeV";   Log << "|\n";
    Log << "|--------------------------------------------------" << endl; 
  }
  //-------------------------------------------------------------------------------



  //-------------------------------------------------------------------------------
  // Event for-loop
  //-------------------------------------------------------------------------------
  Log << "Initiating envent for-loop" << endl;
  for (int evt=0; evt<NEvents; evt++) {
    if (evt%1000==0)
      if (CheckMemoryUsage()==0) {
	Log << "Exiting event for-loop" << endl;
	break;
      }
    
    if (PrintLevel>0) {
      Log << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      Log << "!!!       EVENT " << evt << endl;
      Log << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" << endl;
    }
    ResetBranches();
    
    TraceCan->cd(1);
    HCT->Draw();
    TraceCan->cd(2);
    HPT->Draw();

    // Reset values of TTree branches
    
    // Reset the detector response
    for (int stp=0; stp<AnodeRows; stp++) 
      for (int col=0; col<AnodeCols+1; col++) {
	DeltaEB[stp][col] = 0;
	DeltaEL[stp][col] = 0;
	DeltaEH[stp][col] = 0;
	for (int er=0; er<CurEva; er++) {
	  DeltaE_EvaP[er][stp][col] = 0;
	  DeltaE_EvaR[er][stp][col] = 0;
	}
      }
    
    // 1. Set beam inital conditions (beam energy, position)
    double Ebeam = ctf.Kb;
    // If the user specified a KbFWHM>0 change the beam energy
    // assuming a gaussian distribution.
    if (ctf.KbFWHM>0.0) {
      Gaussian->SetParameters(1.0, 0.0, ctf.KbFWHM/2.355);
      Ebeam += Gaussian->GetRandom();
    }
    SetInitialKinematics(Ebeam);

    int ReacAllowed = 0;
    double Kbr = 0;
    double TOF = 0;
    if (StpID>-1) {
      // 2. Within the selected strip randomly select the position at
      // which the beam particle interacts with the target and calculate
      // the kinetic energy at the reaction point
      this->zr = Rdm->Uniform(MinZ, MaxZ);
      Kbr = Beam->GetFinalEnergy(0, Ebeam, this->zr, 1E-3/*cm*/);
      TOF = Beam->GetTimeOfFlight(0);
      if (PrintLevel>0)
	Log << "Kbr = " << Kbr << "  zr = " << this->zr << "  tof = " << TOF << endl;
      
      // 3. Set the kinematics of all particles at the reaction point
      ReacAllowed = SetReactionKinematics(Kbr, this->zr, TOF);
      // Check conservation of 4-momentum
      if (PrintLevel>0) {
	Log << "Conservation of 4-momentum at reaction point (zr)" 
	    << endl;
	FourVector Pi("initial 4-momemtum (lab)",0,0,0,0);
	Pi += Beam->GetP() + Target->GetP();
	FourVector Pf("final 4-momentum (lab)",0,0,0,0);
	for (int er=0; er<CurEva; er++) {
	  if (!EvaR[er]->DoNotPropagate)
	    Pf += EvaR[er]->GetP();
	  if (!EvaP[er]->DoNotPropagate)
	    Pf += EvaP[er]->GetP();
	}
	Pi.Print(Log);
	Pf.Print(Log);
      }
    }
    else {
      // Update the kinematics label and then draw it
      LabelKine->Clear();
      LabelKine->AddText("Kinematics");
      LabelKine->AddText(Form("beam: K=%.2f MeV", Ebeam));
    }
    
    if (ReacAllowed) {
      // 4. Propagate the beam particle (backwards in time) from the
      // reaction point to the entrance of MUSIC
      Beam->GetX(tf,xf,yf,zf);
      PropagateParticle(Beam, evt, MaxTime, -UserStep, DeltaEB);
      for (int stp=0; stp<AnodeRows; stp++)
	for (int col=0; col<AnodeCols+1; col++) 
	  TraceB[col]->SetPoint(stp, stp, DeltaEB[stp][col]);
      Beam->GetX(ti,xi,yi,zi);                      // <- This is not a mistake
      TrackBeam->SetOrigin(xi,yi,zi);
      TrackBeam->SetVector(xf-xi,yf-yi,zf-zi);
      
      // 5-6. Propagate outgoing particles (evaporation residues)
      for (int er=0; er<CurEva; er++) {
	
	// evaporated (light) particle (p,n,4He)
	EvaP[er]->GetX(ti,xi,yi,zi);
	PropagateParticle(EvaP[er], evt, MaxTime, UserStep, DeltaE_EvaP[er]);
	for (int stp=0; stp<AnodeRows; stp++)
	  for (int col=0; col<AnodeCols+1; col++) 
	    TraceEP[er][col]->SetPoint(stp, stp, DeltaE_EvaP[er][stp][col]);
	EvaP[er]->GetX(tf,xf,yf,zf);
	TrackEvaP[er]->SetOrigin(xi,yi,zi);
	TrackEvaP[er]->SetVector(xf-xi,yf-yi,zf-zi);
	xfl[er] = xf;
	yfl[er] = yf;
	zfl[er] = zf;
	
	// evaporation residue (heavy particle)
	EvaR[er]->GetX(ti,xi,yi,zi);
	PropagateParticle(EvaR[er], evt, MaxTime, UserStep, DeltaE_EvaR[er]);
	for (int stp=0; stp<AnodeRows; stp++)
	  for (int col=0; col<AnodeCols+1; col++) 
	    TraceER[er][col]->SetPoint(stp, stp, DeltaE_EvaR[er][stp][col]);
	EvaR[er]->GetX(tf,xf,yf,zf);
	if (!EvaR[er]->DoNotPropagate) {
	  xfe = xf;
	  yfe = yf;
	  zfe = zf;
	}
	TrackEvaR[er]->SetOrigin(xi,yi,zi);
	TrackEvaR[er]->SetVector(xf-xi,yf-yi,zf-zi);
      }
    }
    // reaction not allowed, only propagate beam
    else {
      //Beam->Copy(BeamInit);
      PropagateParticle(Beam, evt, MaxTime, UserStep, DeltaEB);
      for (int stp=0; stp<AnodeRows; stp++)
	for (int col=0; col<AnodeCols+1; col++) 
	  TraceB[col]->SetPoint(stp, stp, DeltaEB[stp][col]);
      if (StpID>-1) {
	cout << "Warninig: reaction energetically not allowed for event " << evt 
	     << " (Kbr= " << Kbr << " MeV)." << endl;
	Log << "Warninig: reaction energetically not allowed for event " << evt 
	    << " (Kbr= " << Kbr << " MeV)." << endl;
      }
    }
    
    // 7. Compute detector response (i.e. DE for beam + light + heavy)
    // Clone the particle trajectories
    ComputeDetectorResponse(evt, StpID, UpdateVis);
    //    ROOTfile->cd();
    SimTree->Fill();
    
    
    // 8. Display trace and particle trajecories
    if (UpdateVis)
      UpdateVisuals(evt, Kbr, this->zr, TOF, Wait);
    
    // Simple progress monitor
    if (NEvents>99) {
      if ((long double)(evt)>=Frac[FIndex]*NEvents) {
	cout << "\t" << Frac[FIndex]*100 << "% processed (" 
	     << StpWatch.RealTime() << " s)" << endl;
	StpWatch.Start(kFALSE);
	FIndex++;
      }
    }
    
    NTraces++;
  }
  StpWatch.Stop();
  StpWatch.Print();

  cout << "Envent for-loop concluded." << endl;
  CheckMemoryUsage(1);


  if (ROOTfile) {
    ROOTfile->cd();
    SimTree->Write("", TObject::kSingleKey);
    ROOTfile->Close();
  }

#endif
  return;
}



///////////////////////////////////////////////////////////////////////////////////
// With this version of Simulate() one can specify the angular range to be covered.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::Simulate(int StpID, double ThCMMin, double ThCMMax, int ThSteps,
			       double PhiCMMin, double PhiCMMax, int PhiSteps, double MaxTime,
			       double UserStep, int Wait)
{
  // Verify that the anode geometry has been set
  if (VolAnode==0) {
    cout << "Anode geometry not specified. Use SetAnode method." << endl;
    return;
  }
  
  // ROOT file where the traces will be saved
  TFile* ROOTfile = new TFile("sim.root", "recreate");
  
  // Tree similar to the one used for experimental data
  SimTree = InitTree(ROOTfile,"recreate");  

  // Angles in radians
  double theta = 0;
  double phi = 0;

  NEvents = PhiSteps*ThSteps*AnodeRows;
   
  // Create new traces and trajectories (objectrs) for visualizing the
  // detector response
  CreateTracesAndTrajectories(NEvents);
  
  cout << "Simulating MUSIC traces ... " << endl;
  
  SetInitialKinematics(ctf.Kb);   

  // Get the average beam energy loss and print the exc energy of the
  // compound nucleus sampled by each strip.
  Particle* BeamCopy = new Particle("beam copy");
  BeamCopy->Copy(Beam);
  if (PrintLevel>0)
    BeamCopy->Print();
  //  DeltaEB_ave = PropagateParticle(BeamCopy, ctf.Kb, MaxTime, UserStep); 
  PropagateParticle(BeamCopy, ctf.Kb, MaxTime, UserStep, DeltaEB_ave);
  for (int stp=0; stp<AnodeRows; stp++)
    for (int col=0; col<AnodeCols+1; col++) 
      TraceUB[col]->SetPoint(stp, stp, DeltaEB_ave[stp][col]);
  //  PrintCompoundEexc(ctf.Kb, DeltaEB_ave);
  
  //-------------------------------------------------------------------------------
  // Some kinematic variables
  double Kb_min, Kb_max, MinZ, MaxZ, MinT, MaxT;
  int evt = 0; 
  // Converting degrees in radians
  double theta_min = ThCMMin*pi/180;
  double theta_max = ThCMMax*pi/180;
  double phi_min = PhiCMMin*pi/180;
  double phi_max = PhiCMMax*pi/180;

  // Get the beam energy limits in the selected strip (assuming the
  // beam direction is parallel to the z-axis).
  MinZ = 0;
  MaxZ = AnodeDZ[0][0];
  for (int stp=1; stp<AnodeRows; stp++) {
    MinZ += AnodeDZ[stp-1][0];
    MaxZ += AnodeDZ[stp][0];
    if (AnodeStpID[stp][0]==StpID)
      break;
  }
  Kb_max = Beam->GetFinalEnergy(0, ctf.Kb, MinZ, 1E-3/*step size in cm*/);
  MinT = Beam->GetTimeOfFlight(0);
  Kb_min = Beam->GetFinalEnergy(0, ctf.Kb, MaxZ, 1E-3/*step size in cm*/);
  MaxT = Beam->GetTimeOfFlight(0);
  if (PrintLevel>0) {
    Log << "|---- Kinematic constraints for strip ";
    Log.width(3); Log << StpID;
    Log << " ---------\n"
	 << "|     |   In    |   Out   |  Units  |\n"
	 << "| zr  |";
    Log.width(9);    Log << MinZ;   Log << "|";
    Log.width(9);    Log << MaxZ;   Log << "|";
    Log.width(9);    Log << "cm";   Log << "|\n";
    Log << "| tof |";
    Log.width(9);    Log << MinT;   Log << "|";
    Log.width(9);    Log << MaxT;   Log << "|";
    Log.width(9);    Log << "ns";   Log << "|\n";
    Log << "| Kb  |";
    Log.width(9);    Log << Kb_max; Log << "|";
    Log.width(9);    Log << Kb_min; Log << "|";
    Log.width(9);    Log << "MeV";   Log << "|\n";
    Log << "|--------------------------------------------------" << endl; 
  }
  //-------------------------------------------------------------------------------

  // Event for-loop
  for (int ths=0; ths<ThSteps; ths++) {
    theta = ths*(theta_max - theta_min)/(ThSteps - 1) + theta_min;
    for (int phs=0; phs<PhiSteps; phs++) {
      phi = phs*(phi_max - phi_min)/(PhiSteps - 1) + phi_min;
      
      if (PrintLevel>0)
	Log << "\n***************** Event " << evt << "\n" << endl;
    
      TraceCan->cd(1);
      HCT->Draw();
      // TraceCan->cd(2);
      // HPT->Draw();
      
      // Reset the detector response
      for (int stp=0; stp<AnodeRows; stp++) 
	for (int col=0; col<AnodeCols+1; col++) {
	  DeltaEB[stp][col] = 0;
	  DeltaEL[stp][col] = 0;
	  DeltaEH[stp][col] = 0;
	}
      
      // 1. Set beam inital conditions (beam energy, position)
      SetInitialKinematics(ctf.Kb);   
      
      // 2. Within the selected strip randomly select the position at
      // which the beam particle interacts with the target and calculate
      // the kinetic energy at the reaction point
      this->zr = Rdm->Uniform(MinZ, MaxZ);
      double Kbr = Beam->GetFinalEnergy(0, ctf.Kb, this->zr, 1E-3);
      double TOF = Beam->GetTimeOfFlight(0);
      if (PrintLevel>0)
	Log << "Kbr = " << Kbr << "  zr = " << this->zr << "  tof = " << TOF << endl;
      
      // 3. Set the kinematics of all particles at the reaction point
      int ReacAllowed = SetReactionKinematics(Kbr, this->zr, TOF, theta, phi);
      if (ReacAllowed==0) {
	cout << "Warninig: reaction energetically not allowed for event " << evt 
	     << " (Kbr= " << Kbr << " MeV)." << endl;
	continue;
      }
      
      // 4. Propagate the beam particle (backwards in time) from the
      // reaction point to the entrance of MUSIC
      //    DeltaEB = PropagateParticle(Beam, evt, MaxTime, -UserStep);
      PropagateParticle(Beam, evt, MaxTime, -UserStep, DeltaEB);
      
      // 5. Propagate heavy particle and calculate energy loss in the
      // anode elements
      //      DeltaEH = PropagateParticle(Heavy, evt, MaxTime, UserStep);
      PropagateParticle(Heavy, evt, MaxTime, UserStep, DeltaEH);
      
      // 6. Propagate light particle and calculate energy loss in the
      // anode elements
      //      DeltaEL = PropagateParticle(Light, evt, MaxTime, UserStep);
      PropagateParticle(Light, evt, MaxTime, UserStep, DeltaEL);
      
      // 7. Compute detector response (i.e. DE for beam + light + heavy)
      // Clone the particle trajectories
      ComputeDetectorResponse(evt, StpID, 1);
      SimTree->Fill();
      
      // 8. Display trace and particle trajecories   
      UpdateVisuals(evt, Kbr, this->zr, TOF, Wait);
      
      NTraces++;
      evt++;
    }
  }

  if (ROOTfile) {
    ROOTfile->cd();
    SimTree->Write("", TObject::kSingleKey);
    ROOTfile->Close();
  }
  
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Private method, to be used in an event loop.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::UpdateVisuals(int evt, double Kbr, double zr, double TOF, int Wait)
{
#if 1
  if (PrintLevel>0)
    Log << "Update visuals: evt=" << evt << " Kbr=" << Kbr << " MeV zr=" << zr << " cm TOF=" 
	 << TOF << " ns Wait=" << Wait << "\n3D stuff ..." << endl;


  if (Wait) {
    // 3D stuff (only makes sense to do this when Wait=1)
    double tracklength = TrackBeam->GetVector().Mag();
    TrackBeam->SetTubeR(0.1/tracklength);
    TrackBeam->ElementChanged();

    short C,S,W;
    for (int er=0; er<CurEva; er++) {
      if (EvaP[er] && !EvaP[er]->DoNotPropagate) {
	EvaP[er]->GetTrajectoryAtt(C,S,W);
	tracklength = TrackEvaP[er]->GetVector().Mag();
	if (tracklength>0) {
	  TrackEvaP[er]->SetTubeR(0.1/tracklength);
	  TrackEvaP[er]->ElementChanged();
	}
      }
      if (EvaR[er] && !EvaR[er]->DoNotPropagate) {
	tracklength = TrackEvaR[er]->GetVector().Mag();
	if (tracklength>0) {
	  EvaR[er]->GetTrajectoryAtt(C,S,W);
	  TrackEvaR[er]->SetTubeR(0.1/tracklength);
	  TrackEvaR[er]->ElementChanged();
	}
      }
    }
    Eve->Redraw3D();
  }

  // 2D stuff
  if (PrintLevel>0)
    Log << "2D stuff..." << endl;

#if 1
  // Traces of energy loss in columns as a function of the strip
  // number (detector response).
  TraceCan->cd(1);
  LabelKine->Draw();          // Leged with reaction kinematics
  TraceUB[AnodeCols]->Draw("l same");
  for (int col=0; col<AnodeCols; col++)
    Trace[col]->Draw("l same");
  Trace[AnodeCols]->Draw("*l same");
  if (LegCol->GetNRows()==0) {
    LegCol->AddEntry(Trace[AnodeCols],"All columns","l");
    for (int col=0; col<AnodeCols; col++)
      LegCol->AddEntry(Trace[col], Form("Column %d", col),"l");
    LegCol->Draw();
  }

  // Traces of particles' energy loss as a function of the strip
  // number
  TraceCan->cd(2);  
  TraceUB[AnodeCols]->Draw("l same");
  TraceB[AnodeCols]->Draw("l same");
  for (int er=0; er<CurEva; er++) {
    if (TraceER[er][AnodeCols]->GetN()>0)
      TraceER[er][AnodeCols]->Draw("l same");
    if (TraceEP[er][AnodeCols]->GetN()>0)
      TraceEP[er][AnodeCols]->Draw("l same");
  }
  Trace[AnodeCols]->Draw("*l same");  // total detector response
  if (LegPart->GetNRows()==0) {
    LegPart->AddEntry(Trace[AnodeCols], "All particles", "l");
    LegPart->AddEntry(TraceB[AnodeCols], "beam", "l");
    for (int er=0; er<CurEva; er++) {
      if (TraceEP[er][AnodeCols]->GetN()>0)
	LegPart->AddEntry(TraceEP[er][AnodeCols], EvaP[er]->Name.c_str(), "l");
      if (TraceER[er][AnodeCols]->GetN()>0)
	LegPart->AddEntry(TraceER[er][AnodeCols], EvaR[er]->Name.c_str(), "l");
    }
  }
  LegPart->Draw();

#endif

  // Multiplicity 
  TraceCan->cd(3);
  TraceMult->GetYaxis()->SetRangeUser(0,AnodeCols+1);
  TraceMult->Draw("HIST");


  TraceCan->Update();
  if (Wait==1)
    TraceCan->WaitPrimitive();

#endif
  return;
}


///////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::WriteTraces(char* FileName)
{
#if 0
  TFile* Output = new TFile(FileName, "update");
  Output->cd();
  TDirectory* trace_dir = Output->mkdir("traces");
  trace_dir->cd();
  if (Trace!=0) {
    for (int n=0; n<NTraces; n++)
      if (Trace[n]!=0) {
	for (int col=0; col<AnodeCols+1; col++) {
	  if (col==AnodeCols) 
	    Trace[n][col]->Write(Form("Trace%d",n), TObject::kOverwrite);
	  else
	    Trace[n][col]->Write(Form("Trace%dc%d",n,col), TObject::kOverwrite);
	}
      }
  }
  Output->Close();
#endif
  return;
}



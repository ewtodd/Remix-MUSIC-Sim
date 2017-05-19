// Methods for MUSIC_Simulator class.
// See header file for class description and compilation instructions.

#include "MUSIC_Simulator.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////
// Constructor
///////////////////////////////////////////////////////////////////////////////////
MUSIC_Simulator::MUSIC_Simulator()
{
  // Welcome message
  cout << "================================================================================" << endl;
  cout << "|--- MUSIC simulator ----------------------------------------------------------|" << endl;
  cout << "| Written by Daniel Santiago-Gonzalez                                          |" << endl;
  cout << "| ver 2.0 (2016/4)                                                             |" << endl;
  cout << "| To get the latest version type:                                              |" << endl;
  cout << "| git clone https://dasago@bitbucket.org/music_anl/music_simulator.git         |" << endl;
  cout << "================================================================================" << endl;

  Name = "MUSIC_Simulator";

  // Initialize variables.
  CMEMax = 0;
  CMEMin = 0;
  EneSigma = 0;
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
  Gaussian = 0;
  Target = 0;
  Trace = 0;
  TraceH = 0;
  TraceD1 = 0;
  TraceD2 = 0;
  TraceL = 0;
  Compound = 0;
  Light = 0;
  Heavy = 0;
  DeDau1 = 0;
  DeDau2 = 0;

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
  
  // TEveManager for drawing 3D particle trajectories
  Eve = new TEveManager(960, 1018, kTRUE, "V");

  // Canvas and legends for traces
  TraceCan = new TCanvas("TraceCan","Traces", 0, 0, 960, 1018);
  TraceCan->Divide(1,2);
  TraceCan->cd(1)->SetGrid();
  TraceCan->cd(2)->SetGrid();
  //  TraceCan->cd(3)->SetGrid();
  LegCol = new TLegend(0.692,0.616,0.826,0.861);
  LegPart = new TLegend(0.692,0.616,0.826,0.861);
  LabelKine = new TPaveText(0.152,0.679,0.437,0.875,"NDC");

  
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
  double Kb = Kb_after_window;
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
  double Kb = Kb_after_window;
  double Kb_min;
  float TotalLength = 0;
  
 
  for (int i=0; i<AnodeStps; i++) 
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
  for (int i=0; i<AnodeStps; i++) {
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
    if (i+1<AnodeStps) 
      Eexc_beg = Eexc_end;
  }
#endif
  return;
}

///////////////////////////////////////////////////////////////////////////////////
// Simple function that returns 0 if the used RAM is more than the MaxMemory
// (or if the 'System' pointer hasn't been set) and 1 if the used RAM is 
// still less than MaxMemory.
///////////////////////////////////////////////////////////////////////////////////
// int MUSIC_Simulator::CheckMemoryUsage()
// {
//   int status = 0;
//   MemInfo_t* Memory;
//   if (System!=0) {
//     Memory = new MemInfo_t();
//     System->GetMemInfo(Memory);
//     cout << "> Total memory used: "  << Memory->fMemUsed << " MB = "
// 	 << 100.0*Memory->fMemUsed/Memory->fMemTotal << " % of total memory." << endl;
//     if (Memory->fMemUsed < MaxMemory)
//       status = 1;
//     delete Memory;
//   }
//   else
//     cout << "> Warning: System pointer not set." << endl;
//   return status;
// }



///////////////////////////////////////////////////////////////////////////////////
// Private method, to be used in an event loop.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::ComputeDetectorResponse(int evt)
{
#if 1
  double Ethresh = 0.02;  // Right now the threshold for the multiplicity is hard-coded here.
  double DeltaE = 0;
  TraceMult->Reset();
  if (evt<Particle::MaxEvents) {
    TrajH[evt] = (TEveStraightLineSet*)Heavy->AllTraj[evt]->Clone();
    if (DeDau1 && DeDau2) {
      TrajD1[evt] = (TEveStraightLineSet*)DeDau1->AllTraj[evt]->Clone();
      TrajD2[evt] = (TEveStraightLineSet*)DeDau2->AllTraj[evt]->Clone();
    }
    TrajL[evt] = (TEveStraightLineSet*)Light->AllTraj[evt]->Clone();
  }
  
  // Reset the SimTree leaves
  if (SimTree!=0) {
    cathode = 0;
    strip0 = 0;
    strip17 = 0;
    for (int stp=0; stp<ExpAnodeStps; stp++) {
      de_l[stp] = 0;
      de_r[stp] = 0;
    }
  }

  // Loop over the anode's strips and columns
  for (int stp=0; stp<AnodeStps; stp++) {
    DeltaE = 0;
    int mult = 0;
    for (int col=0; col<AnodeCols+1; col++) {
      // Previously we were normalizing to the energy loss of the
      // beam minus the energy loss in the the dead layer, which in
      // this version of the code is typically strip 1.  
      // DeltaE = DeltaEB[n] + DeltaEL[n] + DeltaEH[n] - (DeltaEB_ave[n]-DeltaEB_ave[1]);
      // In this version, we do not normalize to the average energy
      // loss of the beam in each strip.
      if (DeDau1 && DeDau2) {
	TraceD1[evt][col]->SetPoint(stp, stp, DeltaED1[stp][col]);
	TraceD2[evt][col]->SetPoint(stp, stp, DeltaED2[stp][col]);
      }
      else
	TraceH[evt][col]->SetPoint(stp, stp, DeltaEH[stp][col]);
      TraceL[evt][col]->SetPoint(stp, stp, DeltaEL[stp][col]);
      
      DeltaE = DeltaEB[stp][col] + DeltaEL[stp][col] + DeltaEH[stp][col] + 
	DeltaED1[stp][col] + DeltaED2[stp][col];
      
      if (DeltaE>Ethresh && col<AnodeCols)
	mult++;

      // Fill tree leaves
      if (SimTree!=0) {
	int stpid = AnodeStpID[stp][col];
	if (stpid>=0) {
	  cathode += DeltaE;
	  if (stpid==0)
	    strip0 += DeltaE;
	  else if (stpid==17)
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
       if (EneSigma!=0 && Gaussian!=0 && DeltaE>0) {
         Gaussian->SetRange(0.0, 2*DeltaE);
         Gaussian->SetParameters(1.0, DeltaE, EneSigma);
         DeltaE = Gaussian->GetRandom();
       }
      Trace[evt][col]->SetPoint(stp, stp, DeltaE);
      if (PrintLevel>0)
	cout << stp << " " << col << ": " << DeltaEB[stp][col] << " " << DeltaEL[stp][col] 
	     << " " << DeltaEH[stp][col] << " " << DeltaED1[stp][col] << " " << DeltaED2[stp][col] 
	     << endl;
    }
    TraceMult->Fill(stp, mult);
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
  for (int stp=0; stp<AnodeStps; stp++) {
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

  // Initialize beam traces for each column (this trace does not
  // depend on the number of events)
  TraceB = new TGraph*[AnodeCols+1];
  for (int col=0; col<AnodeCols+1; col++) {
    TraceB[col] = new TGraph();
    if (col==AnodeCols) {
      TraceB[col]->SetName("Beam trace");
      TraceB[col]->SetLineColor(kBlack);
      TraceB[col]->SetLineWidth(3);
    }
    else {
      TraceB[col]->SetName(Form("Beam trace col %d", col));
      TraceB[col]->SetLineColor(Chroma[col]);
      TraceB[col]->SetLineStyle(2);
      TraceB[col]->SetLineWidth(2);
    }
  }
  // Initialize the rest of the traces for each column
  Trace = new TGraph**[NEvents];
  TraceH = new TGraph**[NEvents];
  TraceD1 = new TGraph**[NEvents];
  TraceD2 = new TGraph**[NEvents];
  TraceL = new TGraph**[NEvents];
  for (int evt=0; evt<NEvents; evt++) {
    Trace[evt] = new TGraph*[AnodeCols+1];
    TraceH[evt] = new TGraph*[AnodeCols+1];
    TraceD1[evt] = new TGraph*[AnodeCols+1];
    TraceD2[evt] = new TGraph*[AnodeCols+1];
    TraceL[evt] = new TGraph*[AnodeCols+1];
    for (int col=0; col<AnodeCols+1; col++) {
      Trace[evt][col] = new TGraph();
      TraceH[evt][col] = new TGraph();
      TraceD1[evt][col] = new TGraph();
      TraceD2[evt][col] = new TGraph();
      TraceL[evt][col] = new TGraph();
      if (col==AnodeCols) {
	Trace[evt][col]->SetName(Form("evt %d total", evt));
	Trace[evt][col]->SetLineColor(kBlack);
	Trace[evt][col]->SetLineWidth(2);
	// heavy
	TraceH[evt][col]->SetName(Form("evt %d total %s", evt, Heavy->Name.c_str()));
	TraceH[evt][col]->SetLineColor(Heavy->GetColor());
	TraceH[evt][col]->SetLineWidth(2);
	if (DeDau1 && DeDau2) {
	  // decay daughter1
	  TraceD1[evt][col]->SetName(Form("evt %d total %s", evt, DeDau1->Name.c_str()));
	  TraceD1[evt][col]->SetLineColor(DeDau1->GetColor());
	  TraceD1[evt][col]->SetLineWidth(2);
	  // decay daughter2
	  TraceD2[evt][col]->SetName(Form("evt %d total %s", evt, DeDau2->Name.c_str()));
	  TraceD2[evt][col]->SetLineColor(DeDau2->GetColor());
	  TraceD2[evt][col]->SetLineWidth(2);
	}
	// light
	TraceL[evt][col]->SetName(Form("evt %d total %s", evt, Light->Name.c_str()));
	TraceL[evt][col]->SetLineColor(Light->GetColor());
	TraceL[evt][col]->SetLineWidth(2);
      }
      else {
	Trace[evt][col]->SetName(Form("evt %d col %d", evt, col));
	Trace[evt][col]->SetLineColor(Chroma[col]);
	Trace[evt][col]->SetLineStyle(2);
	Trace[evt][col]->SetLineWidth(2);
	// heavy
	TraceH[evt][col]->SetName(Form("evt %d col %d %s", evt, col, Heavy->Name.c_str()));
	TraceH[evt][col]->SetLineColor(Chroma[col]);
	TraceH[evt][col]->SetLineStyle(2);
	TraceH[evt][col]->SetLineWidth(2);
	if (DeDau1 && DeDau2) {
	  // decay daughter1
	  TraceD1[evt][col]->SetName(Form("evt %d col %d %s", evt, col, DeDau1->Name.c_str()));
	  TraceD1[evt][col]->SetLineColor(Chroma[col]);
	  TraceD1[evt][col]->SetLineStyle(2);
	  TraceD1[evt][col]->SetLineWidth(2);
	  // decay daughter1
	  TraceD2[evt][col]->SetName(Form("evt %d col %d %s", evt, col, DeDau2->Name.c_str()));
	  TraceD2[evt][col]->SetLineColor(Chroma[col]);
	  TraceD2[evt][col]->SetLineStyle(2);
	  TraceD2[evt][col]->SetLineWidth(2);
	}
	// light
	TraceL[evt][col]->SetName(Form("evt %d col %d %s", evt, col, Light->Name.c_str()));
	TraceL[evt][col]->SetLineColor(Chroma[col]);
	TraceL[evt][col]->SetLineStyle(2);
	TraceL[evt][col]->SetLineWidth(2);
      }
    }
  }

  TraceMult = new TH1I("TraceMult","Mult.",AnodeStps,-0.5,AnodeStps-0.5);
  TraceMult->GetXaxis()->SetTitle("Strip number");
  TraceMult->GetXaxis()->CenterTitle();
 
  // 3D trajectories
  TrajH = new TEveStraightLineSet*[NEvents];
  TrajD1 = new TEveStraightLineSet*[NEvents];
  TrajD2 = new TEveStraightLineSet*[NEvents];
  TrajL = new TEveStraightLineSet*[NEvents];
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
    Xaxis->SetMainColor(kGreen); 
    gEve->AddElement(Xaxis);
    TEveArrow* Yaxis = new TEveArrow(0,20,0,0,-10,0);
    Yaxis->SetMainColor(kYellow);
    gEve->AddElement(Yaxis);
    TEveArrow* Zaxis = new TEveArrow(0,0,30,0,0,0);
    Zaxis->SetMainTransparency(65);
    gEve->AddElement(Zaxis);
    gEve->Redraw3D(kTRUE);
  }
#endif
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Angles in degrees
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::GenerateTraceDatabase(string FileName, 
					    double ThCMMin, double ThCMMax, int ThSteps,
					    double PhiCMMin, double PhiCMMax, int PhiSteps,
					    double MaxTime, double UserDT, int Update,
					    int Wait)
{
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
  SimTree = new TTree("simt","Simulated MUSIC data");
  SimTree->Branch("de_l",  de_l,     Form("de_l[%d]/F",AnodeStps));
  SimTree->Branch("de_r",  de_r,     Form("de_r[%d]/F",AnodeStps));
  SimTree->Branch("seg",   seg,      Form("seg[%d]/I",AnodeStps));
  SimTree->Branch("stp0",  &strip0,  "stp0/F");
  SimTree->Branch("stp17", &strip17, "stp17/F");
  SimTree->Branch("cath",  &cathode, "cath/F");
  // The following branches are for physical quantities that at the
  // moment can only be obtained from the simulation
  SimTree->Branch("Kb", &Kb, "Kb/F");
  SimTree->Branch("Kl", &Kl, "Kl/F");
  SimTree->Branch("Kh", &Kh, "Kh/F");
  SimTree->Branch("theta_CM", &theta_CM, "theta_CM/F");
  SimTree->Branch("theta_l", &theta_l, "theta_l/F");
  SimTree->Branch("theta_h", &theta_h, "theta_h/F");
  SimTree->Branch("phi_l",   &phi_l,   "phi_l/F");
  SimTree->Branch("phi_h",   &phi_l,   "phi_h/F");
  if (PrintLevel>0)
    SimTree->Print();

  // Angles in radians
  double theta = 0;
  double phi = 0;
  
  NEvents = PhiSteps*ThSteps*AnodeStps;
  NTraces = 0;
    
  // Create new traces and trajectories (objectrs) for visualizing the
  // detector response
  CreateTracesAndTrajectories(NEvents);
  
  cout << "Generating " << NEvents << " MUSIC traces ..." << endl;
  
  SetInitialKinematics(Kb_after_window);   

  // Get the average beam energy loss and print the exc energy of the
  // compound nucleus sampled by each strip.
  Particle* BeamCopy = new Particle("beam copy");
  BeamCopy->Copy(Beam);
  if (PrintLevel>0)
    BeamCopy->Print();
  DeltaEB_ave = PropagateParticle(BeamCopy, Kb_after_window, MaxTime, UserDT); 
  for (int stp=0; stp<AnodeStps; stp++)
    for (int col=0; col<AnodeCols+1; col++) 
      TraceB[col]->SetPoint(stp, stp, DeltaEB_ave[stp][col]);
  // Draw the beam trace in the 3rd pad (over a background histogram)
  // TraceCan->cd(3);
  // HCTB->Draw();
  // for (int col=0; col<AnodeCols; col++)
  //   TraceB[col]->Draw("l same");
  // TraceB[AnodeCols]->Draw("*l same");
  PrintEnergetics(Kb_after_window, DeltaEB_ave);

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
  for (int stp_base=0; stp_base<AnodeStps; stp_base++) {
    if (AnodeStpID[stp_base][0]>=0) {
      
      // Get the beam energy limits in the selected strip (assuming
      // the beam direction is parallel to the z-axis).
      MinZ = 0;
      MaxZ = AnodeDZ[0][0];
      for (int stp=1; stp<AnodeStps; stp++) {
	MinZ += AnodeDZ[stp-1][0];
	MaxZ += AnodeDZ[stp][0];
	if (AnodeStpID[stp][0]==AnodeStpID[stp_base][0])
	  break;
      }
      Kb_max = Beam->GetFinalEnergy(0, Kb_after_window, MinZ, 1E-3/*step size in cm*/);
      MinT = Beam->GetTimeOfFlight(0);
      Kb_min = Beam->GetFinalEnergy(0, Kb_after_window, MaxZ, 1E-3/*step size in cm*/);
      MaxT = Beam->GetTimeOfFlight(0);      
      if (PrintLevel>0) {
	cout << "|---- Kinematic constraints for strip ";
	cout.width(3); cout << AnodeStpID[stp_base][0];
	cout << " ---------\n"
	     << "|     |   In    |   Out   |  Units  |\n"
	     << "| zr  |";
	cout.width(9);    cout << MinZ;   cout << "|";
	cout.width(9);    cout << MaxZ;   cout << "|";
	cout.width(9);    cout << "cm";   cout << "|\n";
	cout << "| tof |";
	cout.width(9);    cout << MinT;   cout << "|";
	cout.width(9);    cout << MaxT;   cout << "|";
	cout.width(9);    cout << "ns";   cout << "|\n";
	cout << "| Kb  |";
	cout.width(9);    cout << Kb_max; cout << "|";
	cout.width(9);    cout << Kb_min; cout << "|";
	cout.width(9);    cout << "MeV";   cout << "|\n";
	cout << "|--------------------------------------------------" << endl; 
      }
      for (int ths=0; ths<ThSteps; ths++) {
	theta = ths*(theta_max - theta_min)/(ThSteps - 1) + theta_min;
	for (int phs=0; phs<PhiSteps; phs++) {
	  phi = phs*(phi_max - phi_min)/(PhiSteps - 1) + phi_min;
	  
	  if (PrintLevel>0)
	    cout << "\n***************** Event " << evt << "\n" << endl;
	  
	  TraceCan->cd(1);
	  HCT->Draw();
	  TraceCan->cd(2);
	  HPT->Draw();
	  
	  // Reset the detector response
	  for (int stp=0; stp<AnodeStps; stp++) 
	    for (int col=0; col<AnodeCols+1; col++) {
	      DeltaEB[stp][col] = 0;
	      DeltaEL[stp][col] = 0;
	      DeltaEH[stp][col] = 0;
	      DeltaED1[stp][col] = 0;
	      DeltaED2[stp][col] = 0;
	    }
	  
	  // 1. Set beam inital conditions (beam energy, position)
	  SetInitialKinematics(Kb_after_window);   
	  
	  // 2. Within the selected strip randomly select the position
	  // at which the beam particle interacts with the target and
	  // calculate the kinetic energy at the reaction point
	  double zr = Rdm->Uniform(MinZ, MaxZ);
	  double Kbr = Beam->GetFinalEnergy(0, Kb_after_window, zr, 1E-3/*cm*/);
	  double TOF = Beam->GetTimeOfFlight(0);
	  if (PrintLevel>0)
	    cout << "Kbr = " << Kbr << "  zr = " << zr << "  tof = " << TOF << endl;
	  
	  // 3. Set the kinematics of all particles at the reaction point
	  int ReacAllowed = SetReactionKinematics(Kbr, zr, TOF, theta, phi);
	  if (ReacAllowed==0) {
	    cout << "Warninig: reaction energetically not allowed for event " << evt 
		 << " (Kbr= " << Kbr << " MeV)." << endl;
	    continue;
	  }
	  
	  // 4. Propagate the beam particle (backwards in time) from
	  // the reaction point to the entrance of MUSIC
	  DeltaEB = PropagateParticle(Beam, evt, MaxTime, -UserDT);
	  
	  // 5. Propagate heavy particle (or decay daughters) and calculate energy 
	  // loss in the anode elements
	  if (DeDau1 && DeDau2) {
	    DeltaED1 = PropagateParticle(DeDau1, evt, MaxTime, UserDT);
	    DeltaED2 = PropagateParticle(DeDau2, evt, MaxTime, UserDT);
	  }
	  else
	    DeltaEH = PropagateParticle(Heavy, evt, MaxTime, UserDT);
	  
	  // 6. Propagate light particle and calculate energy loss in the
	  // anode elements
	  DeltaEL = PropagateParticle(Light, evt, MaxTime, UserDT);
	  
	  // 7. Compute detector response (i.e. DE for beam + light + heavy)
	  // Clone the particle trajectories
	  ComputeDetectorResponse(evt);
	  SimTree->Fill();
	  
	  // 8. Display trace and particle trajecories   
	  if (Update) 
	    UpdateVisuals(evt, Kbr, zr, TOF, Wait);
	  
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
  cout << "Ebmax[MeV]";
  cout.width(15);
  cout << "EbMin[MeV]" << endl;
  cout.width(15);
  cout << "Exmax[MeV]";
  cout.width(15);
  cout << "Exmin[MeV]" << endl;
  cout.width(15);
  cout << "ECMmax[MeV]";
  cout.width(15);
  cout << "ECMmin[MeV]" << endl;

  for (int stp=0; stp<AnodeStps+1; stp++) {
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
    cout << sqrt(Ptot*Ptot) - mc;
    // Initial momentum in the CM reference frame (Ptot*Ptot is an invariant quantity).
    double pCM_max = sqrt((Ptot*Ptot-pow(mt+mb,2))*(Ptot*Ptot-pow(mb-mt,2))/(4*(Ptot*Ptot)));    
    // Lower limit
    Kb -= DeltaEB[stp][AnodeCols];
    cout.width(15);
    if (Kb>0) {
      // Linear momentum and total energy of the beam particle in the lab.
      pb = sqrt(2*mb*Kb*(1 + Kb/2/mb));
      Eb = sqrt(mb*mb + pb*pb);
      // Total four momentum in the lab.
      Ptot.SetCoords(Eb+mt, 0, 0, pb);
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
    else
      cout << "0" << endl;


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
double** MUSIC_Simulator::PropagateParticle(Particle* PO, int Event, double MaxTime, double UserDT)
{
  double** DE = new double*[AnodeStps];
  for (int stp = 0; stp<AnodeStps; stp++) {
    DE[stp] = new double[AnodeCols+1];
    for (int col = 0; col<AnodeCols+1; col++) 
      DE[stp][col] = 0;
  }
#if 1
  TGeoVolume* Vol = 0;
  bool Skip = 0;
  int step = 0;
  // Get the initial conditions from the Particle object.
  double m = PO->Mass;
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
  PO->AllTraj[Event]->SetName(Form("%s evt %d", PO->Name.c_str(), Event));
  
  if (PrintLevel>0)
    cout << "Propagator! " << PO->Name << " " << MaxTime << "  ti=" << ti << endl;
  
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
    double vx = c*p_mag*cos(phi)*sin(theta)/Ene;  // maybe change to px/Ene
    double vy = c*p_mag*sin(phi)*sin(theta)/Ene;  // maybe change to py/Ene
    double vz = c*p_mag*cos(theta)/Ene;  // maybe change to pz/Ene
    double vel = sqrt(vx*vx+vy*vy+vz*vz);
    
    // Very short time step for dense media, equivalent to distance
    // steps of 0.1 um. Let's just hope that the numerator is never
    // zero.
    double Dt = 1E-3/vel;   // from p/c = m*dx/dt -> dt = c*m*dx/p
    if (UserDT<0)
      Dt = -Dt;
    double tf = ti + Dt;
    double xf = xi + vx*Dt;
    double yf = yi + vy*Dt;
    double zf = zi + vz*Dt;

    if (PrintLevel>0)
      if (PrintLevel>1 || step==0) {
	cout << step << " " << PO->Name << " (";
	cout.precision(7);
	cout << tf << ", " << xf << " " << yf << " " << zf << ")  Dt=" << Dt;
      } 
    
    // Exit the while loop if the particl leaves the anode volume.
    if (zf>AnodeDepth || zf<0 || xf>AnodeLength/2 || xf<-AnodeLength/2 || 
	yf>AnodeHeight/2 || yf<-AnodeHeight/2) {
      if (PrintLevel>0)
	cout << PO->Name << ": out!" << endl;
      break;
    }
    
    // Find the anode segment in which the particle is moving.
    int stp = -1;
    int col = -1;
    Vol = Geo->FindNode(xf, yf, zf)->GetVolume(); 
    for (int s=0; s<AnodeStps; s++) 
      for (int c=0; c<AnodeCols; c++) 
	if (Vol==VolAnode[s][c]) {
	  stp = s;
	  col = c;
	  break;
	} 
    // Exit the while loop if the anode segment was not found.
    if (stp==-1 || col==-1) {
      if (PrintLevel>0)
	cout << PO->Name<< ": anode segment not found." << endl;
      break;
    }

    // Get the energy loss for the initial energy (Ki) in the gas
    // medium (0) over a small (differential) path length (dist).
    double dist = sqrt(pow(xf-xi,2) + pow(yf-yi,2) + pow(zf-zi,2));
    double Kf = 0;
    if (UserDT>0) 
      Kf = PO->GetFinalEnergy(0, Ki, dist, dist/10);
    else
      Kf = PO->GetInitialEnergy(0, Ki, dist, dist/10);

    // Exit the while loop if the particle has stopped.
    if (Kf<=0) {
      if (PrintLevel>0)
	cout << "Less than ZERO! " << Kf << endl;
      break;
    }

    DE[stp][col] += fabs(Ki - Kf);
    DE[stp][AnodeCols] += fabs(Ki - Kf);
    
    if (PrintLevel>0) 
      if (PrintLevel>1 || step==0) 
	cout << "\t d=" << dist << " \tKf=" << Kf << " \tDE=" << DE[stp][col] << endl; 
    
    // Get the momentum magnitude using the new (lower) kinetic energy
    // (I'm using a relativistic formula, although it is unlikely for
    // us to use relativistic energies).
    p_mag = sqrt(2*m*Kf*(1+Kf/2/m));
    // Reduce (or increase) the total energy by amount of energy
    // deposited in the medium.
    if (UserDT>0)
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
    if (tf-tt>UserDT) {   
      if (PO->SaveTrajectory) {
	if (!Skip) {
	  PO->AllTraj[Event]->AddLine(xt,yt,zt, xf,yf,zf);
	  Skip = 1;
	}
	else 
	  Skip = 0;
      }
      tt = tf;
      xt = xf;
      yt = yf;
      zt = zf;
      Kt = Kf;
      PO->SetTracePoint((float)tt, (float)xt, (float)yt, (float)zt, (float)Kt);
    }
  } // end 
#endif
  return DE;
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
  AnodeStps = 0;
  AnodeCols = 0;

  GeomFile.open(AnodeGeomFile.c_str());
  if (!GeomFile.is_open()) 
    cout << "ERROR: Anode geometry file \"" << AnodeGeomFile << "\" couldn't be opened." << endl;
  else {
    cout << "Loading anode geometry loaded from \"" << AnodeGeomFile << "\"." << endl;

    // Get the number of strips and columns for the anode segments
    do {
      getline(GeomFile, line);
      if (line.find("strips:")!=string::npos) 
	AnodeStps = atoi(line.substr(line.find(':')+1).c_str());
      if (line.find("columns:")!=string::npos)
	AnodeCols = atoi(line.substr(line.find(':')+1).c_str());     
    } while (!GeomFile.eof());
    GeomFile.close();
    
    cout << "Anode strips: " << AnodeStps << endl;
    cout << "Anode columns: " << AnodeCols << endl;    
    if (AnodeStps>0 && AnodeCols>0) {
      // Initialize the anode segment dimensions
      AnodeDX = new double*[AnodeStps];
      AnodeDY = new double*[AnodeStps];
      AnodeDZ = new double*[AnodeStps];
      AnodeColor = new short*[AnodeStps];
      AnodeSegName = new string*[AnodeStps];
      AnodeStpID = new int*[AnodeStps];
      for (int stp=0; stp<AnodeStps; stp++) {
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
	cout << "Total lines: " << line_counter << endl;
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
	if (PrintLevel>1)
	  cout << nl << "\t" << AnodeStpID[stp][col] << "\t" << AnodeSegName[stp][col] << "\t" 
	       << AnodeDX[stp][col] << "\t" << AnodeDY[stp][col] << "\t" << AnodeDZ[stp][col] 
	       << "\t" << AnodeColor[stp][col] << endl;
      }
      GeomFile.close();
    
      // The total anode depth (distance along the z axis) is the sum of
      // all segment depths for the first column
      for (int stp=0; stp<AnodeStps; stp++)
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
      VolAnode = new TGeoVolume**[AnodeStps];
      for (int stp=0; stp<AnodeStps; stp++) {
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
      DeltaEB_ave = new double*[AnodeStps];// average beam energy loss
      DeltaEB = new double*[AnodeStps];    // beam
      DeltaEL = new double*[AnodeStps];    // light
      DeltaEH = new double*[AnodeStps];    // heavy
      DeltaED1 = new double*[AnodeStps];   // decay daughter1
      DeltaED2 = new double*[AnodeStps];   // decay daughter2
      for (int stp=0; stp<AnodeStps; stp++) {
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
    } // end if (AnodeStps>0 && AnodeCols>0)
    
    HCTB = new TH2F("HCTB","Beam", AnodeStps,-0.5, AnodeStps-0.5, ELossBins,0,MaxELoss);
    HCTB->GetXaxis()->SetTitle("Strip number");
    HCTB->GetXaxis()->CenterTitle();
    HCTB->GetYaxis()->SetTitle("Energy loss [MeV]");
    HCTB->GetYaxis()->CenterTitle(); 
    
    HCT = new TH2F("HCT","Column traces", AnodeStps,-0.5, AnodeStps-0.5, ELossBins,0,MaxELoss);
    HCT->GetXaxis()->SetTitle("Strip number");
    HCT->GetXaxis()->CenterTitle();
    HCT->GetYaxis()->SetTitle("Energy loss [MeV]");
    HCT->GetYaxis()->CenterTitle(); 
    
    HPT = new TH2F("HPT","Particle traces", AnodeStps,-0.5, AnodeStps-0.5, ELossBins,0,MaxELoss);
    HPT->GetXaxis()->SetTitle("Strip number");
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
void MUSIC_Simulator::SetBeamParticle(string Name, int Color, string ELossFile, double K)
{
#if 1
  // Use the nuclide finder object to determine the mass and atomic
  // number of this particle. Mass must be in MeV/c^2) and Z (in e).
  double m = NuF->GetMass(Name, "MeV/c^2");
  int Z = NuF->GetZ(Name);
  Beam = new Particle(Name, m, Z, true /*SaveTrajectories on*/);
  Beam->SetTrajectoryAtt((short)Color);
  // Currently, this simulation is restricted to one medium (gas) in MUSIC.
  Beam->SetMedium(ELossFile);
  Kb_after_window = K;
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
  Compound = new Particle(Name, m, Z, true /*SaveTrajectories on*/);
  //  Compound->SetTrajectoryAtt((short)Color);
  // Need excitation energy based on kinematics

  //  Compound->SetExcEnergies(NEexc, Eexc);
  // I would be good to print the excitation energy of the compound 
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
  DeDau1 = new Particle(Name, m, Z, true /*SaveTrajectories on*/);
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
#if 1
  // Use the nuclide finder object to determine the mass and atomic
  // number of this particle. Mass must be in MeV/c^2) and Z (in e).
  double m = NuF->GetMass(Name, "MeV/c^2");
  int Z = NuF->GetZ(Name);
  DeDau2 = new Particle(Name, m, Z, true /*SaveTrajectories on*/);
  DeDau2->SetTrajectoryAtt((short)Color);
  // Currently, this simulation is restricted to one medium (gas) in MUSIC.
  DeDau2->SetMedium(ELossFile);
#endif
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
  Heavy = new Particle(Name, m, Z, true /*SaveTrajectories on*/);
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
  if (PrintLevel>0)
    Beam->Print();
#endif
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Create the light particle (evaporation residue) object. It is assumed that the 
// light particle (e.g. p, n, alpha) is in its ground state.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetLightParticle(string Name, int Color, string ELossFile)
{
#if 1
  // Use the nuclide finder object to determine the mass and atomic
  // number of this particle. Mass must be in MeV/c^2) and Z (in e).
  double m = NuF->GetMass(Name, "MeV/c^2");
  int Z = NuF->GetZ(Name);
  Light = new Particle(Name, m, Z, true /*SaveTrajectories on*/);
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
    cout << "Information printing level = " << Level;
    PrintLevel = Level;
    Log.open("music_sim.log");
    Log << "================================================================================" << endl;
    Log << "|--- MUSIC simulator log file -------------------------------------------------|" << endl;
  }
#endif
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Private method, to be used in an event loop. Establish the kinematics of the 
// particles at the reaction point. If theta_CM and phi_CM are both -1 (default 
// values) then their values are assigned randomly.
///////////////////////////////////////////////////////////////////////////////////
int MUSIC_Simulator::SetReactionKinematics(double Kbr/*MeV*/, double zr/*cm*/, double tof/*ns*/,
					   double theta_CM, double phi_CM)
{
  int ReactionAllowed = 1;
#if 1
  if (PrintLevel>0)
    cout << "\n*** Reaction kinematics ********************************************" << endl;
  double mb = Beam->Mass;
  double mt = Target->Mass;
  double mc = Compound->Mass;
  double ml = Light->Mass;
  double mh = Heavy->Mass;
  double Ex = Heavy->GetEexc();

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
    cout << "Center-of-mass velocity (v/c):" << endl;
    cout << "\tBetaX=" << BetaX << "  BetaY=" << BetaY << "  BetaZ=" << BetaZ << endl; 
    cout << "--- Beam and target particles --------------------------------------" << endl;
  }

  // Four-momentum of the beam.
  Beam->SetP(Eb, pb*sin(theta_b)*cos(phi_b), pb*sin(theta_b)*sin(phi_b), pb*cos(theta_b));
  Beam->SetX(tof, 0, 0, zr);
  if (PrintLevel>0)
    Beam->Print();
  
  // Four-momentum of the target.
  Target->SetP(mt, 0, 0, 0);
  Target->SetX(tof, 0, 0, zr);
  if (PrintLevel>0)
    Target->Print();
  
  // Total four momentum in the lab.
  FourVector Ptot = Beam->GetP() + Target->GetP();
  //   Ptot.SetName("Total four-mom. in the lab");
  
  // Exc. energy of compound particle.
  Compound->SetP(Ptot);
  if (PrintLevel>0) {
    cout << "--- Compound particle ----------------------------------------------" << endl;
    Compound->Print();
  }
  if (PrintLevel>0)
    cout << "Eexc(" << Compound->Name << ") = " << sqrt(Ptot*Ptot) - mc << " MeV" << endl;

  // If the user did not specify the value of theta and phi, randomly
  // select the scattering angle in the center of mass.
  if (theta_CM==-1 && phi_CM==-1) {
    theta_CM = acos(Rdm->Uniform(-1,1));
    phi_CM = Rdm->Uniform(-pi,pi);
  }
  if (PrintLevel>0) {
    cout << "--- Outgoing particles (heavy and light) ---------------------------" << endl;
    cout << "theta_cm=" << theta_CM*180/pi << "   phi_cm=" << phi_CM*180/pi << endl;
  }
  
  // Verify that the reaction can occur.
  // In this case, the reaction is NOT energetically allowed.
  if (Ptot*Ptot<pow(ml+mh+Ex,2)) {
    Light->SetP(ml, 0, 0, 0);
    Light->SetX(tof, 0, 0, zr);
    Heavy->SetP(mh+Ex, 0, 0, 0);
    Heavy->SetX(tof, 0, 0, zr);
    ReactionAllowed = 0;
  }
  // In this case, the reaction IS energetically allowed.
  else {
    // Final momentum in the CM reference frame (Ptot*Ptot is an invariant quantity).
    double pf_CM = sqrt((Ptot*Ptot - pow(ml+mh+Ex,2))*(Ptot*Ptot - pow(ml-mh-Ex,2))/(4*(Ptot*Ptot)));
    
    // Set the four-momentum components of the light particle in the center of mass.
    double plxCM = -pf_CM*sin(theta_CM)*cos(phi_CM);
    double plyCM = -pf_CM*sin(theta_CM)*sin(phi_CM);
    double plzCM = -pf_CM*cos(theta_CM);
    double ElCM = sqrt(ml*ml + pf_CM*pf_CM);
    Light->SetP(ElCM, plxCM, plyCM, plzCM);
    
    // Do a Lorentz transformation (boost) into the lab reference frame.
    // I've double checked the sign of the boost and is correct (-Beta).
    Light->Boost(-BetaX, -BetaY, -BetaZ);
    
    // Initial position of the light particle (at the target).
    Light->SetX(tof, 0, 0, zr);
    if (PrintLevel>0)
      Light->Print();
    
    // Four-momentum of the heavy recoil (lab).
    Heavy->SetP(Beam->GetP() + Target->GetP() - Light->GetP());
    Heavy->SetX(tof, 0, 0, zr);
    if (PrintLevel>0)
      Heavy->Print();

    // If the excitation energy allows the decay channel d1+d2 (provided that the user defined it)
    // then make the decay happen.
    double md1=0;      double md2=0;
    if (DeDau1 && DeDau2) {
      md1 = DeDau1->Mass;
      md2 = DeDau2->Mass;
      if ((md1+md2)<(mh+Ex)) {
	// Final momentum of the daughter particles in the reference frame of the heavy particle.
	FourVector Ph = Heavy->GetP();
	double pfd = sqrt((Ph*Ph - pow(md1+md2,2))*(Ph*Ph - pow(md1-md2,2))/(4*(Ph*Ph)));
	double theta_d1 = acos(Rdm->Uniform(-1,1));
	double phi_d1 = Rdm->Uniform(-pi,pi);
	DeDau1->SetP(sqrt(md1*md1 + pfd*pfd), pfd*sin(theta_d1)*cos(phi_d1),
		     pfd*sin(theta_d1)*sin(phi_d1), -pfd*cos(theta_d1));
	DeDau2->SetP(sqrt(md2*md2 + pfd*pfd), -pfd*sin(theta_d1)*cos(phi_d1),
		     -pfd*sin(theta_d1)*sin(phi_d1), pfd*cos(theta_d1));
	// Boost the two daughters from the ref. frame of the heavy particle to the lab.
	double BetaHX;	double BetaHY;	double BetaHZ;
	Heavy->GetBeta(BetaHX, BetaHY, BetaHZ);	
	DeDau1->Boost(-BetaHX, -BetaHY, -BetaHZ);
	DeDau1->SetX(tof, 0, 0, zr);
	DeDau2->Boost(-BetaHX, -BetaHY, -BetaHZ);
	DeDau2->SetX(tof, 0, 0, zr);
      }
      else 
	ReactionAllowed = 0;
    }

    // Fill the leaves related to the reaction kinematics
    Kb = Beam->GetKE();
    Kl = Light->GetKE();
    Kh = Heavy->GetKE();
    this->theta_CM = theta_CM*180/pi;
    this->phi_CM = phi_CM*180/pi;
    theta_l = (Light->GetTheta())*180/pi;
    phi_l = (Light->GetPhi())*180/pi;
    theta_h = (Heavy->GetTheta())*180/pi;
    phi_h = (Heavy->GetPhi())*180/pi;
  }
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
  Gaussian = new TF1("Gaussian","gaus",0, 100);
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
// 2. Randomly select a reaction point (based on the beam energy range in MUSIC)
// 3. Propagate beam to reaction point and calculate energy loss in the anode elements 
// 4. Set initial conditions for heavy and light particles
// 5. Propagate heavy particle and calculate energy loss in the anode elements 
// 6. Propagate light particle and calculate energy loss in the anode elements
// 7. Compute detector response (i.e. DE for beam + light + heavy)
// 8. Display trace and particle trajecories
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::Simulate(int StpID, int NEvents, double MaxTime, double UserDT, int Wait)
{
#if 1
  // Verify that the anode geometry has been set
  if (VolAnode==0) {
    cout << "Anode geometry not specified. Use SetAnode method." << endl;
    return;
  }
  
  this->NEvents = NEvents;

  // Create new traces and trajectories (objectrs) for visualizing the
  // detector response
  CreateTracesAndTrajectories(NEvents);
  
  cout << "Simulating MUSIC traces ... " << endl;
  
  SetInitialKinematics(Kb_after_window);   

  // Get the average beam energy loss and print the exc energy of the
  // compound nucleus sampled by each strip.
  Particle* BeamCopy = new Particle("beam copy");
  BeamCopy->Copy(Beam);
  if (PrintLevel>0)
    BeamCopy->Print();
  DeltaEB_ave = PropagateParticle(BeamCopy, Kb_after_window, MaxTime, UserDT); 
  for (int stp=0; stp<AnodeStps; stp++)
    for (int col=0; col<AnodeCols+1; col++) 
      TraceB[col]->SetPoint(stp, stp, DeltaEB_ave[stp][col]);
  // Draw the beam trace in the 3rd pad (over a background histogram)
  // TraceCan->cd(3);
  // HCTB->Draw();
  // for (int col=0; col<AnodeCols; col++)
  //   TraceB[col]->Draw("l same");
  // TraceB[AnodeCols]->Draw("*l same");
  PrintEnergetics(Kb_after_window, DeltaEB_ave);

  //-------------------------------------------------------------------------------
  // Some kinematic variables
  double Kb_min, Kb_max, MinZ, MaxZ, MinT, MaxT;
  
  // Get the beam energy limits in the selected strip (assuming the
  // beam direction is parallel to the z-axis).
  MinZ = 0;
  MaxZ = AnodeDZ[0][0];
  for (int stp=1; stp<AnodeStps; stp++) {
    MinZ += AnodeDZ[stp-1][0];
    MaxZ += AnodeDZ[stp][0];
    if (AnodeStpID[stp][0]==StpID)
      break;
  }
  Kb_max = Beam->GetFinalEnergy(0, Kb_after_window, MinZ, 1E-3/*step size in cm*/);
  MinT = Beam->GetTimeOfFlight(0);
  Kb_min = Beam->GetFinalEnergy(0, Kb_after_window, MaxZ, 1E-3/*step size in cm*/);
  MaxT = Beam->GetTimeOfFlight(0);
  if (PrintLevel>0) {
    cout << "|---- Kinematic constraints for strip ";
    cout.width(3); cout << StpID;
    cout << " ---------\n"
	 << "|     |   In    |   Out   |  Units  |\n"
	 << "| zr  |";
    cout.width(9);    cout << MinZ;   cout << "|";
    cout.width(9);    cout << MaxZ;   cout << "|";
    cout.width(9);    cout << "cm";   cout << "|\n";
    cout << "| tof |";
    cout.width(9);    cout << MinT;   cout << "|";
    cout.width(9);    cout << MaxT;   cout << "|";
    cout.width(9);    cout << "ns";   cout << "|\n";
    cout << "| Kb  |";
    cout.width(9);    cout << Kb_max; cout << "|";
    cout.width(9);    cout << Kb_min; cout << "|";
    cout.width(9);    cout << "MeV";   cout << "|\n";
    cout << "|--------------------------------------------------" << endl; 
  }
  //-------------------------------------------------------------------------------

  // Event for-loop
  for (int evt=0; evt<NEvents; evt++) {
    if (PrintLevel>0)
      cout << "\n***************** Event " << evt << "\n" << endl;
    
    TraceCan->cd(1);
    HCT->Draw();
    TraceCan->cd(2);
    HPT->Draw();
    
    // Reset the detector response
    for (int stp=0; stp<AnodeStps; stp++) 
      for (int col=0; col<AnodeCols+1; col++) {
	DeltaEB[stp][col] = 0;
	DeltaEL[stp][col] = 0;
	DeltaEH[stp][col] = 0;
      }
    
    // 1. Set beam inital conditions (beam energy, position)
    SetInitialKinematics(Kb_after_window);   

    // 2. Within the selected strip randomly select the position at
    // which the beam particle interacts with the target and calculate
    // the kinetic energy at the reaction point
    double zr = Rdm->Uniform(MinZ, MaxZ);
    double Kbr = Beam->GetFinalEnergy(0, Kb_after_window, zr, 1E-3);
    double TOF = Beam->GetTimeOfFlight(0);
    if (PrintLevel>0)
      cout << "Kbr = " << Kbr << "  zr = " << zr << "  tof = " << TOF << endl;
    
    // 3. Set the kinematics of all particles at the reaction point
    int ReacAllowed = SetReactionKinematics(Kbr, zr, TOF);
    if (ReacAllowed==0) {
      cout << "Warninig: reaction energetically not allowed for event " << evt 
	   << " (Kbr= " << Kbr << " MeV)." << endl;
      continue;
    }

    // 4. Propagate the beam particle (backwards in time) from the
    // reaction point to the entrance of MUSIC
    DeltaEB = PropagateParticle(Beam, evt, MaxTime, -UserDT);

    // 5. Propagate heavy particle (or decay daughters) and calculate energy 
    // loss in the anode elements
    if (DeDau1 && DeDau2) {
      DeltaED1 = PropagateParticle(DeDau1, evt, MaxTime, UserDT);
      DeltaED2 = PropagateParticle(DeDau2, evt, MaxTime, UserDT);
    }
    else
      DeltaEH = PropagateParticle(Heavy, evt, MaxTime, UserDT);
    
    // 6. Propagate light particle and calculate energy loss in the
    // anode elements
    DeltaEL = PropagateParticle(Light, evt, MaxTime, UserDT);

    // 7. Compute detector response (i.e. DE for beam + light + heavy)
    // Clone the particle trajectories
    ComputeDetectorResponse(evt);
    
    // 8. Display trace and particle trajecories   
    UpdateVisuals(evt, Kbr, zr, TOF, Wait);
 
    NTraces++;
  }
#endif
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Private method, to be used in an event loop.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::UpdateVisuals(int evt, double Kbr, double zr, double TOF, int Wait)
{
#if 1
  // 3D stuff
  short C,S,W;
  if (Light!=0 && Light->SaveTrajectory) {
    Light->GetTrajectoryAtt(C,S,W);
    TrajL[evt]->SetLineColor(C);
    TrajL[evt]->SetLineStyle(S);
    TrajL[evt]->SetLineWidth(W);
    Eve->AddElement(TrajL[evt]);
  }
  if (Heavy!=0 && Heavy->SaveTrajectory) {
    Heavy->GetTrajectoryAtt(C,S,W);
    TrajH[evt]->SetLineColor(C);
    TrajH[evt]->SetLineStyle(S);
    TrajH[evt]->SetLineWidth(W);
    Eve->AddElement(TrajH[evt]);
  } 
  if (DeDau1 && DeDau1->SaveTrajectory) {
    DeDau1->GetTrajectoryAtt(C,S,W);
    TrajD1[evt]->SetLineColor(C);
    TrajD1[evt]->SetLineStyle(S);
    TrajD1[evt]->SetLineWidth(W);
    Eve->AddElement(TrajD1[evt]);
  } 
  if (DeDau2 && DeDau2->SaveTrajectory) {
    DeDau2->GetTrajectoryAtt(C,S,W);
    TrajD2[evt]->SetLineColor(C);
    TrajD2[evt]->SetLineStyle(S);
    TrajD2[evt]->SetLineWidth(W);
    Eve->AddElement(TrajD2[evt]);
  } 
  Eve->Redraw3D();

  // 2D stuff
  // Traces of particles' energy loss and total.
  TraceCan->cd(1);
  // Legend
  if (LegPart->GetNRows()==0) {
    LegPart->AddEntry(Trace[evt][AnodeCols], "All particles", "l");
    if (DeDau1 && DeDau2) {
      LegPart->AddEntry(TraceD1[evt][AnodeCols], Form("%s",DeDau1->Name.c_str()), "l");
      LegPart->AddEntry(TraceD2[evt][AnodeCols], Form("%s",DeDau2->Name.c_str()), "l");
    }
    else
      LegPart->AddEntry(TraceH[evt][AnodeCols], Form("%s",Heavy->Name.c_str()), "l");
    LegPart->AddEntry(TraceL[evt][AnodeCols], Form("%s",Light->Name.c_str()), "l");
  }
  LegPart->Draw();
  // Update the kinematics label and then draw it
  LabelKine->Clear();
  LabelKine->AddText("Kinematics");
  LabelKine->AddLine(0.0,0.76,1.0,0.76);
  LabelKine->AddText(Form("beam: K=%.2f MeV  z_{r}=%.2f cm  tof=%.1f ns", Kbr, zr, TOF));
  if (DeDau1 && DeDau2) {
    LabelKine->AddText(Form("%s: K=%.2f MeV  #theta_{lab}=%.1f deg  #phi_{lab}=%.1f deg",
			    DeDau1->Name.c_str(), DeDau1->GetKE(), DeDau1->GetTheta()*180/pi,
			    DeDau1->GetPhi()*180/pi));
    LabelKine->AddText(Form("%s: K=%.2f MeV  #theta_{lab}=%.1f deg  #phi_{lab}=%.1f deg",
			    DeDau2->Name.c_str(), DeDau2->GetKE(), DeDau2->GetTheta()*180/pi,
			    DeDau2->GetPhi()*180/pi));
  }
  else 
    LabelKine->AddText(Form("%s: K=%.2f MeV  #theta_{lab}=%.1f deg  #phi_{lab}=%.1f deg",
			    Heavy->Name.c_str(), Heavy->GetKE(), Heavy->GetTheta()*180/pi,
			    Heavy->GetPhi()*180/pi));
  LabelKine->AddText(Form("%s: K=%.2f MeV  #theta_{lab}=%.1f deg  #phi_{lab}=%.1f deg",
			  Light->Name.c_str(), Light->GetKE(), Light->GetTheta()*180/pi,
			  Light->GetPhi()*180/pi));
  LabelKine->AddText(Form("#theta_{c.m.}=%.1f deg", theta_CM));

  LabelKine->Draw();
  // Draw traces
  if (DeDau1 && DeDau2) {
    TraceD1[evt][AnodeCols]->Draw("l same");
    TraceD2[evt][AnodeCols]->Draw("l same");
  }
  else
    TraceH[evt][AnodeCols]->Draw("l same");
  TraceL[evt][AnodeCols]->Draw("l same");
  Trace[evt][AnodeCols]->Draw("*l same");

  // Traces of energy loss in columns as a function of the strip
  // number.
  TraceCan->cd(2);
  TraceMult->GetYaxis()->SetRangeUser(0,AnodeCols+1);
  TraceMult->Draw();
#if 0
  for (int col=0; col<AnodeCols; col++)
    Trace[evt][col]->Draw("l same");
  Trace[evt][AnodeCols]->Draw("*l same");
  // Legend
  if (LegCol->GetNRows()==0) {
    LegCol->AddEntry(Trace[evt][AnodeCols],"All columns","l");
    for (int col=0; col<AnodeCols; col++)
      LegCol->AddEntry(Trace[evt][col], Form("Column %d", col),"l");
    // TraceCan->cd(3);
    // LegCol->Draw();
    TraceCan->cd(2);      
  }
  LegCol->Draw();
#endif

  TraceCan->Update();
  if (Wait==1)
    TraceCan->WaitPrimitive();

  // Remove the 3D trajecories. Make space for the trajectories of
  // the next event, but don't remove the ones of the last event.
  if (evt<NEvents-1) {
    if (Light!=0 && Light->SaveTrajectory)
      Eve->RemoveElement(TrajL[evt], (TEveElement*)Eve->GetCurrentEvent());
    if (Heavy!=0 && Heavy->SaveTrajectory)
      Eve->RemoveElement(TrajH[evt], (TEveElement*)Eve->GetCurrentEvent());
    if (DeDau1!=0 && DeDau1->SaveTrajectory)
      Eve->RemoveElement(TrajD1[evt], (TEveElement*)Eve->GetCurrentEvent());
    if (DeDau2!=0 && DeDau2->SaveTrajectory)
      Eve->RemoveElement(TrajD2[evt], (TEveElement*)Eve->GetCurrentEvent());
  }
#endif
  return;
}


///////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::WriteTraces(char* FileName)
{
#if 1
  TFile* Output = new TFile(FileName, "RECREATE");
  if (Trace!=0)
    for (int n=0; n<NTraces; n++)
      if (Trace[n]!=0) {
	for (int col=0; col<AnodeCols+1; col++) {
	  if (col==AnodeCols) 
	    Trace[n][col]->Write(Form("Trace%d",n), TObject::kOverwrite);
	  else
	    Trace[n][col]->Write(Form("Trace%dc%d",n,col), TObject::kOverwrite);
	}
      }
  Output->Close();
#endif
  return;
}



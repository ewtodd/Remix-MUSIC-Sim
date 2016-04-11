
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
  EneSigma = 0;
  EexcMax = 0;
  EexcMin = 0;
  Kb_after_window = 0;
  MaxParticles = 10;   // I have arbitrarily set the maximum number of particles equal to 10.
                       // Feel free to increase or decrease this limit.
  MaxPrevCS = 10;
  NFusedParticles = 0;
  NHeavyParticles = 0;
  NLightParticles = 0;
  NPrevCS = 0;
  NSegments = 0;
  NTraces = 0;
  ParamDirectory = "";

  NuF = new NuclideFinder();

  // Initialize selected pointers to zero (non-zero pointers initialized below).
  Beam = 0;
  Target = 0;
  Trace = 0;
  // Initialize some array of pointers.
  Fused = new Particle*[MaxParticles];
  Heavy = new Particle*[MaxParticles];
  Light = new Particle*[MaxParticles];
  FusedInTgt = new EnergyLoss*[MaxParticles];
  HeavyInTgt = new EnergyLoss*[MaxParticles];
  LightInTgt = new EnergyLoss*[MaxParticles];
  for (int i=0; i<MaxParticles; i++) {
    Fused[i] = 0;
    Light[i] = 0;
    Heavy[i] = 0;
    FusedInTgt[i] = 0;
    LightInTgt[i] = 0;
    HeavyInTgt[i] = 0;
  } 
  PrevCS = new TGraph*[MaxPrevCS];
  for (int i=0; i<MaxPrevCS; i++) 
    PrevCS[i] = 0;
  SegLength = 0;
  SegCMERange = 0;
  SegEexcRange = 0;

  // Geometry manager
  Geo = new TGeoManager("Geo", "MUSIC geometry manager");  

  // Define some materials
  MatVacuum = new TGeoMaterial("Vac", 0, 0, 0);
  MatAl = new TGeoMaterial("Al", 26.9815, 13, 2.7);
  MatSi = new TGeoMaterial("Si", 28.0855, 14, 2.329);
  // NOTE: Not sure about units of arguments

  // Define some media
  Vacuum = new TGeoMedium("Vacuum", 1, MatVacuum);
  Al = new TGeoMedium("Aluminium", 1, MatAl); 
  Si = new TGeoMedium("Silicon", 1, MatSi); 
  
  // Special media (empty, just using their address)
  CD2 = new TGeoMedium("CD2", 1, MatVacuum); 
  CF4 = new TGeoMedium("CF4", 1, MatVacuum); 
  D2 = new TGeoMedium("D2", 1, MatVacuum); 
  He3 = new TGeoMedium("He3", 1, MatVacuum); 
  He4 = new TGeoMedium("He4", 1, MatVacuum); 
  Kapton = new TGeoMedium("Kapton", 1, MatVacuum); 
  LiF = new TGeoMedium("LiF", 1, MatVacuum); 
  Mylar = new TGeoMedium("Mylar", 1, MatVacuum); 
  Ti = new TGeoMedium("Ti", 1, MatVacuum); 


  // Make the top container volume
  VolTop = Geo->MakeBox("VolTop", Vacuum, 300., 300., 300.);
  Geo->SetTopVolume(VolTop);


  // Zero other volume pointers
  VolAnode = 0;
  VolIC = 0;
  VolICBkFlange = 0;
  VolICFlange = 0;
  VolICPSGFrame = 0;
  VolICSec = 0;
  VolICWin = 0;
  VolRDBody = 0;
  VolRDFace = 0;
  VolRDSi1 = 0;
  VolRDSi2 = 0;
  VolRDSi3 = 0;
  VolSD = 0;
  VolSolDSDoor = 0;
  VolSolUSDoor = 0;
  VolTgt = 0;
  VolTgtFrame = 0;
  VolTgtWinDS = 0;
  VolTgtWinUS = 0;


  // Pseudo-random number generator.
  Rdm = new TRandom3();
  Rdm->SetSeed();     // Provide a seed that depends on the time.
}



///////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::CalculateCMEnergyRange()
{
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
  Kb = BeamInTgt->GetFinalEnergy(Kb,3.522, 0.001);

  if (BeamInTgt!=0) { 
    for (int i=0; i<NSegments; i++) 
      TotalLength += SegLength[i];
    
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
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////////
double* MUSIC_Simulator::CalculateELoss(Particle* P, EnergyLoss* PInTgt, int Event)
{
  double t, xi,yi,zi, xf,yf,zf;
  double x1,y1,z1, x2,y2,z2;
  double vx,vy,vz;
  double E,px,py,pz;
  double Q = P->Q;
  double m = P->Mass;
  double Ki = 0, K1 = 0, K2 = 0;
  double theta = P->GetTheta();
  double phi = P->GetPhi();
  double TOF, TrackLength;
  double* TrackLengthInSeg = new double[AnodeStps];
  double* DeltaE = new double[AnodeStps];

  for (int n=0; n<AnodeStps; n++)
    DeltaE[n] = 0;

  P->AllTraj[Event]->SetName(Form("%s evt %d", P->Name.c_str(),Event));
  
  // Get the initial conditions from the Particle object.
  //  cout << "\nTraj of " << P->Name << endl;
  P->GetX(t, xi, yi, zi);
  //  P->X.Print();
  P->GetP(E, px, py, pz);
  //  P->P.Print();
  Ki = E - m;

  //  cout << "Ki = " << Ki << " MeV   theta = " << theta*180/TMath::Pi() << " deg   phi = "
  //     << phi*180/TMath::Pi() << " deg" << endl;
  // Exit if the energy loss object has not been created or if the charge of the
  // particle is zero.
  if (PInTgt==0 || Q==0)
    return DeltaE;

  // Calculate the total track lenght (until the particle has nearly lost all its energy)
  TrackLength = PInTgt->GetPathLength(Ki, 0.005/*MeV*/, 1.0/*ns*/);
  if (TrackLength==0) 
    TrackLength = AnodeDepth;
  TOF = PInTgt->GetTimeOfFlight(Ki, TrackLength, 0.01);
  
  //  cout << "TrackLength = " << TrackLength << " cm   TOF = " << TOF << " ns" << endl;
  
  if (TrackLength>0) {
    // Get the final position (where the particle stops).
    xf = xi + TrackLength*sin(theta)*cos(phi);
    yf = yi + TrackLength*sin(theta)*sin(phi);
    zf = zi + TrackLength*cos(theta);
    // cout << "ri = " << xi << "," << yi << "," << zi << endl;
    // cout << "rf = " << xf << "," << yf << "," << zf << endl;
 
    // Check whether the final point is within the detector volume.
    if (fabs(xf)>AnodeLength/2 || fabs(yf)>AnodeHeight/2 || zf>AnodeDepth || zf<0)
      cout << "Out!" << endl;
    else 
      cout << "In!" << endl;    
    
    x1 = xi;
    y1 = yi;
    z1 = zi;
    K1 = Ki;
    z2 = 0;

    // Determine energy loss in each strip
    for (int n=0; n<AnodeStps; n++) {
      TrackLengthInSeg[n] = 0;
      z2 += AnodeDZ[n][0];

      if (z2<zi)
	continue;

      y2 = (yf-yi)*(z2-zi)/(zf-zi) + yi;
      x2 = (xf-xi)*(z2-zi)/(zf-zi) + xi;

      if (fabs(x2)>AnodeLength/2 || fabs(y2)>AnodeHeight/2)
	break;

      TrackLengthInSeg[n] = sqrt(pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2));
      K2 = PInTgt->GetFinalEnergy(K1, TrackLengthInSeg[n], TrackLengthInSeg[n]/100);
      if (K2<0) 
	break;

      DeltaE[n] = K1 - K2;
      //  cout << "S" << n << "  EL=" << K1-K2 << endl;
           
      x1 = x2;
      y1 = y2;
      z1 = z2;
      K1 = K2;
    }
    P->AllTraj[Event]->AddLine(xi,yi,zi, x1,y1,z1);
  }
  delete TrackLengthInSeg;
  return DeltaE;
}



///////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::CalculateExcEnergyRange()
{
  double pb, Eb;
  double Eexc_beg, Eexc_end;
  double mb = Beam->Mass;
  double mt = Target->Mass;
  double mf = Fused[0]->Mass;
  FourVector Pb("Pb");
  FourVector Pt("Pt", mt, 0, 0, 0);
  FourVector Ptot("Total four-mom. in the lab");
  double Kb = Kb_after_window;
  double Kb_min;
  float TotalLength = 0;
  
  if (BeamInTgt!=0) { 
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
    Kb_min = BeamInTgt->GetFinalEnergy(Kb, TotalLength, 0.001);
    pb = sqrt(2*mb*Kb_min*(1 + Kb_min/(2*mb)));
    Eb = sqrt(mb*mb + pb*pb);
    Pb.SetCoords(Eb, 0, 0, pb);
    Ptot = Pb + Pt;
    EexcMin = Eexc_end = sqrt(Ptot*Ptot) - mf;
    cout << "   Eexc(end) = " << Eexc_end << " MeV" << endl;
    cout << "Exc. energy range in each segment:" << endl;
    for (int i=0; i<AnodeStps; i++) {
      Kb = BeamInTgt->GetFinalEnergy(Kb, AnodeDZ[i][0], 0.001);
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
  }
  else {
    cout << "Warning: energy loss file for beam has not been loaded.  Calculations were not made.\n";
  }
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Show a 3D image of MUSIC with an given transparency level into the specified
// TEveManager object pointer.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::DrawMUSIC(TEveManager* gEve, short Transparency /*From 0 to 100*/)
{
  if (VolAnode!=0) {
    if (Transparency<0 || Transparency>100) {
      cout << "Warning: Transparency level must be from 0 to 100." << endl;
      Transparency = 0;
    }
    
    /* From: http://root.cern.ch/root/html/TGeoManager.html#TGeoManager:CloseGeometry
      Closing geometry implies checking the geometry validity, fixing shapes
      with negative parameters (run-time shapes) building the cache manager,
      voxelizing all volumes, counting the total number of physical nodes and
      registring the manager class to the browser.
    */
    Geo->CloseGeometry();
    
    TEveGeoTopNode* TopNode = new TEveGeoTopNode(Geo, Geo->GetTopNode());
    gEve->AddGlobalElement(TopNode);
    gEve->Redraw3D(kTRUE);
  }
  return;
}



///////////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::DrawTrajecotries(TEveManager* gEve)
{
  short C, W, S;
  TEveStraightLineSet* Traj;
  if (Beam!=0 && Beam->SaveTrajectory) {
    Beam->GetTrajectoryAtt(C,S,W);
    Traj = Beam->Trajectory;
    Traj->SetLineColor(C);
    Traj->SetLineStyle(S);
    Traj->SetLineWidth(W);
    gEve->AddElement(Traj);
  }
  
  if (Light[0]!=0 && Light[0]->SaveTrajectory) {
    Light[0]->GetTrajectoryAtt(C,S,W);
    for (int e=0; e<NEvents; e++) {
      TrajL[e]->SetLineColor(C);
      TrajL[e]->SetLineStyle(S);
      TrajL[e]->SetLineWidth(W);
      gEve->AddElement(TrajL[e]);
    }
  }
  if (Heavy[0]!=0 && Heavy[0]->SaveTrajectory) {
    Heavy[0]->GetTrajectoryAtt(C,S,W);
    for (int e=0; e<NEvents; e++) {
      TrajH[e]->SetLineColor(C);
      TrajH[e]->SetLineStyle(S);
      TrajH[e]->SetLineWidth(W);
      gEve->AddElement(TrajH[e]);
    }
  } 

  gEve->FullRedraw3D(kTRUE);
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Establish the dimensions of the MUSIC components (anode, cathode, etc).
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetAnode(string AnodeGeomFile, short Trans)
{

  ifstream GeomFile;
  string line;
  int line_counter = 0;

  AnodeDepth = 0;
  AnodeLength = 0;
  AnodeHeight = 0;
  AnodeStps = 0;
  AnodeCols = 0;

  GeomFile.open(AnodeGeomFile.c_str());
  if (!GeomFile.is_open()) 
    cout << "> ERROR: Anode geometry file \"" << AnodeGeomFile << "\" couldn't be opened." << endl;
  else {
    cout << "> Loading anode geometry loaded from \"" << AnodeGeomFile << "\"." << endl;

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
      for (int stp=0; stp<AnodeStps; stp++) {
	AnodeDX[stp] = new double[AnodeCols];
	AnodeDY[stp] = new double[AnodeCols];
	AnodeDZ[stp] = new double[AnodeCols];
	AnodeColor[stp] = new short[AnodeCols];
	AnodeSegName[stp] = new string[AnodeCols];
	for (int col=0; col<AnodeCols; col++) {
	  AnodeDX[stp][col] = 0;
	  AnodeDY[stp][col] = 0;
	  AnodeDZ[stp][col] = 0;
	  AnodeColor[stp][col] = 632; // red
	  AnodeSegName[stp][col] = "";
	}
      }
 
      // Get the name and dimensions for each anode segments(strip and
      // column)
      GeomFile.open(AnodeGeomFile.c_str());
      do {
	getline(GeomFile, line);
	if (line.find("Stp	Col	Name	Dx	Dy	Dz	Color	Comment")<line.npos)
	  break;
      } while (!GeomFile.eof());
      // First we need to count how many (not empty) lines are listed
      do {
	getline(GeomFile, line);
	if (!line.empty())
	  line_counter++;
      } while (!GeomFile.eof());
      GeomFile.close();
      //     cout << "Total lines: " << line_counter << endl;
      // Now that the number of lines has been established reopen the
      // text file to load the parameters
      GeomFile.open(AnodeGeomFile.c_str());
      do {
	getline(GeomFile, line);
	if (line.find("Stp	Col	Name	Dx	Dy	Dz	Color	Comment")<line.npos)
	  break;
      } while (!GeomFile.eof());
      // Loading parameters and printing them to confirm they have been
      // read correctly
      for (int nl=0; nl<line_counter; nl++) {
	int stp = 0;
	int col = 0;      
	string name;
	double dx, dy, dz;
	short color;
	GeomFile >> stp >> col >> name >> dx >> dy >> dz >> color;
	getline(GeomFile, line); // The last column is for comments
	AnodeSegName[stp][col] = name;
	AnodeDX[stp][col] = dx;
	AnodeDY[stp][col] = dy;
	AnodeDZ[stp][col] = dz;
	AnodeColor[stp][col] = color;
 	// cout << nl << "\t" << AnodeSegName[stp][col] << "\t" << AnodeDX[stp][col] << "\t" 
	//      << AnodeDY[stp][col] << "\t" << AnodeDZ[stp][col] << "\t" << AnodeColor[stp][col] 
	//      << endl;
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

      // Define anode volumes
      double z0 = 0;
      VolAnode = new TGeoVolume**[AnodeStps];
      for (int stp=0; stp<AnodeStps; stp++) {
	VolAnode[stp] = new TGeoVolume*[AnodeCols];
	z0 += AnodeDZ[stp][0]/2;
	double x0 = -AnodeLength/2;
	for (int col=0; col<AnodeCols; col++) {
	  if (AnodeDX[stp][col]>0) {
	    VolAnode[stp][col] = Geo->MakeBox(Form("VolAnode%d%d",stp,col), Vacuum, 
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
    } // end if (AnodeStps>0 && AnodeCols>0)
    
    HELoss = new TH2F("HELoss","HELoss", AnodeStps,-0.5, AnodeStps-0.5, 400,-4,6);
    HELoss->GetXaxis()->SetTitle("Strip number");
    HELoss->GetXaxis()->CenterTitle();
    HELoss->GetYaxis()->SetTitle("Energy loss - average beam energy loss (middle of strip) [MeV]");
    HELoss->GetYaxis()->CenterTitle(); 
  }
  
  cout << "Anode dimensions: " << AnodeLength << "x" << AnodeHeight << "x" 
       << AnodeDepth << "cm^3" << endl;
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Create the beam particle object.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetBeamParticle(string ParticleName, double M, int Q, int Color, double KineticE)
{
  Beam = new Particle(ParticleName, M, Q, true);
  Beam->SetTrajectoryAtt((short)Color);
  Kb_after_window = KineticE;
  return;
}


////////////////////////////////////////////////////////////////////////////
// Create new EnergyLoss objects for the corresponding particle type. As 
// with many other file names in this class, only the actual name is expected
// in this function's arguments, then the parameter directory will be
// inserted before the name.
////////////////////////////////////////////////////////////////////////////
bool MUSIC_Simulator::SetEnergyLossFile(string ParticleName, string TgtELossFile)
{
  bool MatchFound = 0;
  // Check if the particle name is "beam".
  if (ParticleName=="beam") {
    MatchFound = 1;
    if (TgtELossFile!="") {
      BeamInTgt = new EnergyLoss();
      BeamInTgt->LoadSRIMFile(ParamDirectory + TgtELossFile);
      BeamInTgt->SetIonMass(Beam->Mass);
      if (!BeamInTgt->GoodELossFile) {
	delete BeamInTgt;
	BeamInTgt = 0;
      }
    }
  }
  // If it doesn't match the string "beam" then look for a match with other particle type.
  else {
    // First look in the light particles.
    for (int p=0; p<NLightParticles; p++) 
      // If the name matches then create new the EnergyLoss objects and break.
      if (ParticleName==Light[p]->Name) {
	MatchFound = 1;
	if (TgtELossFile!="") {
	  LightInTgt[p] = new EnergyLoss();
	  LightInTgt[p]->LoadSRIMFile(ParamDirectory + TgtELossFile);
	  LightInTgt[p]->SetIonMass(Light[p]->Mass);
	  if (!LightInTgt[p]->GoodELossFile) {
	    delete LightInTgt[p];
	    LightInTgt[p] = 0;
	  }
	}
	break;
      } 
    // If a match has not been found, look in the heavy particles.
    if (!MatchFound) {
      for (int p=0; p<NHeavyParticles; p++)
	// If the name matches then create new the EnergyLoss objects and break.
	if (ParticleName==Heavy[p]->Name) {
	  MatchFound = 1;
	  if (TgtELossFile!="") {
	    HeavyInTgt[p] = new EnergyLoss();
	    HeavyInTgt[p]->LoadSRIMFile(ParamDirectory + TgtELossFile);
	    HeavyInTgt[p]->SetIonMass(Heavy[p]->Mass);
	    if (!HeavyInTgt[p]->GoodELossFile) {
	      delete HeavyInTgt[p];
	      HeavyInTgt[p] = 0;
	    }
	  } 
	  break;
	}
    }
    // Finally, look in the fused particles.
    if (!MatchFound) {
      for (int p=0; p<NFusedParticles; p++)
	// If the name matches then create new the EnergyLoss objects and break.
	if (ParticleName==Fused[p]->Name) {
	  MatchFound = 1;
	  if (TgtELossFile!="") {
	    FusedInTgt[p] = new EnergyLoss();
	    FusedInTgt[p]->LoadSRIMFile(ParamDirectory + TgtELossFile);
	    FusedInTgt[p]->SetIonMass(Fused[p]->Mass);
	    if (!FusedInTgt[p]->GoodELossFile) {
	      delete FusedInTgt[p];
	      FusedInTgt[p] = 0;
	    }
	  } 
	  break;
	}
    }
    // If no match whatsoever was found, then print warning message.
    if (!MatchFound) {
      cout << "> WARNING: Particle '" << ParticleName << "' was not found in the particles' list."  
	   << " Energy loss file(s) was (were) not loaded." << endl;
    }
  } // end of else (not beam particle).
  return MatchFound;
}


///////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetFusedParticle(string ParticleName, double M, int Q, int NEexc, double* Eexc)
{
  int p = NFusedParticles;
  if (p<MaxParticles) {
    Fused[p] = new Particle(ParticleName, M, Q);
    Fused[p]->SetExcEnergies(NEexc, Eexc);
    NFusedParticles++;
  }
  else
    cout << "> Warning: Maximum number of particles has been exceeded." << endl;
  return;
}


///////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetHeavyParticle(string ParticleName, double M, int Q, int Color, 
				       int NEexc, double* Eexc)
{
  int p = NHeavyParticles;
  if (p<MaxParticles) {
    Heavy[p] = new Particle(ParticleName, M, Q, true);
    Heavy[p]->SetTrajectoryAtt((short)Color);
    Heavy[p]->SetExcEnergies(NEexc, Eexc);
    NHeavyParticles++;
  }
  else
    cout << "> Warning: Maximum number of particles has been exceeded." << endl;
 return;
}


///////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetLightParticle(string ParticleName, double M, int Q, int Color)
{
  int p = NLightParticles;
  if (p<MaxParticles) {
    Light[p] = new Particle(ParticleName, M, Q, true);
    Light[p]->SetTrajectoryAtt((short)Color);
    NLightParticles++;
  }
  else
    cout << "> Warning: Maximum number of particles has been exceeded." << endl;
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Especify the directory where the parameters (such as eneryg loss files) are
// located.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetParamDirectory(string Dir)
{
  ParamDirectory = Dir;
  return;
}


///////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetPreviousCrossSection(string FileName, string Format, short Marker,
					      short Color)
{
  int i = NPrevCS;
  if (i<MaxPrevCS) {
    FileName = ParamDirectory + FileName;
    PrevCS[i] = new TGraph(FileName.c_str(), Format.c_str());
    PrevCS[i]->SetMarkerStyle(Marker);
    PrevCS[i]->SetMarkerColor(Color);
    NPrevCS++;
  }
  return;
}



///////////////////////////////////////////////////////////////////////////////////
// Taken from PhysicalEventProcessor
// Assuming the reaction string (name) look like "20Ne(12C,32S,p)31P", where 
// 20Ne is the beam, 12C is the target, 32S is the fused particle, p is the light
// outgoing particle and 31P is the heavy recoil.
///////////////////////////////////////////////////////////////////////////////////
#if 0
void MUSIC_Simulator::SetPossibleReaction(string ReactionName)
{
  
  bool FoundLightPart = 0;
  int PosTarget, Kl_bins=0, Th_bins=0;
  string LightPName = "";
  float Kl_max=0, Kl_min=0, Kb_min=0, Kb_max=0, Th_max=0, Th_min=0, zr_min=0, zr_max=0;
  this->ReactionName[NReactions] = ReactionName;

  // Determine the which of the light particles was involved in this reaction.
  // 1st find the target name position
  PosTarget = ReactionName.find(Target->Name);
  // 2nd look for the light particle name after the target name position.
  for (int p=0; p<NLightParticles; p++)
    if (ReactionName.find(Light[p]->Name, PosTarget+1) < ReactionName.length()) {
      LightPName = Light[p]->Name;
      ml[NReactions] = Light[p]->Mass;
      FoundLightPart = 1;
      IsParticleInReaction[p][NReactions] = 1;
      Kl_bins = HP_Kl_bins[p];
      Kl_max = HP_Kl_max[p];
      Kl_min = HP_Kl_min[p];
      Th_bins = HP_Theta_bins[p];
      Th_max = HP_Theta_max[p];
      Th_min = HP_Theta_min[p];
    }
  if (FoundLightPart) {   
    // Create the histograms for this reaction.    
    ReactionReady[NReactions] = 1;
    cout << "> Reaction " << ReactionName << " ready!" << endl;
    NReactions++;
  }
  else 
  cout << "> WARNING: Could not find light particle for reaction " << ReactionName << endl;
  return;
}
#endif

///////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetSegmentLength(int NSegments, float* SegLength /*cm*/)
{
  this->NSegments = NSegments;
  this->SegLength = new double[NSegments];
  SegCMERange = new double[NSegments];
  SegEexcRange = new double[NSegments];
  for (int i=0; i<NSegments; i++) {
    this->SegLength[i] = SegLength[i];
    SegCMERange[i] = 0;
    SegEexcRange[i] = 0;
  }
  HELoss = new TH2F("HELoss","HELoss", NSegments,-0.5, NSegments-0.5, 400,0,6);
  HELoss->GetXaxis()->SetTitle("Segment number");
  HELoss->GetXaxis()->CenterTitle();
  HELoss->GetYaxis()->SetTitle("Energy loss [MeV]");
  HELoss->GetYaxis()->CenterTitle();
  return;
}


///////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetStripEnergyResolution(float Sigma)
{
  EneSigma = Sigma;
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Create the target particle object (currently the charge is not used).
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetTargetParticle(string ParticleName, double M, int Q)
{ 
  Target = new Particle(ParticleName, M, Q);
  return;
}

///////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::ShowCMEnergyRange()
{
  double CME_beg = CMEMax, CME_end = CMEMin, DeltaCME;
  TH2F* Bkgnd;
  TCanvas* Can;
  TLine** RangeInSeg;

  if (BeamInTgt!=0) { 
    RangeInSeg = new TLine*[AnodeStps];
    for (int i=0; i<AnodeStps; i++) {
      RangeInSeg[i] = new TLine();
      RangeInSeg[i]->SetLineColor(kGray+1);
      RangeInSeg[i]->SetLineWidth(3);
    }

    // Change in exc. energy
    DeltaCME = CMEMax - CMEMin;

    RangeInSeg[0]->SetX2(CME_beg);
    RangeInSeg[0]->SetY2(0);
    Bkgnd = new TH2F("Bkgnd","C.M. Energy Range",
		     10,CME_end-0.05*DeltaCME,CME_beg+0.05*DeltaCME,
		     AnodeStps,-0.5, AnodeStps+0.5);
    Bkgnd->GetXaxis()->SetTitle("C.M. Energy [MeV]");
    Bkgnd->GetXaxis()->CenterTitle();
    Bkgnd->GetYaxis()->SetTitle("Strip number");
    Bkgnd->GetYaxis()->CenterTitle();

    for (int i=0; i<AnodeStps; i++) {
      CME_end = CME_beg - SegCMERange[i];
      RangeInSeg[i]->SetX1(CME_end);
      RangeInSeg[i]->SetY1(i);
      // The excitation energy at the end of this segment is the excitation energy at the beginnig
      // of the next segment.
      if (i+1<AnodeStps) {
	CME_beg = CME_end;
	RangeInSeg[i+1]->SetX2(CME_beg);
	RangeInSeg[i+1]->SetY2(i+1);
      }
    }
    Can = new TCanvas("Can","CM energy range",1400,900);
    Can->SetGrid();
    Bkgnd->Draw();
    for (int i=0; i<AnodeStps; i++)
      RangeInSeg[i]->Draw();
  }
  else {
    cout << "Warning: energy loss file for beam has not been loaded.  Calculations were not made.\n";
  }
  return;
}


///////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::ShowPreviousCrossSection(float XMin, float XMax, float YMin, float YMax)
{
  TH2F* Bkgnd;
  TCanvas* Can;
  double CurrentCMEMax, CurrentCMEMin; 
  TPolyLine** CMERegion;

  Bkgnd = new TH2F("Bkgnd2", 
		   Form("CME range for each segment (K_{b}=%.2f MeV, P=X Torr)",Kb_after_window),
		   10,XMin,XMax, 10,YMin,YMax);
  Bkgnd->GetXaxis()->SetTitle("E_{CM} [MeV]");
  Bkgnd->GetXaxis()->CenterTitle();
  Bkgnd->GetYaxis()->SetTitle("Cross section [mb]");
  Bkgnd->GetYaxis()->CenterTitle();
   
  CMERegion = new TPolyLine*[AnodeStps];
  CurrentCMEMax = CMEMax;
  for (int i=0; i<AnodeStps; i++) {
    CurrentCMEMin = CurrentCMEMax - SegCMERange[i];
    CMERegion[i] = new TPolyLine(5);    
    if (CurrentCMEMin>0) {      
      CMERegion[i]->SetPoint(0, CurrentCMEMax, YMin);
      CMERegion[i]->SetPoint(1, CurrentCMEMax, YMax);
      CMERegion[i]->SetPoint(2, CurrentCMEMin, YMax);
      CMERegion[i]->SetPoint(3, CurrentCMEMin, YMin);
      CMERegion[i]->SetPoint(4, CurrentCMEMax, YMin);
      CMERegion[i]->SetFillStyle(3002);
      CMERegion[i]->SetFillColor(kRed - 5 - 2*(i%2));
    }
    CurrentCMEMax = CurrentCMEMin;
  }

  Can = new TCanvas("Can1","Can1",1200,800);
  Can->SetGrid();
  Bkgnd->Draw();
  for (int i=0; i<AnodeStps; i++)
    CMERegion[i]->Draw("f");
  for (int i=0; i<NPrevCS; i++)
    PrevCS[i]->Draw("p same");
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::Simulate(int Reaction, int SegNum, int NEvents)
{
  this->NEvents = NEvents;
  double Kbr, Klr, theta_CM, phi_CM, pf_CM;
  double Eb, pb, theta_b, phi_b;
  double El, pl, theta_l, phi_l;
  double BetaX, BetaY, BetaZ;
  FourVector Ptot;
  double mb = Beam->Mass;
  double mt = Target->Mass;
  double mf = Fused[0]->Mass;
  double ml = Light[Reaction]->Mass;
  double mh = Heavy[Reaction]->Mass;
  double Kb_min, Kb_max, MinZ, MaxZ, zr, Ex=0, Mh;
  FourVector Pb("Pb");
  FourVector Pt("Pt");
  FourVector Pf("Pf");
  FourVector Pl("Pl");
  FourVector Ph("Ph");
  TrajH = new TEveStraightLineSet*[NEvents];
  TrajL = new TEveStraightLineSet*[NEvents];
  double* DeltaEB = new double[AnodeStps];
  double* DeltaEB_mid = new double[AnodeStps];
  double* DeltaEL;
  double* DeltaEH;
  double DeltaE;
  double zb1, zb2, Kb1, Kb2;
  TF1* Gaussian = 0;
  if (EneSigma!=0) 
    Gaussian = new TF1("Gaussian","gaus",0, 100);

  cout << "Simulating reactions ... " << endl;
  cout << "Gas volume = " << AnodeLength << "x"<< AnodeHeight << "x"<< AnodeDepth << endl;
  
  TCanvas* Can = new TCanvas("Can","Traces", 0, 0, 1918, 630);
  HELoss->Draw();

  Trace = new TGraph*[NEvents];

  Kb1 = Kb_after_window;
  zb1 = 0; 
  zb2 = 0;
  // Get the beam energy loss for all strips. This quantity will be
  // subtracted from the total energy loss of the residual particles
  // (as it is done for the data analysis).
  for (int n=0; n<AnodeStps; n++) {
    zb2 += AnodeDZ[n][0];
    DeltaEB_mid[n] = 0;
    Kb2 = BeamInTgt->GetFinalEnergy(Kb1, zb2-zb1, 0.01);
    DeltaEB_mid[n] = Kb1 - Kb2;
    cout << "Beam  stp=" << AnodeSegName[n][0] << "  DeltaEB_mid=" << DeltaEB_mid[n] 
	 << " MeV  dz=" << AnodeDZ[n][0] << " cm" << endl;
    zb1 = zb2;
    Kb1 = Kb2;
  }

  
  // Get the beam energy limits in the selected strip.
  MinZ = 0;
  MaxZ = AnodeDZ[0][0];
  for (int n=1; n<=SegNum; n++) {
    MinZ += AnodeDZ[n-1][0];
    MaxZ += AnodeDZ[n][0];
  }
  cout << "\nMinZ = " << MinZ << "    MaxZ = "<< MaxZ << endl;
  Kb_max = BeamInTgt->GetFinalEnergy(Kb_after_window, MinZ, 0.01);
  Kb_min = BeamInTgt->GetFinalEnergy(Kb_after_window, MaxZ, 0.01);
  cout << "MaxKb = " << Kb_max << "    MinKb = "<< Kb_min << endl;

  for (int e=0; e<NEvents; e++) {
    //    cout << "\n***************** Event " << e << "\n" << endl;
   
    // Within the selected strip select randomly the energy at which the beam particle
    // interacts with the target and from this value calculate the reaction z.
    Kbr = Rdm->Uniform(Kb_min, Kb_max);
    zr = BeamInTgt->GetPathLength(Kb_after_window, Kbr, 0.01/*ns*/);
    cout << "Kbr = " << Kbr << "  zr = " << zr << endl;
    
    Kb1 = Kb_after_window;
    Kb1 = BeamInTgt->GetFinalEnergy(Kb1,3.522, 0.001);
    zb1 = 0; 
    zb2 = 0;
    bool skip = 0;
    for (int n=0; n<AnodeStps; n++) {
      zb2 += AnodeDZ[n][0];
      DeltaEB[n] = 0;
      if (zb2<zr)
	Kb2 = BeamInTgt->GetFinalEnergy(Kb1, zb2-zb1, 0.01);
      else if (!skip) {
	Kb2 = BeamInTgt->GetFinalEnergy(Kb1, zr-zb1, 0.01);
	DeltaEB[n] = Kb1 - Kb2;
	//	cout << "Beam S=" << n << "  DeltaE = " << DeltaEB[n] << endl;  
	skip = 1;
      }
      if (!skip) { 
	DeltaEB[n] = Kb1 - Kb2;
	//	cout << "Beam S=" << n << "  DeltaE = " << DeltaEB[n] << endl;  
	zb1 = zb2;
	Kb1 = Kb2;
      }
    }
    
    // Linear momentum and total energy of the beam particle in the lab.
    pb = sqrt(2*mb*Kbr*(1 + Kbr/(2*mb)));
    theta_b = 0;
    phi_b = 0;
    Eb = sqrt(mb*mb + pb*pb);

    // Center of mass beta (v/c)
    BetaX = pb*sin(theta_b)*cos(phi_b)/(Eb+mt);
    BetaY = pb*sin(theta_b)*sin(phi_b)/(Eb+mt);
    BetaZ = pb*cos(theta_b)/(Eb+mt);
    //    cout << "BetaX=" << BetaX << "  BetaY=" << BetaY << "  BetaZ=" << BetaZ << endl; 

    // Four-momentum of the beam.
    Pb.SetCoords(Eb, pb*sin(theta_b)*cos(phi_b), pb*sin(theta_b)*sin(phi_b), pb*cos(theta_b));
    Beam->SetP(Pb);
    Beam->SetX(0, 0, 0, zr);
    //cout << "\nBOOST to CM!" << endl;
    //Beam->Print();

    // Four-momentum of the target.
    Pt.SetCoords(mt, 0, 0, 0);
    Target->SetP(Pt);    
    Target->SetX(0, 0, 0, zr);       
    //Target->Print();

    //cout << "\nReturn to Lab!" << endl;
    //Beam->Boost(0,0,-BetaZ);
    //Beam->Print();
    // Target->Boost(0,0,-BetaZ);
    //Target->Print();


    // Total four momentum in the lab.
    Ptot = Pb + Pt;
    //   Ptot.SetName("Total four-mom. in the lab");
  
    // Exc. energy of fused particle.
    //    cout << "Eexc(21Ne) = " << sqrt(Ptot*Ptot) - mf << " MeV" << endl;

    // Randomly select the scattering angle in the center of mass, then go to the laboratory
    // reference frame.
    theta_CM = acos(Rdm->Uniform(-1,1));
    phi_CM = Rdm->Uniform(-pi,pi);
    //   cout << "theta_cm=" << theta_CM*180/pi << "   phi_cm=" << phi_CM*180/pi << endl;

    // Select randomly the energy of excitation of the heavy particle from 0 to the 
    // maximum possible value (skip in case of elastic scattering).
    //cout << "(mb+mt) = " << mb+mt << "\n(ml+mh) = " << ml+mh << endl;
    //  if ((mb+mt)!=(mh+ml)) 
    //Ex = Rdm->Uniform(0, sqrt(Ptot*Ptot) - ml - mh);
    //else 
    Ex = 0;
    Mh = mh + Ex;

    // Final momentum in the CM reference frame (Ptot*Ptot is an invariant quantity).
    pf_CM = sqrt((Ptot*Ptot - pow(ml+Mh,2))*(Ptot*Ptot - pow(ml-Mh,2))/(4*(Ptot*Ptot)));
    //   cout << "Ex = " << Ex << " MeV    pf_CM = " << pf_CM << endl;
    
    // Set the four-momentum components of the light particle in the center of mass.
    Pl.SetCoords(sqrt(ml*ml + pf_CM*pf_CM), -pf_CM*sin(theta_CM)*cos(phi_CM),
		 -pf_CM*sin(theta_CM)*sin(phi_CM), -pf_CM*cos(theta_CM));
    // Do a Lorentz transformation (boost) into the lab reference frame.
    // I've double checked the sign of the boost and is correct (-Beta).
    Pl.Boost(-BetaX, -BetaY, -BetaZ);
    Light[Reaction]->SetP(Pl);
    // Initial position of the light particle (at the target).
    Light[Reaction]->SetX(0, 0, 0, zr);
    //Light[Reaction]->Print();

    // Four-momentum of the heavy recoil (lab).
    Ph = Pb + Pt - Pl;
    Heavy[Reaction]->SetP(Ph);
    Heavy[Reaction]->SetX(0, 0, 0, zr);
    //Heavy[Reaction]->Print();
    
    DeltaEH = CalculateELoss(Heavy[Reaction], HeavyInTgt[Reaction], e);
    DeltaEL = CalculateELoss(Light[Reaction], LightInTgt[Reaction], e);
    if (e<Particle::MaxEvents) {
      TrajH[e] = (TEveStraightLineSet*)Heavy[Reaction]->AllTraj[e]->Clone();
      TrajL[e] = (TEveStraightLineSet*)Light[Reaction]->AllTraj[e]->Clone();
    }
    float total=0;
    float average=0;

    Trace[e] = new TGraph();
    Trace[e]->SetName(Form("evt %d",e));
    for (int n=0; n<AnodeStps; n++) {

      // Previously we were normalizing to the energy loss of the beam
      // minus the energy loss in the the dead layer, which in this
      // version of the code is typically strip 1.
      //  DeltaE = DeltaEB[n] + DeltaEL[n] + DeltaEH[n] - (DeltaEB_mid[n]-DeltaEB_mid[1]);

      // In this version, we normalize to the average energy loss of
      // the beam in the middle of the strip.
      DeltaE = DeltaEB[n] + DeltaEL[n] + DeltaEH[n] - DeltaEB_mid[n];

      //cout<< "Energy loss L: "<<DeltaEL[n]<<  "Energy loss H: "<<DeltaEH[n]<<endl;
      if (EneSigma!=0 && Gaussian!=0 && DeltaE>0) {
	Gaussian->SetRange(0.0, 2*DeltaE);
	Gaussian->SetParameters(1.0, DeltaE, EneSigma);
	DeltaE = Gaussian->GetRandom();
      }
      
      Trace[e]->SetPoint(n, n, DeltaE);
      if(n>=SegNum && n<=(SegNum +6)) total=DeltaE+total;  
    }
   
     
   

    average=total/7;
    //    cout<<"Average= "<<average<<endl;
    Trace[e]->Draw("*l same");
    NTraces++;
    Can->Update();
    //Can->WaitPrimitive();
  }
  

  return;
}



///////////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::WriteTraces(char* FileName)
{
  TFile* Output = new TFile(FileName, "RECREATE");
  if (Trace!=0)
    for (int n=0; n<NTraces; n++)
      Trace[n]->Write(Form("Trace%d",n), TObject::kOverwrite);
  return;
}



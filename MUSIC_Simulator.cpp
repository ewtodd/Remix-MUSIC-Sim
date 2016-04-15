
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
  cout << "| ver 1.1 (2016/4)                                                             |" << endl;
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
  SegLength = 0;
  SegCMERange = 0;
  SegEexcRange = 0;

  // Initialize selected pointers to zero (non-zero pointers initialized below).
  Beam = 0;
  Target = 0;
  Trace = 0;
  Compound = 0;
  Light = 0;
  Heavy = 0;

  // Nuclide finder object
  NuF = new NuclideFinder();

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
  Kb = Beam->GetFinalEnergy(0, Kb,3.522, 0.001);

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
#endif
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Calculate and array with the beam energy loss for all strips and return it. 
// This quantity will be subtracted from the total energy loss of the residual 
// particles (as it is done in the data analysis). This method assumes the beam
// momentum is parallel to the z-axis and it does not account for the anode
// x splitting.
///////////////////////////////////////////////////////////////////////////////////
double* MUSIC_Simulator::CalcAverageBeamELoss()
{
  double* DeltaEB_ave = new double[AnodeStps];  // average beam energy loss
  double Kb1 = Kb_after_window;
  double Kb2 = 0;
  double zb1 = 0; 
  double zb2 = 0;
  for (int n=0; n<AnodeStps; n++) {
    zb2 += AnodeDZ[n][0];
    DeltaEB_ave[n] = 0;
    // Calculate the finial energy starting from Kb1 for a distance of
    // (zb2-zb1) in 100 steps.
    Kb2 = Beam->GetFinalEnergy(0, Kb1, zb2-zb1, (zb2-zb1)/100);
    DeltaEB_ave[n] = Kb1 - Kb2;
    zb1 = zb2;
    Kb1 = Kb2;
  }
  return DeltaEB_ave;
}


///////////////////////////////////////////////////////////////////////////////////
// 
///////////////////////////////////////////////////////////////////////////////////
double* MUSIC_Simulator::CalculateELoss(Particle* P, int Event)
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
  if (Q==0)
    return DeltaE;

  // Calculate the total track lenght (until the particle has nearly lost all its energy)
  TrackLength = P->GetPathLength(0, Ki, 0.005/*MeV*/, 1.0/*ns*/);
  if (TrackLength==0) 
    TrackLength = AnodeDepth;
  TOF = P->GetTimeOfFlight(0, Ki/*MeV*/, TrackLength/*cm*/, 0.01/*cm*/);
  
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
      K2 = P->GetFinalEnergy(0, K1, TrackLengthInSeg[n], TrackLengthInSeg[n]/100);
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
  // if (Beam!=0 && Beam->SaveTrajectory) {
  //   Beam->GetTrajectoryAtt(C,S,W);
  //   for (int e=0; e<NEvents; e++) {
  //     TrajB[e]->SetLineColor(C);
  //     TrajB[e]->SetLineStyle(S);
  //     TrajB[e]->SetLineWidth(W);
  //     gEve->AddElement(TrajB[e]);
  //   }
  // }
  
  if (Light!=0 && Light->SaveTrajectory) {
    Light->GetTrajectoryAtt(C,S,W);
    for (int e=0; e<NEvents; e++) {
      TrajL[e]->SetLineColor(C);
      TrajL[e]->SetLineStyle(S);
      TrajL[e]->SetLineWidth(W);
      gEve->AddElement(TrajL[e]);
    }
  }
  if (Heavy!=0 && Heavy->SaveTrajectory) {
    Heavy->GetTrajectoryAtt(C,S,W);
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
// This method transports the particle from its initial position until the maximum
// transport time has been reached or if the particle stops or if it leaves the 
// detector volume. As the particle is propagated through the gaseous medium, 
// energy loses in each anode segment are saved in a bidimensional array (DE) which 
// is retured when this method ends.
///////////////////////////////////////////////////////////////////////////////////
double** MUSIC_Simulator::PropagateParticle(Particle* PO, int Event, double MaxTime, double UserDT)
{
  double** DE = new double*[AnodeStps];
  for (int stp = 0; stp<AnodeStps; stp++) {
    DE[stp] = new double[AnodeCols];
    for (int col = 0; col<AnodeCols; col++) 
      DE[stp][col] = 0;
  }
  TGeoVolume* Vol = 0;
  bool Skip = 0;
  int step = 0;
  // Get the initial conditions from the Particle object.
  double Z = PO->Z;
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
    double vx = c*p_mag*cos(phi)*sin(theta)/Ene;
    double vy = c*p_mag*sin(phi)*sin(theta)/Ene;
    double vz = c*p_mag*cos(theta)/Ene;

    // Very short time step for dense media, equivalent to distance
    // steps of 1 um. Let's just hope that the numerator is never
    // zero.
    double Dt = c*m*1E-4/p_mag;   // from p/c = m*dx/dt -> dt = c*m*dx/p
    double tf = ti + Dt;
    double xf = xi + vx*Dt;
    double yf = yi + vy*Dt;
    double zf = zi + vz*Dt;

    cout << step << " (" << tf << ", " << xf << " " << yf << " " << zf << ")  Dt=" << Dt;

    // Exit the while loop if the particl leaves the anode volume.
    if (zf>AnodeDepth || zf<0 || xf>AnodeLength/2 || xf<-AnodeLength/2 || 
	yf>AnodeHeight/2 || yf<-AnodeHeight/2) {
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
      cout << PO->Name<< ": anode segment not found." << endl;
      break;
    }

    // Get the energy loss for the initial energy (Ki) in the gas
    // medium (0) over a small (differential) path length (dist).
    double dist = sqrt(pow(xf-xi,2) + pow(yf-yi,2) + pow(zf-zi,2));
    double Kf = Ki - PO->GetEnergyLoss(0, Ki, dist);
    // Exit the while loop if the particle has stopped.
    if (Kf<=0) {
      cout << "Less than ZERO! " << Kf << endl;
      break;
    }

    DE[stp][col] += Ki - Kf;

    cout << "\t d=" << dist << " \tKf=" << Kf << " \tDE=" << DE[stp][col] << endl; 

    // Get the momentum magnitude using the new (lower) kinetic energy
    // (I'm using a relativistic formula, although it is unlikely for
    // us to use relativistic energies).
    p_mag = sqrt(2*m*Kf*(1+Kf/2/m));
    // Reduce the total energy by amount of energy deposited in the
    // medium.
    Ene -= DE[stp][col];

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
  }

  return DE;
}


///////////////////////////////////////////////////////////////////////////////////
// Establish the dimensions of the MUSIC anode segments by reading a text file with
// the following structre
//------------------------------------------------------------------------------|
// Title                                                                        |
// Number of strips: #                                                          |
// Number of columns: #                                                         |
// Stp	Col	Name	Dx	Dy	Dz	Color	Comment                 |
// #	#	string	#f	#f	#f	#	string (spaces allowed) |
// ... (until all segments have been specified)                                 |
//------------------------------------------------------------------------------|
// , where '#' is meant to be replaced by an integer and '#f' by a floating point 
// number.
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
void MUSIC_Simulator::SetBeamParticle(string Name, int Color, string ELossFile, double K)
{
  // Use the nuclide finder object to determine the mass and atomic
  // number of this particle. Mass must be in MeV/c^2) and Z (in e).
  double m = NuF->GetMass(Name, "MeV/c^2");
  int Z = NuF->GetZ(Name);
  Beam = new Particle(Name, m, Z, true /*SaveTrajectories on*/);
  Beam->SetTrajectoryAtt((short)Color);
  // Currently, this simulation is restricted to one medium (gas) in MUSIC.
  Beam->SetMedium(ELossFile);
  Kb_after_window = K;
  double E = K + m;

  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Create the fused (a.k.a. compound) particle object.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetCompoundParticle(string Name)
{
  // Use the nuclide finder object to determine the mass and atomic
  // number of this particle. Mass must be in MeV/c^2) and Z (in e).
  double m = NuF->GetMass(Name, "MeV/c^2");
  int Z = NuF->GetZ(Name);
  Compound = new Particle(Name, m, Z, true /*SaveTrajectories on*/);
  //  Compound->SetTrajectoryAtt((short)Color);
  // Need excitation energy based on kinematics

  //  Compound->SetExcEnergies(NEexc, Eexc);
  // I would be good to print the excitation energy of the compound 
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Create the heavy particle (evaporation residue) object.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetHeavyParticle(string Name, int Color, string ELossFile, int NEexc,
				       double* Eexc)
{
  // Use the nuclide finder object to determine the mass and atomic
  // number of this particle. Mass must be in MeV/c^2) and Z (in e).
  double m = NuF->GetMass(Name, "MeV/c^2");
  int Z = NuF->GetZ(Name);
  Heavy = new Particle(Name, m, Z, true /*SaveTrajectories on*/);
  Heavy->SetTrajectoryAtt((short)Color);
  // Currently, this simulation is restricted to one medium (gas) in MUSIC.
  Heavy->SetMedium(ELossFile);
  Heavy->SetExcEnergies(NEexc, Eexc);
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Establish the kinematics of the particles right after the window.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetInitialKinematics(double Kbi/*MeV*/)
{
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
  Beam->Print();
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Create the light particle (evaporation residue) object. It is assumed that the 
// light particle (e.g. p, n, alpha) is in its ground state.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetLightParticle(string Name, int Color, string ELossFile)
{
  // Use the nuclide finder object to determine the mass and atomic
  // number of this particle. Mass must be in MeV/c^2) and Z (in e).
  double m = NuF->GetMass(Name, "MeV/c^2");
  int Z = NuF->GetZ(Name);
  Light = new Particle(Name, m, Z, true /*SaveTrajectories on*/);
  Light->SetTrajectoryAtt((short)Color);
  // Currently, this simulation is restricted to one medium (gas) in MUSIC.
  Light->SetMedium(ELossFile);
  return;
}




///////////////////////////////////////////////////////////////////////////////////
// Establish the kinematics of the particles at the reaction point.
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetReactionKinematics(double Kbr/*MeV*/, double zr/*cm*/)
{
  cout << "\nReaction kine" << endl;
  double mb = Beam->Mass;
  double mt = Target->Mass;
  double mc = Compound->Mass;
  double ml = Light->Mass;
  double mh = Heavy->Mass;

  // Linear momentum and total energy of the beam particle in the lab.
  double pb = sqrt(2*mb*Kbr*(1 + Kbr/2/mb));
  double theta_b = 0;
  double phi_b = 0;
  double Eb = sqrt(mb*mb + pb*pb);

  // Center of mass beta (v/c)
  double BetaX = pb*sin(theta_b)*cos(phi_b)/(Eb+mt);
  double BetaY = pb*sin(theta_b)*sin(phi_b)/(Eb+mt);
  double BetaZ = pb*cos(theta_b)/(Eb+mt);
  cout << "BetaX=" << BetaX << "  BetaY=" << BetaY << "  BetaZ=" << BetaZ << endl; 

  // Four-momentum of the beam.
  Beam->SetP(Eb, pb*sin(theta_b)*cos(phi_b), pb*sin(theta_b)*sin(phi_b), pb*cos(theta_b));
  Beam->SetX(0, 0, 0, zr);
  Beam->Print();
  
  // Four-momentum of the target.
  Target->SetP(mt, 0, 0, 0);
  Target->SetX(0, 0, 0, zr);
  Target->Print();
  
  // Total four momentum in the lab.
  FourVector Ptot = Beam->GetP() + Target->GetP();
  //   Ptot.SetName("Total four-mom. in the lab");
  
  // Exc. energy of compound particle.
  Compound->SetP(Ptot);
  Compound->Print();
  cout << "Eexc(" << Compound->Name << ") = " << sqrt(Ptot*Ptot) - mc << " MeV" << endl;

  // Randomly select the scattering angle in the center of mass, then
  // go to the laboratory reference frame.
  double theta_CM = acos(Rdm->Uniform(-1,1));
  double phi_CM = Rdm->Uniform(-pi,pi);
  cout << "theta_cm=" << theta_CM*180/pi << "   phi_cm=" << phi_CM*180/pi << endl;

  // Final momentum in the CM reference frame (Ptot*Ptot is an invariant quantity).
  double pf_CM = sqrt((Ptot*Ptot - pow(ml+mh,2))*(Ptot*Ptot - pow(ml-mh,2))/(4*(Ptot*Ptot)));
  //   cout << "Ex = " << Ex << " MeV    pf_CM = " << pf_CM << endl;
    
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
  Light->SetX(0, 0, 0, zr);
  Light->Print();
  
  // Four-momentum of the heavy recoil (lab).
  Heavy->SetP(Beam->GetP() + Target->GetP() - Light->GetP());
  Heavy->SetX(0, 0, 0, zr);
  Heavy->Print();

  return;
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
// Create the target particle object (currently the charge is not used).
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::SetTargetParticle(string Name)
{ 
  // Use the nuclide finder object to determine the mass and atomic
  // number of this particle. Mass must be in MeV/c^2) and Z (in e).
  double m = NuF->GetMass(Name, "MeV/c^2");
  int Z = NuF->GetZ(Name);
  Target = new Particle(Name, m, Z);
  return;
}


///////////////////////////////////////////////////////////////////////////////////
// Core methode of this class, where the simulation takes place. Logic:
// Loop over events
// 1. Set beam inital conditions (beam energy, position)
// 2. Randomly select a reaction point (based on the beam energy range in MUSIC)
// 3. Set initial conditions for heavy and light particles
// 4. Propagate heavy particle and calculate energy loss in the anode elements 
// 5. Propagate light particle and calculate energy loss in the anode elements
// 6. Compute detector response (i.e. DE for beam + light + heavy)
// 7. Display trace and particle trajecories
///////////////////////////////////////////////////////////////////////////////////
void MUSIC_Simulator::Simulate(int SegNum, int NEvents, double MaxTime, double UserDT)
{
  this->NEvents = NEvents;
  // Masses (MeV/c^2)
  double mb = Beam->Mass;
  double mt = Target->Mass;
  double mf = Compound->Mass;
  double ml = Light->Mass;
  double mh = Heavy->Mass;
  // Some kinematic variables
  double Kb_min, Kb_max, MinZ, MaxZ;
  double zb1, zb2, Kb1, Kb2;
  // 3D trajectories
  TrajH = new TEveStraightLineSet*[NEvents];
  TrajL = new TEveStraightLineSet*[NEvents];
  // Arrays where the energy loss in each strip will be saved
  //  double* DeltaEB_ave = CalcAverageBeamELoss(); // average beam energy loss
  double** DeltaEB = new double*[AnodeStps];   // beam
  double** DeltaEL = new double*[AnodeStps];   // light
  double** DeltaEH = new double*[AnodeStps];   // heavy
  for (int stp=0; stp<AnodeStps; stp++) {
    DeltaEB[stp] = new double[AnodeCols];
    DeltaEL[stp] = new double[AnodeCols];
    DeltaEH[stp] = new double[AnodeCols];
  }
  double DeltaE;
  // For randomizing the detector response
  TF1* Gaussian = 0;
  if (EneSigma!=0) 
    Gaussian = new TF1("Gaussian","gaus",0, 100);
  // Canvas for displaying the traces
  TCanvas* Can = new TCanvas("Can","Traces", 0, 0, 1918, 630);
  HELoss->Draw();
  Trace = new TGraph*[NEvents];
  
  cout << "Simulating MUSIC traces ... " << endl;
  
  // Get the beam energy limits in the selected strip (assuming the
  // beam direction is parallel to the z-axis).
  MinZ = 0;
  MaxZ = AnodeDZ[0][0];
  for (int n=1; n<=SegNum; n++) {
    MinZ += AnodeDZ[n-1][0];
    MaxZ += AnodeDZ[n][0];
  }
  cout << "\nMinZ = " << MinZ << "    MaxZ = "<< MaxZ << endl;
  Kb_max = Beam->GetFinalEnergy(0, Kb_after_window, MinZ, 0.001/*step size in cm*/);
  Kb_min = Beam->GetFinalEnergy(0, Kb_after_window, MaxZ, 0.001/*step size in cm*/);
  cout << "MaxKb = " << Kb_max << "    MinKb = "<< Kb_min << endl;

  //-------------------------------------------------------------------------------
  // Event for-loop
  for (int evt=0; evt<NEvents; evt++) {
    cout << "\n***************** Event " << evt << "\n" << endl;

    // Reset the detector response
    for (int stp=0; stp<AnodeStps; stp++) 
      for (int col=0; col<AnodeCols; col++) {
	DeltaEB[stp][col] = 0;
	DeltaEL[stp][col] = 0;
	DeltaEH[stp][col] = 0;
      }
    

    // 1. Set beam inital conditions (beam energy, position)
    SetInitialKinematics(Kb_after_window);   

    // Within the selected strip randomly select the energy at which the beam particle
    // interacts with the target and from this value calculate the reaction z.
    double Kbr = Rdm->Uniform(Kb_min, Kb_max);
    double zr = Beam->GetPathLength(0, Kb_after_window, Kbr, 0.01/*time step in ns*/);
    double TOF = Beam->GetTimeOfFlight(0, Kb_after_window, zr/*cm*/, zr/100/*cm*/);
    cout << "Kbr = " << Kbr << "  zr = " << zr << "  tof = " << TOF << endl;
    
    // Propagate the beam particle from the origin to the reaction point
    DeltaEB = PropagateParticle(Beam, evt, TOF, UserDT);

    // Set the kinematics of all particles at the reaction point
    SetReactionKinematics(Kbr, zr);

    // Propagate the particles and get the energy deposited in each
    // anode segment.
    DeltaEH = PropagateParticle(Heavy, evt, MaxTime, UserDT);
    DeltaEL = PropagateParticle(Light, evt, MaxTime, UserDT);

    // Clone the particle trajectories
    if (evt<Particle::MaxEvents) {
      TrajH[evt] = (TEveStraightLineSet*)Heavy->AllTraj[evt]->Clone();
      TrajL[evt] = (TEveStraightLineSet*)Light->AllTraj[evt]->Clone();
    }

    // Draw the energy loss traces (detector response)
    Trace[evt] = new TGraph();
    Trace[evt]->SetName(Form("evt %d",evt));
    for (int stp=0; stp<AnodeStps; stp++) {
      DeltaE = 0;
      for (int col=0; col<AnodeCols; col++) {	
	// Previously we were normalizing to the energy loss of the beam
	// minus the energy loss in the the dead layer, which in this
	// version of the code is typically strip 1.
	//  DeltaE = DeltaEB[n] + DeltaEL[n] + DeltaEH[n] - (DeltaEB_ave[n]-DeltaEB_ave[1]);
	// In this version, we normalize to the average energy loss of
	// the beam in each strip.
	DeltaE += DeltaEB[stp][col] + DeltaEL[stp][col] + DeltaEH[stp][col];// - DeltaEB_ave[n];
	cout << stp << " " << col << ": " << DeltaEB[stp][col] << " " << DeltaEL[stp][col] 
	     << " " << DeltaEH[stp][col] << endl;
      }
      //cout<< "Energy loss L: "<<DeltaEL[n]<<  "Energy loss H: "<<DeltaEH[n]<<endl;
      if (EneSigma!=0 && Gaussian!=0 && DeltaE>0) {
	Gaussian->SetRange(0.0, 2*DeltaE);
	Gaussian->SetParameters(1.0, DeltaE, EneSigma);
	DeltaE = Gaussian->GetRandom();
      }      
      Trace[evt]->SetPoint(stp, stp, DeltaE);
    }
    Trace[evt]->Draw("*l same");
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
  Output->Close();
  return;
}



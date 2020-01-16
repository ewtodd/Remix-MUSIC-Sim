// This file will be automatically read when starting a ROOT session.
{
  // Get rid of the -Wshadow flag for the compiler.
  TString CompilerString = gSystem->GetMakeSharedLib();
  CompilerString.ReplaceAll("-Wshadow", "");
  gSystem->SetMakeSharedLib(CompilerString);

  /////////////////////////////////////////////////////////////////////////////
  // Load the necessary libraries for the script to run.
  /////////////////////////////////////////////////////////////////////////////
  gStyle->SetOptStat("");  
  gSystem->Load("../../../physics-tools/EnergyLoss.so"); 
  gSystem->Load("../../../physics-tools/FourVector.so"); 
  gSystem->Load("../../../physics-tools/Particle.so"); 
  gSystem->Load("../../../physics-tools/NuclideFinder_cpp.so"); 
  // Special lib
  //  gSystem->Load("../../../physics-tools/SRIM_Table_Maker_cpp.so");
  gSystem->Load("../../MUSIC_Simulator_cpp.so");
}

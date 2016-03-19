{
  gSystem->Load("../../PhysicsTools/EnergyLoss.so"); 
  gSystem->Load("../../PhysicsTools/FourVector.so"); 
  gSystem->Load("../../PhysicsTools/Particle.so"); 
  gSystem->Load("../../PhysicsTools/SRIM_Table_Maker.so");
  gSystem->Load("../../NuclideFinder/NuclideFinder_cpp.so"); 
  // gROOT->ProcessLine(".L EnergyLoss.cpp+");
  // gROOT->ProcessLine(".L FourVector.cpp+");
  // gROOT->ProcessLine(".L Particle.cpp+");
  // gROOT->ProcessLine(".L SRIM_Table_Maker.cpp+");
  // gROOT->ProcessLine(".L NuclideFinder.cpp+");
  gROOT->ProcessLine(".L MUSIC_Simulator.cpp++");
}

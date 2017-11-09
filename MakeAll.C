{
  string IncludePath = "../physicstools/";
  int LibStatus = 0;
  LibStatus += gSystem->Load((IncludePath + "EnergyLoss.so").c_str());
  LibStatus += gSystem->Load((IncludePath + "FourVector.so").c_str());
  LibStatus += gSystem->Load((IncludePath + "Particle.so").c_str());
  LibStatus += gSystem->Load((IncludePath + "NuclideFinder_cpp.so").c_str());
  if (LibStatus==0) 
    gROOT->ProcessLine(".L MUSIC_Simulator.cpp++");
}

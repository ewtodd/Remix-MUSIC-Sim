{
  string IncludePath = "../physics-tools/";
  const int numlibs = 4;
  string mylib[numlibs];
  mylib[0] = IncludePath + "EnergyLoss.so";
  mylib[1] = IncludePath + "FourVector.so";
  mylib[2] = IncludePath + "Particle.so";
  mylib[3] = IncludePath + "NuclideFinder_cpp.so";
  
  int LibStatus = 0;
  for (int l=0; l<numlibs; l++)
    LibStatus += gSystem->Load(mylib[l].c_str());
  if (LibStatus==0) {
    for (int l=0; l<numlibs; l++)
      cout << mylib[l] << " loaded" << endl;
    gROOT->ProcessLine(".L MUSIC_Simulator.cpp++");
  }
  else 
    cout << "MUSIC Simulator not compiled! Physics Tools not loaded." << endl;
}

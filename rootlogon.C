// ROOT Logon file for ~/Dropbox/DataAnalisysProject/ANASEN/
// This file will be automatically read when starting a ROOT session.
{
  // Get rid of the -Wshadow flag for the compiler.
  TString CompilerString = gSystem->GetMakeSharedLib();
  CompilerString.ReplaceAll("-Wshadow", "");
  CompilerString.ReplaceAll("-Wunused-parameter", "");
  gSystem->SetMakeSharedLib(CompilerString);
  //  cout << CompilerString << endl;
}

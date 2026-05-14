///////////////////////////////////////////////////////////////////////////
// musicsim main.cpp                                                     //
// Copyright (C) 2021  MUSICSIM (pending)                                //
//                                                                       //
// This program is free software: you can redistribute it and/or modify  //
// it under the terms of the GNU General Public License as published by  //
// the Free Software Foundation, either version 3 of the License, or     //
// (at your option) any later version.                                   //
//                                                                       //
// This program is distributed in the hope that it will be useful,       // 
// but WITHOUT ANY WARRANTY; without even the implied warranty of        //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         //
// GNU General Public License for more details.                          //
//                                                                       //
// You should have received a copy of the GNU General Public License     //
// along with this program.  If not, see <http://www.gnu.org/licenses/>. //
//                                                                       //
///////////////////////////////////////////////////////////////////////////
//                                                                       //
//  __  __ _    _  _____ _____ _____                                     //
// |  \/  | |  | |/ ____|_   _/ ____|                                    //
// | \  / | |  | | (___   | || |                                         //
// | |\/| | |  | |\___ \  | || |                                         //
// | |  | | |__| |____) |_| || |____                                     //
// |_|  |_|\____/|_____/|_____\_____|                                    //
//   _____ _____ __  __                                                  //
//  / ____|_   _|  \/  |                                                 //
// | (___   | | | \  / |                                                 //
//  \___ \  | | | |\/| |                                                 //
//  ____) |_| |_| |  | |                                                 //
// |_____/|_____|_|  |_|                                                 //
//                                                                       //
//                                                                       //
// Created by Daniel Santiago-Gonzalez (dsg@anl.gov)                     //
// For documentation, software updates and license details, please visit //
// https://gitlab.phy.anl.gov/music/sim or see README.md file.           //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "MUSIC_Simulator.hpp"
#include <TApplication.h>

using namespace std;

int main(int argc, char* argv[])
{

#ifndef MUSICSIM_VERSION
#  define MUSICSIM_VERSION "dev"
#endif
  cout << "==========================================================================" << endl;
  cout << "|--- MUSIC simulator (musicsim) version " << MUSICSIM_VERSION << endl;
  cout << "| Usage: ./musicsim control.toml                                         |" << endl;
  cout << "| See README.md for installation and usage.                              |" << endl;
  cout << "| Fork of https://gitlab.phy.anl.gov/music/sim (D. Santiago-Gonzalez)     |" << endl;
  cout << "==========================================================================" << endl;

  
  if (argc>2) {
    cout << "musicsim error: only one argument is expected." << endl;
    return 0;
  }
  
  auto MS = new MUSIC_Simulator(); 
  if (argc==1) {
    cout << "musicsim warning: no control file specified. Using default parameters." << endl;
  }
  else if (argc==2) {
    cout << "Loading control file: " << argv[1] << endl;
    if (MS->loadCtrlFile(argv[1])==0)
      cout << "musicsim warning: invalid control file (check address)." << endl;
  }
  TApplication rootApp("musicsim", &argc, argv);
  if (MS->run()) {
    cout << "To quit musicsim you have 3 options:\n"
	 << " 1) In the Eve Main Window, click 'Browser', then click 'Quit ROOT'\n"
	 << " 2) In the Chart window, click 'File', then click 'Quit ROOT'\n"
	 << " 3) In this terminal, type CTRL+C" << endl;
    rootApp.Run(kTRUE);
    rootApp.HandleException(kSigSegmentationViolation);
  }
    
  return 1;
}

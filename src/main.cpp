///////////////////////////////////////////////////////////////////////////
// andes.cpp                                                             //
// Copyright (C) 2020  ANDES                                             //
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
//     /\                                                                //
//    /  \   _ __ __ _  ___  _ __  _ __   ___                            //
//   / /\ \ | '__/ _` |/ _ \| '_ \| '_ \ / _ \                           //
//  / ____ \| | | (_| | (_) | | | | | | |  __/                           //
// /_/    \_\_|  \__, |\___/|_| |_|_| |_|\___|                           //
//                __/ |                                                  //
//  _   _        |___/                                                   //
// | \ | |          | |                                                  //
// |  \| |_   _  ___| | ___  __ _ _ __                                   //
// | . ` | | | |/ __| |/ _ \/ _` | '__|                                  //
// | |\  | |_| | (__| |  __/ (_| | |                                     //
// |_|_\_|\__,_|\___|_|\___|\__,_|_|                                     //
// |  __ \      | |                                                      //
// | |  | | __ _| |_ __ _                                                //
// | |  | |/ _` | __/ _` |                                               //
// | |__| | (_| | || (_| |                                               //
// |_____/ \__,_|\__\__,_|              _   _                            //
// |  ____|          | |               | | (_)                           //
// | |__  __  ___ __ | | ___  _ __ __ _| |_ _  ___  _ __                 //
// |  __| \ \/ / '_ \| |/ _ \| '__/ _` | __| |/ _ \| '_ \                //
// | |____ >  <| |_) | | (_) | | | (_| | |_| | (_) | | | |               //
// |______/_/\_\ .__/|_|\___/|_|  \__,_|\__|_|\___/|_| |_|               //
//   _____     | |__ _                                                   //
//  / ____|    |_/ _| |                                                  //
// | (___   ___ | |_| |___      ____ _ _ __ ___                          //
//  \___ \ / _ \|  _| __\ \ /\ / / _` | '__/ _ \                         //
//  ____) | (_) | | | |_ \ V  V / (_| | | |  __/                         //
// |_____/ \___/|_|  \__| \_/\_/ \__,_|_|  \___|                         //
//           _   _ _____  ______  _____                                  //
//     /\   | \ | |  __ \|  ____|/ ____|                                 //
//    /  \  |  \| | |  | | |__  | (___                                   //
//   / /\ \ | . ` | |  | |  __|  \___ \                                  //
//  / ____ \| |\  | |__| | |____ ____) |                                 //
// /_/    \_\_| \_|_____/|______|_____/                                  //
//                                                                       //
//                                                                       //
// Created by Daniel Santiago-Gonzalez (dsg@anl.gov)                     //
// For documentation, software updates and license details, please visit //
// https://gitlab.phy.anl.gov/andes/ndviz or see README.md file.         //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "MUSIC_Simulator.hpp"
#include <TApplication.h>

using namespace std;

int main(int argc, char* argv[])
{

  cout << "==========================================================================" << endl;
  cout << "|--- MUSIC simulator ----------------------------------------------------|" << endl;
  cout << "| Created by Daniel Santiago-Gonzalez (Argonne National Laboratory)      |" << endl;
  cout << "| ver 3.0 (2020/1)                                                       |" << endl;
  cout << "| Usage: ./andes control.file                                            |" << endl;
  cout << "| See README.md file for basic installation and usage.                   |" << endl;
  cout << "| For documentation, software updates and license details, please visit: |" << endl;
  cout << "| https://gitlab.phy.anl.gov/music/sim                                   |" << endl;
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
  MS->run();
  cout << "To quit musicsim you have 3 options:\n"
       << " 1) In the Eve Main Window, click 'Browser', then click 'Quit ROOT'\n"
       << " 2) In the Chart window, click 'File', then click 'Quit ROOT'\n"
       << " 3) In this terminal, type CTRL+Q" << endl;
  rootApp.Run(kTRUE);
  rootApp.HandleException(kSigSegmentationViolation);
    
  return 1;
}

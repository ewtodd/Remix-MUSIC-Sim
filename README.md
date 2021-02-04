# MUSIC simulator

Created by Daniel Santiago-Gonzalez (Argonne National Laboratory)
ver 3.0 (2021/1)


## Linux requirements

1. linux based operating system (tested in Ubuntu 20, Debian 9 and CentOS 7)
2. ROOT 6 with EVE utilities for 3D visualization (see https://root.cern.ch/eve). Tested with ROOT 6.18/02 and 6.22/06 prepackaged binary files downloaded from [this website](https://root.cern/install/all_releases/).
3. g++ compiler (tested with g++ 4.8.4)
4. git (tested with version 1.9.1)
5. sudo apt install libopengl-dev

## Installation instructions LINUX (for MAC see below)

1. Get the source code. `cd` to your preferred working directory then type

`git clone https://gitlab.phy.anl.gov/music/sim`

2. Compile it by typing 

`make`

3. Test it
cd Examples/Fusion16C
root -l RunMusicSim_16C.C

If everything went well you should see two windows, one with a 2D plot of the MUSIC traces and
another one with a 3D view of the event.
* End of LINUX installation instructions *


## Installation instructions MACOS (ver2, needs update for ver3)

1. Get the codes. `cd` to your preferred working directory, then type

`git clone https://gitlab.phy.anl.gov/dasago/music-simulator.git`

`git clone https://gitlab.phy.anl.gov/andes/physics-tools.git`

2. Compile it (ignore all warnings but if there are errors send them to dasago@anl.gov)

`cd physics-tools`

`root -l`

`.L EnergyLoss.cpp++ `

`.L FourVector.cpp++`

`.L NuclideFinder.cpp++`

`.L Particle.cpp++`

`.L SRIM_Table_Maker.cpp++`

`.q`

`cd ../music-simulator`

Edit the names of all the .so library files in the MakeAll.C file. (Ex: EnergyLoss.so -> EnergyLoss_cpp.so)

Open MUSIC_Simulator.cpp 
Change line 834 from `...stp<AnodeStps+1;...` to `...stp<AnodeStps;...` (Delete the +1)

`root -l MakeAll.C`

`.q`


3. Test it

`cd Examples/Fusion16C`

Change the paths to the .so files, and edit the library file names from .so to _cpp.so

`root -l RunMusicSim_16C.C`

If everything went well you should see two windows, one with a 2D plot of the MUSIC traces and
another one with a 3D view of the event.
End of MACOS installation instructions.

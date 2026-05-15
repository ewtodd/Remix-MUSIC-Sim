#include <iostream>

#include <TApplication.h>

#include "Simulator.hpp"

#ifndef MUSICSIM_VERSION
#define MUSICSIM_VERSION "dev"
#endif

Int_t main(Int_t argc, char *argv[]) {
  std::cout
      << "==========================================================================\n"
      << "|--- MUSIC simulator (musicsim) version " << MUSICSIM_VERSION << "\n"
      << "| Usage: ./musicsim control.toml                                         |\n"
      << "| See README.md for installation and usage.                              |\n"
      << "| Fork of https://gitlab.phy.anl.gov/music/sim (D. Santiago-Gonzalez)    |\n"
      << "==========================================================================\n";

  if (argc > 2) {
    std::cout << "musicsim error: only one argument is expected." << std::endl;
    return 0;
  }

  auto *MS = new Simulator();
  if (argc == 1) {
    std::cout
        << "musicsim warning: no control file specified. Using default parameters."
        << std::endl;
  } else {
    std::cout << "Loading control file: " << argv[1] << std::endl;
    if (MS->loadCtrlFile(argv[1]) == 0)
      std::cout << "musicsim warning: invalid control file (check address)."
                << std::endl;
  }

  TApplication rootApp("musicsim", &argc, argv);
  if (MS->run()) {
    std::cout << "To quit musicsim:\n"
              << "  1) In the Eve Main Window, click 'Browser' -> 'Quit ROOT'\n"
              << "  2) In the Chart window, click 'File' -> 'Quit ROOT'\n"
              << "  3) Hit CTRL+C in this terminal" << std::endl;
    rootApp.Run(kTRUE);
    rootApp.HandleException(kSigSegmentationViolation);
  }
  return 1;
}

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include <TApplication.h>
#include <TString.h>

#include "Simulator.hpp"

#ifndef MUSICSIM_VERSION
#define MUSICSIM_VERSION "dev"
#endif
#ifndef MUSICSIM_ASSETS_DIR
#define MUSICSIM_ASSETS_DIR "assets"
#endif

// Wordmark, printed once before anything else (including the multi-threaded
// simulate path, which is only ever reached after main() has already run).
// Read from assets/remix.txt so the same art can be pasted into the README.
static void PrintLogo() {
  std::ifstream in(TString(MUSICSIM_ASSETS_DIR) + "/remix.txt");
  if (!in)
    return; // art is decoration; never let a missing asset stop a run
  std::cout << std::endl << in.rdbuf() << std::endl;
}

Int_t main(Int_t argc, char *argv[]) {
  PrintLogo();
  std::cout
      << "==========================================================================\n"
      << "|--- MUSIC simulator (musicsim) version " << MUSICSIM_VERSION << "\n"
      << "| Usage: ./musicsim [--check] control.toml                               |\n"
      << "| See README.md for installation and usage.                              |\n"
      << "| Fork of https://gitlab.phy.anl.gov/music/sim (D. Santiago-Gonzalez)    |\n"
      << "==========================================================================\n";

  const Bool_t checkOnly = argc == 3 && std::string(argv[1]) == "--check";
  if ((!checkOnly && argc != 2) ||
      (argc >= 2 && std::string(argv[1]) == "--check" && argc != 3)) {
    std::cerr << "musicsim error: usage: musicsim [--check] control.toml"
              << std::endl;
    return EXIT_FAILURE;
  }

  auto simulator = std::make_unique<Simulator>();
  const char *controlPath = argv[checkOnly ? 2 : 1];
  std::cout << "Loading control file: " << controlPath << std::endl;
  if (simulator->loadCtrlFile(controlPath) == 0)
    return EXIT_FAILURE;
  if (checkOnly) {
    std::cout << "Control file is valid." << std::endl;
    return EXIT_SUCCESS;
  }

  TApplication rootApp("musicsim", &argc, argv);
  const Bool_t visualize = simulator->WantsVisualization();
  if (simulator->run() == 0)
    return EXIT_FAILURE;
  if (visualize) {
    std::cout << "To quit musicsim:\n"
              << "  1) In the Eve Main Window, click 'Browser' -> 'Quit ROOT'\n"
              << "  2) In the Chart window, click 'File' -> 'Quit ROOT'\n"
              << "  3) Hit CTRL+C in this terminal" << std::endl;
    rootApp.Run(kTRUE);
    rootApp.HandleException(kSigSegmentationViolation);
  }
  return EXIT_SUCCESS;
}

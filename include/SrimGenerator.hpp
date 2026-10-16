#ifndef SRIMGENERATOR_HPP
#define SRIMGENERATOR_HPP

// Generates missing SRIM stopping tables with SRIM-nix's make-srim-table.
// srim-cache builds the beam's and the steps' tables; talys-cache builds those
// of a talys_class's products, which only its TALYS tables name. The flake
// compiles the generator's store path in, SRIM_TABLE_BIN in the environment
// overrides it, and with neither the name is looked up on PATH.
//
// Table naming and location must match EnergyLoss::BuildTables exactly:
//   <dir of run.output>/<ion>_in_<gas>_<P>Torr_<T>K.srim

#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <set>
#include <spawn.h>
#include <string>
#include <sys/stat.h>
#include <sys/wait.h>
#include <vector>

extern char **environ;

#ifndef MUSICSIM_SRIM_TABLE_BIN
#define MUSICSIM_SRIM_TABLE_BIN ""
#endif

namespace music {

inline std::string SrimTableGenerator() {
  const char *env = std::getenv("SRIM_TABLE_BIN");
  if (env && env[0] != '\0')
    return env;
  if (MUSICSIM_SRIM_TABLE_BIN[0] != '\0')
    return MUSICSIM_SRIM_TABLE_BIN;
  return "make-srim-table";
}

inline bool SrimFileExists(const std::string &path) {
  struct stat st;
  return ::stat(path.c_str(), &st) == 0;
}

// SRIM needs a real ion. Reaction steps also name things like "n" or "gamma",
// which deposit nothing and have no stopping table; skip anything that is not
// mass-number-then-symbol.
inline bool IsSrimIon(const std::string &name) {
  size_t i = 0;
  while (i < name.size() && std::isdigit(static_cast<unsigned char>(name[i])))
    i++;
  if (i == 0 || i == name.size())
    return false;
  for (size_t k = i; k < name.size(); k++)
    if (!std::isalpha(static_cast<unsigned char>(name[k])))
      return false;
  return true;
}

inline std::string SrimNumberArg(double value) {
  char text[64];
  std::snprintf(text, sizeof(text), "%g", value);
  return text;
}

inline bool RunSrimGenerator(const std::string &tool,
                             const std::string &generator,
                             const std::string &ion, const std::string &gas,
                             double pressure, double temperature,
                             const std::string &output) {
  std::vector<std::string> args = {generator,
                                   "--ion",
                                   ion,
                                   "--gas",
                                   gas,
                                   "--pressure",
                                   SrimNumberArg(pressure),
                                   "--temp",
                                   SrimNumberArg(temperature),
                                   "--out",
                                   output};
  std::vector<char *> argv;
  argv.reserve(args.size() + 1);
  for (size_t a = 0; a < args.size(); ++a)
    argv.push_back(args[a].data());
  argv.push_back(nullptr);

  pid_t child = -1;
  const int spawnError = ::posix_spawnp(&child, generator.c_str(), nullptr,
                                        nullptr, argv.data(), environ);
  if (spawnError != 0) {
    std::cerr << tool << ": cannot start generator '" << generator
              << "': " << std::strerror(spawnError) << std::endl;
    return false;
  }
  int status = 0;
  while (::waitpid(child, &status, 0) < 0) {
    if (errno == EINTR)
      continue;
    std::cerr << tool << ": waitpid failed: " << std::strerror(errno)
              << std::endl;
    return false;
  }
  return WIFEXITED(status) && WEXITSTATUS(status) == 0;
}

struct SrimCounts {
  int built = 0, kept = 0, failed = 0;
};

// Generates each ion's table in `dir` unless it is already there.
inline SrimCounts EnsureSrimTables(const std::string &tool,
                                   const std::set<std::string> &ions,
                                   const std::string &gas, double pressure,
                                   double temperature, const std::string &dir) {
  char tag[128];
  std::snprintf(tag, sizeof(tag), "%s_%gTorr_%gK", gas.c_str(), pressure,
                temperature);
  SrimCounts counts;
  for (std::set<std::string>::const_iterator ion = ions.begin();
       ion != ions.end(); ++ion) {
    const std::string out = dir + "/" + *ion + "_in_" + tag + ".srim";
    if (SrimFileExists(out)) {
      counts.kept++;
      continue;
    }
    std::cout << tool << ": generating " << out << std::endl;
    const std::string generator = SrimTableGenerator();
    if (!RunSrimGenerator(tool, generator, *ion, gas, pressure, temperature,
                          out) ||
        !SrimFileExists(out)) {
      std::cerr << tool << ": FAILED to generate " << out << std::endl
                << "  (generator: " << generator
                << "; set SRIM_TABLE_BIN or build through the flake)"
                << std::endl;
      counts.failed++;
      continue;
    }
    counts.built++;
  }
  return counts;
}

} // namespace music

#endif

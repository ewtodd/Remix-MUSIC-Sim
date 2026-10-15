// talys-cache <control.toml>: runs TALYS for each per-strip level table the
// control needs that is not cached yet; musicsim only reads the tables.
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

#include "Simulator.hpp"
#include "TalysTables.hpp"

#ifndef MUSICSIM_TALYS_BIN
#define MUSICSIM_TALYS_BIN ""
#endif

namespace {

// TALYS_BIN, else the path the build compiled in, else `talys` on PATH.
std::string TalysBinary() {
  const char *env = std::getenv("TALYS_BIN");
  if (env && env[0] != '\0')
    return env;
  if (MUSICSIM_TALYS_BIN[0] != '\0')
    return MUSICSIM_TALYS_BIN;
  return "talys";
}

// Runs TALYS in `work`, reading talys.inp and writing talys.out there.
Bool_t RunTalys(const std::string &talys, const std::string &work,
                std::string &error) {
  const pid_t child = fork();
  if (child < 0) {
    error = std::string("cannot fork: ") + std::strerror(errno);
    return kFALSE;
  }
  if (child == 0) {
    const int input = open((work + "/talys.inp").c_str(), O_RDONLY);
    const int output =
        open((work + "/talys.out").c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (input < 0 || output < 0 || chdir(work.c_str()) != 0 ||
        dup2(input, 0) < 0 || dup2(output, 1) < 0 || dup2(output, 2) < 0)
      _exit(126);
    execlp(talys.c_str(), talys.c_str(), static_cast<char *>(nullptr));
    _exit(127);
  }
  int status = 0;
  while (waitpid(child, &status, 0) < 0) {
    if (errno != EINTR) {
      error = std::string("waitpid failed: ") + std::strerror(errno);
      return kFALSE;
    }
  }
  if (WIFEXITED(status) && WEXITSTATUS(status) == 127) {
    error = "cannot start TALYS '" + talys +
            "' (set TALYS_BIN or build through the flake)";
    return kFALSE;
  }
  if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
    error = "TALYS '" + talys + "' failed";
    return kFALSE;
  }
  return kTRUE;
}

// The one file TALYS wrote in `work` named prefix*suffix, or empty if none.
Bool_t OutputFile(const std::string &work, const std::string &prefix,
                  const std::string &suffix, std::string &path,
                  std::string &error) {
  path.clear();
  std::error_code fileError;
  for (std::filesystem::directory_iterator entry(work, fileError);
       !fileError && entry != std::filesystem::directory_iterator();
       entry.increment(fileError)) {
    const std::string name = entry->path().filename().string();
    if (name.rfind(prefix, 0) != 0 || name.size() < suffix.size() ||
        name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0)
      continue;
    if (!path.empty()) {
      error = "more than one TALYS " + prefix + "*" + suffix + " file";
      return kFALSE;
    }
    path = entry->path().string();
  }
  if (fileError) {
    error = "cannot list '" + work + "': " + fileError.message();
    return kFALSE;
  }
  return kTRUE;
}

// TALYS's version line, e.g. "TALYS-2.2 (Version: December 18, 2025)".
std::string TalysVersion(const std::string &talysOut) {
  std::istringstream lines(talysOut);
  std::string line;
  while (std::getline(lines, line)) {
    const size_t start = line.find("TALYS-");
    if (start != std::string::npos)
      return line.substr(start);
  }
  return "TALYS (version not reported)";
}

// Runs TALYS for one request and writes its table atomically. Work names carry
// the PID, so concurrent runs cannot clash; the work dir is kept on failure.
Bool_t BuildTable(const std::string &talys, const music::TalysRequest &request,
                  std::string &error) {
  const std::filesystem::path table(request.path);
  const std::string pid = std::to_string(getpid());
  const std::string work = request.path + ".work." + pid;
  std::error_code fileError;
  std::filesystem::create_directories(table.parent_path(), fileError);
  std::filesystem::remove_all(work, fileError);
  if (!std::filesystem::create_directory(work, fileError)) {
    error = "cannot create '" + work + "'";
    return kFALSE;
  }
  {
    std::ofstream input(work + "/talys.inp");
    input << request.deck;
    if (!input) {
      error = "cannot write '" + work + "/talys.inp'";
      return kFALSE;
    }
  }
  if (!request.potential.empty()) {
    std::filesystem::create_symlink(request.potential,
                                    work + "/alphaomp9real.gnu", fileError);
    if (fileError) {
      error = "cannot link '" + request.potential + "': " + fileError.message();
      return kFALSE;
    }
  }
  if (!RunTalys(talys, work, error)) {
    error += "; see " + work;
    return kFALSE;
  }

  // The population file is absent when no binary channel is open; the (a,a)
  // elastic file is aa*ang.L00 (aX*ang.L00 holds other channels).
  std::string talysOut, populationPath, population, elasticPath, angular;
  if (!music::ReadFileText(work + "/talys.out", talysOut) ||
      !OutputFile(work, "populationE", ".out", populationPath, error) ||
      (!populationPath.empty() &&
       !music::ReadFileText(populationPath, population)) ||
      !OutputFile(work, "aa", "ang.L00", elasticPath, error) ||
      elasticPath.empty() || !music::ReadFileText(elasticPath, angular)) {
    if (error.empty())
      error = "cannot read TALYS output";
    error += "; see " + work;
    return kFALSE;
  }
  std::string rows, elasticRows;
  Double_t sigma = 0.0;
  Bool_t rutherfordOnly = kFALSE;
  // TALYS's projectile is the alpha, Z = 2.
  if (!music::TalysLevelRows(talysOut, population, request.Z, request.A, rows,
                             sigma, error) ||
      !music::TalysElasticRows(angular, request.Z, 2, request.ecm, elasticRows,
                               rutherfordOnly, error)) {
    error += "; see " + work;
    return kFALSE;
  }
  if (sigma == 0.0)
    std::cout << "talys-cache: strip " << request.strip
              << ": no (a,a') channel open at E_cm = " << request.ecm
              << " MeV, so its reactions will all be elastic" << std::endl;
  if (rutherfordOnly)
    std::cout << "talys-cache: strip " << request.strip
              << ": TALYS gave no elastic angular distribution at E_cm = "
              << request.ecm << " MeV; elastic will be pure Rutherford"
              << std::endl;

  char sigmaLine[96];
  std::snprintf(sigmaLine, sizeof(sigmaLine), "# sigma_inelastic_mb = %.9e",
                sigma);
  std::ostringstream text;
  text << "# talys-cache level table: residue excitation after binary alpha "
          "emission, weights in mb"
       << std::endl
       << "# " << TalysVersion(talysOut) << std::endl
       << sigmaLine << std::endl
       << (rutherfordOnly
               ? "# elastic: TALYS gave no angular distribution; "
                 "pure Rutherford"
               : "# elastic: TALYS shape + compound elastic over Rutherford "
                 "by theta_cm [deg]")
       << std::endl
       << music::TalysDeckComment(request.deck) << rows << elasticRows;
  const std::string staged = request.path + ".tmp." + pid;
  {
    std::ofstream output(staged, std::ios::binary);
    output << text.str();
    if (!output) {
      error = "cannot write '" + staged + "'";
      return kFALSE;
    }
  }
  std::filesystem::rename(staged, table, fileError);
  if (fileError) {
    error = "cannot publish '" + request.path + "': " + fileError.message();
    return kFALSE;
  }
  std::filesystem::remove_all(work, fileError);
  return kTRUE;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "usage: talys-cache <control.toml>" << std::endl;
    return 2;
  }
  Simulator simulator;
  if (simulator.loadCtrlFile(argv[1]) == 0)
    return 1;
  std::vector<music::TalysRequest> requests;
  if (!simulator.PlanTalysTables(requests))
    return 1;
  if (requests.empty()) {
    std::cout << "talys-cache: " << argv[1]
              << " has no talys branch without level_table; nothing to do"
              << std::endl;
    return 0;
  }

  const std::string talys = TalysBinary();
  Int_t built = 0, kept = 0, failed = 0, stopped = 0;
  for (size_t r = 0; r < requests.size(); ++r) {
    const music::TalysRequest &request = requests[r];
    if (request.stopped) {
      std::cout << "talys-cache: strip " << request.strip
                << ": the mean beam stops before the strip centre; no table, "
                   "and any reaction there will be elastic"
                << std::endl;
      ++stopped;
      continue;
    }
    std::error_code fileError;
    if (std::filesystem::exists(request.path, fileError)) {
      ++kept;
      continue;
    }
    std::cout << "talys-cache: strip " << request.strip
              << " (E_cm = " << request.ecm << " MeV): generating "
              << request.path << std::endl;
    std::string error;
    if (!BuildTable(talys, request, error)) {
      std::cerr << "talys-cache: FAILED " << request.path << ": " << error
                << std::endl;
      ++failed;
      continue;
    }
    ++built;
  }
  std::cout << "talys-cache: " << built << " generated, " << kept
            << " already cached, " << failed << " failed";
  if (stopped > 0)
    std::cout << ", " << stopped << " beyond the beam's range";
  std::cout << std::endl;
  return failed ? 1 : 0;
}

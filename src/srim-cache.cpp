// srim-cache — generate the SRIM stopping tables one control file needs.
//
// musicsim reads SRIM tables but never creates them: generating one drives the
// SRIM SR-Module under wine, which is slow and needs a toolchain the simulator
// has no business carrying. This tool does that step once, up front, so the
// simulation itself only ever reads a table that is already there.
//
//   srim-cache <control.toml>
//
// It is a no-op unless [physics].stopping is "srim" or "mean", and a no-op per
// table that already exists, so running it over a whole batch of control files
// is cheap:
// each distinct (ion, gas, pressure, temperature) is built once and then shared
// by every control file that asks for it.
//
// Table naming and location must match EnergyLoss::BuildTables exactly:
//   <dir of run.output>/<ion>_in_<gas>_<P>Torr_<T>K.srim
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <set>
#include <string>
#include <sys/stat.h>
#include <toml++/toml.hpp>
#include <vector>

namespace {

bool FileExists(const std::string &p) {
  struct stat st;
  return ::stat(p.c_str(), &st) == 0;
}

// SRIM needs a real ion. Reaction steps also name things like "n" or "gamma",
// which deposit nothing and have no stopping table; skip anything that is not
// mass-number-then-symbol.
bool IsIonName(const std::string &s) {
  size_t i = 0;
  while (i < s.size() && std::isdigit(static_cast<unsigned char>(s[i])))
    i++;
  if (i == 0 || i == s.size())
    return false;
  for (size_t k = i; k < s.size(); k++)
    if (!std::isalpha(static_cast<unsigned char>(s[k])))
      return false;
  return true;
}

std::string DirOf(const std::string &path) {
  size_t slash = path.find_last_of('/');
  return (slash == std::string::npos) ? std::string(".")
                                      : path.substr(0, slash);
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "usage: srim-cache <control.toml>" << std::endl;
    return 2;
  }
  const std::string ctrl = argv[1];

  toml::table tbl;
  try {
    tbl = toml::parse_file(ctrl);
  } catch (const toml::parse_error &err) {
    std::cerr << "srim-cache: cannot parse " << ctrl << ": "
              << err.description() << std::endl;
    return 1;
  }

  const std::string model =
      tbl["physics"]["stopping"].value_or(std::string("catima"));
  // "mean" averages catima with the SRIM tables, so it needs them cached too.
  if (model != "srim" && model != "mean") {
    std::cout << "srim-cache: [physics].stopping is neither \"srim\" nor "
                 "\"mean\" in "
              << ctrl << "; nothing to do" << std::endl;
    return 0;
  }

  const std::string gas = tbl["gas"]["species"].value_or(std::string("4He"));
  const double pressure = tbl["gas"]["pressure"].value_or(760.0);
  const double temperature = tbl["gas"]["temperature"].value_or(293.0);
  const std::string output = tbl["run"]["output"].value_or(std::string(""));
  if (output.empty()) {
    std::cerr << "srim-cache: " << ctrl << " has no [run].output, so there is "
              << "nowhere to put the tables" << std::endl;
    return 1;
  }

  // Every ion that will traverse the gas: the beam, the target nucleus it can
  // recoil off, and each reaction step's ejectile and residue.
  std::set<std::string> ions;
  auto add = [&ions](const std::string &n) {
    if (IsIonName(n))
      ions.insert(n);
  };
  add(tbl["beam"]["species"].value_or(std::string("")));
  add(tbl["target"]["species"].value_or(std::string("")));
  if (auto *steps = tbl["reaction"]["step"].as_array()) {
    for (auto &&s : *steps) {
      if (auto *st = s.as_table()) {
        add((*st)["evap"]["name"].value_or(std::string("")));
        add((*st)["res"]["name"].value_or(std::string("")));
      }
    }
  }
  if (ions.empty()) {
    std::cerr << "srim-cache: no ion names found in " << ctrl << std::endl;
    return 1;
  }

  char tag[128];
  snprintf(tag, sizeof(tag), "%s_%gTorr_%gK", gas.c_str(), pressure,
           temperature);
  const std::string dir = DirOf(output);

  int built = 0, kept = 0, failed = 0;
  for (const std::string &ion : ions) {
    const std::string out = dir + "/" + ion + "_in_" + tag + ".srim";
    if (FileExists(out)) {
      kept++;
      continue;
    }
    char cmd[1024];
    snprintf(cmd, sizeof(cmd),
             "make-srim-table --ion %s --gas %s --pressure %g --temp %g "
             "--out '%s'",
             ion.c_str(), gas.c_str(), pressure, temperature, out.c_str());
    std::cout << "srim-cache: generating " << out << std::endl;
    if (std::system(cmd) != 0 || !FileExists(out)) {
      std::cerr << "srim-cache: FAILED to generate " << out
                << "\n  (needs make-srim-table on PATH, from SRIM-nix)"
                << std::endl;
      failed++;
      continue;
    }
    built++;
  }
  std::cout << "srim-cache: " << built << " generated, " << kept
            << " already cached, " << failed << " failed" << std::endl;
  return failed ? 1 : 0;
}

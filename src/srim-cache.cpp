// srim-cache — generate the SRIM tables one control file needs.
//
// musicsim reads SRIM tables but never creates them: generating one drives the
// SRIM SR-Module under wine, which is slow and needs a toolchain the simulator
// has no business carrying. This tool does that step once, up front, so the
// simulation itself only ever reads a table that is already there. The
// generator is SRIM-nix's make-srim-table (see SrimGenerator.hpp).
//
//   srim-cache <control.toml>
//
// It is a no-op unless [physics].stopping is "srim" or "mean", and a no-op per
// table that already exists, so running it over a whole batch of control files
// is cheap:
// each distinct (ion, gas, pressure, temperature) is built once and then shared
// by every control file that asks for it. A talys_class's products are only
// known from its TALYS tables, so talys-cache builds theirs.
#include <iostream>
#include <set>
#include <string>
#include <toml++/toml.hpp>

#include "SrimGenerator.hpp"

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
  std::vector<std::string> names;
  names.push_back(tbl["beam"]["species"].value_or(std::string("")));
  names.push_back(tbl["target"]["species"].value_or(std::string("")));
  if (const toml::array *steps = tbl["reaction"]["step"].as_array()) {
    for (size_t s = 0; s < steps->size(); ++s) {
      if (const toml::table *step = (*steps)[s].as_table()) {
        names.push_back((*step)["evap"]["name"].value_or(std::string("")));
        names.push_back((*step)["res"]["name"].value_or(std::string("")));
      }
    }
  }
  std::set<std::string> ions;
  for (size_t n = 0; n < names.size(); ++n)
    if (music::IsSrimIon(names[n]))
      ions.insert(names[n]);
  if (ions.empty()) {
    std::cerr << "srim-cache: no ion names found in " << ctrl << std::endl;
    return 1;
  }

  const size_t slash = output.find_last_of('/');
  const std::string dir =
      slash == std::string::npos ? std::string(".") : output.substr(0, slash);
  const music::SrimCounts counts = music::EnsureSrimTables(
      "srim-cache", ions, gas, pressure, temperature, dir);
  std::cout << "srim-cache: " << counts.built << " generated, " << counts.kept
            << " already cached, " << counts.failed << " failed" << std::endl;
  return counts.failed ? 1 : 0;
}

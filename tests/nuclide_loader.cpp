#include "NuclideFinder.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace {

std::string FixedRecord() {
  std::string line(100, ' ');
  line.replace(6, 3, "  1");
  line.replace(11, 3, "  1");
  line.replace(16, 3, "  2");
  line.replace(20, 2, "H ");
  line.replace(29, 12, "       0.000");
  return line;
}

std::string NubaseRecord() {
  std::string line(100, ' ');
  line.replace(0, 3, "002");
  line.replace(4, 3, "001");
  line.replace(7, 1, "0");
  line.replace(11, 5, "2H   ");
  line.replace(79, 14, "1+            ");
  return line;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2)
    return 2;
  const std::filesystem::path directory = argv[1];
  const auto ame = directory / "valid.ame";
  const auto badAme = directory / "bad.ame";
  const auto nubase = directory / "valid.nubase";
  const auto badNubase = directory / "bad.nubase";
  {
    std::ofstream out(ame);
    out << "MASS LIST\nheader 1\nheader 2\nheader 3\nheader 4\n"
        << FixedRecord() << '\n';
  }
  {
    std::ofstream out(badAme);
    out << "MASS LIST\nheader 1\nheader 2\nheader 3\nheader 4\nshort\n";
  }
  {
    std::ofstream out(nubase);
    out << "NUBASE\n#-------------------------------------------------\ncolumns\n"
        << NubaseRecord() << '\n';
  }
  {
    std::ofstream out(badNubase);
    out << "NUBASE\n#-------------------------------------------------\ncolumns\nshort\n";
  }

  NuclideFinder finder;
  if (!finder.LoadAMEFile(ame.string()) || finder.GetArraySize() != 1 ||
      !finder.Contains("2H") || std::fabs(finder.GetMass(0, "u") - 2.0) > 1e-12)
    return 1;
  if (finder.LoadAMEFile(badAme.string()) != 0 || finder.GetArraySize() != 1 ||
      !finder.Contains("2H"))
    return 1;
  if (!finder.LoadNubaseFile(nubase.string()) ||
      std::fabs(finder.GetGSspin(0) - 1.0f) > 1e-6f ||
      finder.GetGSparity(0) != 1)
    return 1;
  if (finder.LoadNubaseFile(badNubase.string()) != 0 ||
      std::fabs(finder.GetGSspin(0) - 1.0f) > 1e-6f ||
      finder.GetGSparity(0) != 1)
    return 1;
  return 0;
}

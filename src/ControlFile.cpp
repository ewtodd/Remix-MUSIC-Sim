#include "EnergyLoss.hpp"
#include "Simulator.hpp"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <sstream>
#include <type_traits>
#include <unordered_set>

namespace {

Int_t ParseZeffName(const std::string &name) {
  if (name == "none")
    return catima::z_eff_type::none;
  if (name == "pierce_blann")
    return catima::z_eff_type::pierce_blann;
  if (name == "anthony_landorf")
    return catima::z_eff_type::anthony_landorf;
  if (name == "hubert")
    return catima::z_eff_type::hubert;
  if (name == "winger")
    return catima::z_eff_type::winger;
  if (name == "schiwietz")
    return catima::z_eff_type::schiwietz;
  if (name == "global")
    return catima::z_eff_type::global;
  if (name == "atima14")
    return catima::z_eff_type::atima14;
  return -1;
}

// atima14 and global need catima's optional GLOBAL library, which this build
// lacks; without it catima returns z_eff = 0 or a wrong charge, not an error.
Bool_t NeedsCatimaGlobal(Int_t code) {
  return code == catima::z_eff_type::atima14 ||
         code == catima::z_eff_type::global;
}

std::string CatimaGlobalError(const std::string &key,
                              const std::string &value) {
  return "'" + key + "' = '" + value +
         "' requires catima's GLOBAL library, which this build lacks; use "
         "'winger', which is what atima14 computes below 30 MeV/u";
}

Int_t ParseLowEnergyName(const std::string &name) {
  if (name == "srim_85")
    return catima::low_energy_types::srim_85;
  if (name == "srim_95")
    return catima::low_energy_types::srim_95;
  return -1;
}

void CheckKeys(const toml::table &table,
               std::initializer_list<const char *> allowed,
               const std::string &path, std::vector<std::string> &errors) {
  std::unordered_set<std::string> names;
  for (const char *name : allowed)
    names.emplace(name);
  for (const auto &[key, node] : table) {
    const std::string name(key.str());
    if (!names.count(name))
      errors.push_back("unknown key '" + path + name + "'");
  }
}

// Z of an element symbol ("Y"), from the mass table's names; -1 if unknown.
Int_t ElementZ(NuclideFinder &finder, const std::string &symbol) {
  for (Int_t i = 0; i < finder.GetArraySize(); ++i) {
    const std::string name = finder.GetName(i);
    const size_t letters = name.find_first_not_of("0123456789");
    if (letters != std::string::npos && name.substr(letters) == symbol)
      return finder.GetZ(i);
  }
  return -1;
}

std::string JoinErrors(const std::vector<std::string> &errors) {
  std::ostringstream out;
  for (const auto &error : errors)
    out << "\n  - " << error;
  return out.str();
}

Int_t ParseExcitationName(const std::string &name) {
  if (name == "forced")
    return 0;
  if (name == "ground")
    return 1;
  if (name == "uniform")
    return 2;
  if (name == "talys")
    return 3;
  return -1;
}

Int_t ParseAngularName(const std::string &name) {
  if (name == "isotropic")
    return 0;
  if (name == "rutherford")
    return 1;
  return -1;
}

// `level_table` is required exactly when the model is talys. A relative path
// is resolved against the directory holding the control file.
void ReadLevelTable(const toml::table &table, const std::string &path,
                    Int_t residueExc, const std::filesystem::path &controlDir,
                    std::string &levelTable, std::vector<std::string> &errors) {
  const Bool_t present = table.contains("level_table");
  const std::optional<std::string> value =
      table["level_table"].value<std::string>();
  if (present && !value) {
    errors.push_back("'" + path + ".level_table' must be a string");
    return;
  }
  if (present && residueExc != 3) {
    errors.push_back("'" + path +
                     ".level_table' requires residue_excitation='talys'");
    return;
  }
  if (!present) {
    if (residueExc == 3)
      errors.push_back("missing required key '" + path +
                       ".level_table' (for TALYS's channels per strip, use "
                       "reaction.talys_class)");
    return;
  }
  const std::filesystem::path given(*value);
  levelTable =
      given.is_absolute() ? given.string() : (controlDir / given).string();
}

// Reads one [[reaction.step.branch]] table. 'forced' is refused here: the
// minimum-excitation chain it relies on assumes one model per step.
void ParseBranch(const toml::table &table, const std::string &path,
                 Double_t stepThetaMin, const std::filesystem::path &controlDir,
                 Simulator::ExcitationBranch &branch,
                 std::vector<std::string> &errors) {
  CheckKeys(table,
            {"residue_excitation", "angular_distribution", "theta_cm_min_deg",
             "weight", "level_table"},
            path + ".", errors);

  const std::optional<std::string> excitation =
      table["residue_excitation"].value<std::string>();
  if (!excitation)
    errors.push_back("missing or invalid string '" + path +
                     ".residue_excitation'");
  else
    branch.residueExc = ParseExcitationName(*excitation);
  if (excitation && branch.residueExc < 1)
    errors.push_back("'" + path +
                     ".residue_excitation' must be 'ground', 'uniform', or "
                     "'talys'");

  const std::optional<std::string> angular =
      table["angular_distribution"].value<std::string>();
  if (!angular)
    errors.push_back("missing or invalid string '" + path +
                     ".angular_distribution'");
  else
    branch.angularDist = ParseAngularName(*angular);
  if (angular && branch.angularDist < 0)
    errors.push_back("'" + path +
                     ".angular_distribution' must be 'isotropic' or "
                     "'rutherford'");

  const std::optional<Double_t> weight = table["weight"].value<Double_t>();
  if (!weight)
    errors.push_back("missing or invalid number '" + path + ".weight'");
  else if (!std::isfinite(*weight) || *weight < 0.0)
    errors.push_back("'" + path + ".weight' must be finite and >= 0");
  else
    branch.weight = *weight;

  branch.thetaCmMinDeg = stepThetaMin;
  if (table.contains("theta_cm_min_deg")) {
    const std::optional<Double_t> theta =
        table["theta_cm_min_deg"].value<Double_t>();
    if (!theta)
      errors.push_back("'" + path + ".theta_cm_min_deg' must be a number");
    else if (!std::isfinite(*theta) || *theta <= 0.0 || *theta >= 180.0)
      errors.push_back("'" + path + ".theta_cm_min_deg' must be in (0, 180)");
    else
      branch.thetaCmMinDeg = *theta;
  }

  ReadLevelTable(table, path, branch.residueExc, controlDir, branch.levelTable,
                 errors);
}

// Reads a step's branch tables and normalizes their weights. Branch tables
// replace the step-level model keys, so those must then be absent.
void ParseStepBranches(const toml::table &step, const std::string &base,
                       Double_t stepThetaMin,
                       const std::filesystem::path &controlDir,
                       std::vector<Simulator::ExcitationBranch> &branches,
                       std::vector<std::string> &errors) {
  const toml::array *tables = step["branch"].as_array();
  if (!tables || !tables->is_array_of_tables()) {
    errors.push_back("'" + base + ".branch' must be an array of tables");
    return;
  }
  if (step.contains("residue_excitation") ||
      step.contains("angular_distribution"))
    errors.push_back("'" + base +
                     "' cannot set residue_excitation or "
                     "angular_distribution alongside [[reaction.step.branch]]");
  if (tables->size() > static_cast<size_t>(Simulator::kMaxBranchesPerStep)) {
    errors.push_back("'" + base + "' has more than " +
                     std::to_string(Simulator::kMaxBranchesPerStep) +
                     " branches");
    return;
  }

  Double_t totalWeight = 0.0;
  for (size_t i = 0; i < tables->size(); ++i) {
    const std::string path = base + ".branch[" + std::to_string(i) + "]";
    Simulator::ExcitationBranch branch;
    ParseBranch(*(*tables)[i].as_table(), path, stepThetaMin, controlDir,
                branch, errors);
    if (branch.weight > 0.0)
      totalWeight += branch.weight;
    branches.push_back(branch);
  }
  if (!(totalWeight > 0.0)) {
    errors.push_back("'" + base + ".branch' needs at least one weight > 0");
    return;
  }
  for (size_t i = 0; i < branches.size(); ++i)
    if (branches[i].weight > 0.0)
      branches[i].weight /= totalWeight;
}

// Reads the optional [talys] block: the alpha optical model, extra TALYS
// keywords for talys-cache tables, and the class channels' thresholds.
void ParseTalysBlock(const toml::table &talys, music::TalysSettings &settings,
                     std::vector<std::string> &errors) {
  CheckKeys(talys,
            {"alpha_omp", "keywords", "min_channel_fraction", "min_channel_mb"},
            "talys.", errors);
  if (talys.contains("alpha_omp")) {
    const std::optional<int64_t> omp = talys["alpha_omp"].value<int64_t>();
    if (!omp || *omp < 1 || *omp > 9)
      errors.push_back("'talys.alpha_omp' must be an integer in [1, 9]");
    else
      settings.alphaOmp = static_cast<Int_t>(*omp);
  }
  if (talys.contains("min_channel_fraction")) {
    const std::optional<Double_t> fraction =
        talys["min_channel_fraction"].value<Double_t>();
    if (!fraction || !std::isfinite(*fraction) || *fraction < 0.0 ||
        *fraction >= 1.0)
      errors.push_back("'talys.min_channel_fraction' must be a number in "
                       "[0, 1)");
    else
      settings.minChannelFraction = *fraction;
  }
  if (talys.contains("min_channel_mb")) {
    const std::optional<Double_t> mb =
        talys["min_channel_mb"].value<Double_t>();
    if (!mb || !std::isfinite(*mb) || *mb < 0.0)
      errors.push_back("'talys.min_channel_mb' must be a finite number >= 0");
    else
      settings.minChannelMb = *mb;
  }
  if (!talys.contains("keywords"))
    return;
  const toml::array *keywords = talys["keywords"].as_array();
  if (!keywords) {
    errors.push_back("'talys.keywords' must be an array of strings");
    return;
  }
  for (size_t k = 0; k < keywords->size(); ++k) {
    const std::string where = "'talys.keywords[" + std::to_string(k) + "]'";
    const std::optional<std::string> keyword =
        (*keywords)[k].value<std::string>();
    if (!keyword) {
      errors.push_back(where + " must be a string");
      continue;
    }
    const std::string problem = music::TalysKeywordProblem(*keyword);
    if (!problem.empty())
      errors.push_back(where + " " + problem);
    else
      settings.keywords.push_back(*keyword);
  }
}

} // namespace

Int_t Simulator::loadCtrlFile(const char *fileName) {
  if (!fileName || fileName[0] == '\0') {
    std::cerr << "musicsim ERROR: a control-file path is required."
              << std::endl;
    return 0;
  }
  std::string text;
  if (!music::ReadFileText(fileName, text)) {
    std::cerr << "musicsim ERROR: cannot read control file '" << fileName
              << "'." << std::endl;
    return 0;
  }
  return loadCtrlText(text, fileName);
}

Int_t Simulator::loadCtrlText(
    const std::string &text, const char *sourcePath,
    const std::map<std::string, std::string> *levelTables) {
  toml::table tbl;
  try {
    tbl = toml::parse(text, std::string_view(sourcePath));
  } catch (const toml::parse_error &error) {
    std::cerr << "musicsim ERROR: failed to parse TOML control file '"
              << sourcePath << "': " << error.description() << " at "
              << error.source().begin << std::endl;
    return 0;
  }

  const std::filesystem::path controlDir =
      std::filesystem::path(sourcePath).parent_path();
  controlFileParams next;
  music::PhysicsConfig nextPhysics;
  nextPhysics.mean.z_effective = catima::z_eff_type::winger;
  nextPhysics.mean.low_energy = catima::low_energy_types::srim_95;
  nextPhysics.straggling = nextPhysics.mean;
  nextPhysics.straggling.z_effective = catima::z_eff_type::pierce_blann;
  std::vector<std::string> errors;

  CheckKeys(tbl,
            {"gas", "beam", "target", "windows", "detector", "reaction",
             "physics", "run", "talys"},
            "", errors);
  if (const toml::table *talys = tbl["talys"].as_table()) {
    next.hasTalysBlock = kTRUE;
    ParseTalysBlock(*talys, next.talys, errors);
  } else if (tbl["talys"]) {
    errors.push_back("'talys' must be a table");
  }

  auto nodeAt = [&](const std::string &path) { return tbl.at_path(path); };
  auto readString = [&](const std::string &path, TString &dest,
                        Bool_t required = kFALSE) {
    const auto node = nodeAt(path);
    if (!node) {
      if (required)
        errors.push_back("missing required key '" + path + "'");
      return;
    }
    if (auto value = node.value<std::string>())
      dest = *value;
    else
      errors.push_back("'" + path + "' must be a string");
  };
  auto readDouble = [&](const std::string &path, auto &dest,
                        Bool_t required = kFALSE) {
    const auto node = nodeAt(path);
    if (!node) {
      if (required)
        errors.push_back("missing required key '" + path + "'");
      return;
    }
    if (auto value = node.value<Double_t>()) {
      using Value = std::remove_reference_t<decltype(dest)>;
      dest = static_cast<Value>(*value);
    } else {
      errors.push_back("'" + path + "' must be a number");
    }
  };
  auto readInt = [&](const std::string &path, Int_t &dest,
                     Bool_t required = kFALSE) {
    const auto node = nodeAt(path);
    if (!node) {
      if (required)
        errors.push_back("missing required key '" + path + "'");
      return;
    }
    if (auto value = node.value<int64_t>()) {
      if (*value < std::numeric_limits<Int_t>::min() ||
          *value > std::numeric_limits<Int_t>::max())
        errors.push_back("'" + path +
                         "' is outside the supported integer range");
      else
        dest = static_cast<Int_t>(*value);
    } else {
      errors.push_back("'" + path + "' must be an integer");
    }
  };
  auto readBool = [&](const std::string &path, Bool_t &dest) {
    const auto node = nodeAt(path);
    if (!node)
      return;
    if (auto value = node.value<bool>())
      dest = *value;
    else
      errors.push_back("'" + path + "' must be true or false");
  };

  auto requireTable = [&](const char *name) -> const toml::table * {
    const auto node = tbl[name];
    if (!node) {
      errors.push_back(std::string("missing required table '") + name + "'");
      return nullptr;
    }
    const auto *table = node.as_table();
    if (!table)
      errors.push_back(std::string("'") + name + "' must be a table");
    return table;
  };

  if (const auto *gas = requireTable("gas"))
    CheckKeys(*gas, {"species", "pressure", "temperature"}, "gas.", errors);
  readString("gas.species", next.gas, kTRUE);
  readDouble("gas.pressure", next.pressure, kTRUE);
  readDouble("gas.temperature", next.temperature, kTRUE);

  if (const auto *beam = requireTable("beam"))
    CheckKeys(*beam, {"species", "energy", "energy_fwhm", "dedx_scale"},
              "beam.", errors);
  readString("beam.species", next.beamName, kTRUE);
  readDouble("beam.energy", next.BeamEnergy, kTRUE);
  readDouble("beam.energy_fwhm", next.KbFWHM);
  readDouble("beam.dedx_scale", next.dEdxScaleBeam);

  if (const auto *target = requireTable("target"))
    CheckKeys(*target, {"species", "compound"}, "target.", errors);
  readString("target.species", next.target, kTRUE);
  readString("target.compound", next.compound, kTRUE);

  if (const auto *windows = tbl["windows"].as_table()) {
    CheckKeys(*windows, {"entrance", "exit", "degrader"}, "windows.", errors);
    auto loadLayer = [&](const char *name, TString &material,
                         Double_t &thickness, Bool_t &byLength,
                         Double_t &scale) {
      const std::string base = std::string("windows.") + name;
      const auto node = tbl.at_path(base);
      if (!node)
        return;
      const auto *layer = node.as_table();
      if (!layer) {
        errors.push_back("'" + base + "' must be a table");
        return;
      }
      CheckKeys(*layer,
                {"material", "thickness_mg_cm2", "thickness_um", "dedx_scale"},
                base + ".", errors);
      const Bool_t hasAreal = bool(tbl.at_path(base + ".thickness_mg_cm2"));
      const Bool_t hasLinear = bool(tbl.at_path(base + ".thickness_um"));
      if (hasAreal == hasLinear) {
        errors.push_back("'" + base +
                         "' must specify exactly one of thickness_mg_cm2 or "
                         "thickness_um");
        return;
      }
      if (hasAreal) {
        readDouble(base + ".thickness_mg_cm2", thickness);
        byLength = kFALSE;
      } else {
        readDouble(base + ".thickness_um", thickness);
        byLength = kTRUE;
      }
      readString(base + ".material", material, thickness > 0.0);
      readDouble(base + ".dedx_scale", scale);
      if (!std::isfinite(thickness))
        errors.push_back("'" + base + "' thickness must be finite");
      if (!std::isfinite(scale) || scale <= 0.0)
        errors.push_back("'" + base + ".dedx_scale' must be finite and > 0");
    };
    loadLayer("entrance", next.entranceMaterial, next.entranceThickness,
              next.entranceByLength, next.entranceScale);
    loadLayer("exit", next.exitMaterial, next.exitThickness, next.exitByLength,
              next.exitScale);
    loadLayer("degrader", next.degraderMaterial, next.degraderLength,
              next.degraderByLength, next.degraderScale);
  } else if (tbl["windows"] && !tbl["windows"].is_table()) {
    errors.push_back("'windows' must be a table");
  }

  if (const auto *detector = requireTable("detector"))
    CheckKeys(*detector,
              {"eloss_bins", "max_eloss", "strip", "strip_first", "strip_last",
               "noise_sigma_mev", "resolution_fwhm_percent"},
              "detector.", errors);
  readInt("detector.eloss_bins", next.ELossBins, kTRUE);
  readDouble("detector.max_eloss", next.MaxELoss, kTRUE);
  readInt("detector.strip", next.strip);
  readInt("detector.strip_first", next.stripFirst);
  readInt("detector.strip_last", next.stripLast);

  auto loadNoise = [&](const char *key, Int_t mode) {
    const std::string path = std::string("detector.") + key;
    const auto node = tbl.at_path(path);
    if (!node)
      return;
    if (next.noiseMode != 0) {
      errors.push_back(
          "detector.noise_sigma_mev and "
          "detector.resolution_fwhm_percent are mutually exclusive");
      return;
    }
    next.noiseMode = mode;
    auto validNoise = [&](Double_t value, const std::string &where) {
      if (!std::isfinite(value) || (value < 0.0 && value != -1.0))
        errors.push_back("'" + where +
                         "' must be finite and >= 0, or -1 to disable");
    };
    if (auto value = node.value<Double_t>()) {
      validNoise(*value, path);
      next.Noise.fill(*value);
      return;
    }
    const auto *channels = node.as_table();
    if (!channels) {
      errors.push_back("'" + path + "' must be a number or channel table");
      return;
    }
    for (const auto &[rawKey, valueNode] : *channels) {
      const std::string channel(rawKey.str());
      const auto value = valueNode.value<Double_t>();
      if (!value) {
        errors.push_back("'" + path + "." + channel + "' must be a number");
        continue;
      }
      validNoise(*value, path + "." + channel);
      if (channel == "Cathode") {
        next.NoiseCathode = *value;
        continue;
      }
      Int_t strip = -1;
      Int_t column = -1;
      if (channel == "S0") {
        strip = 0;
        column = 0;
      } else if (channel == "S17") {
        strip = 17;
        column = 0;
      } else if (channel.size() >= 2 &&
                 (channel.front() == 'L' || channel.front() == 'R')) {
        try {
          size_t consumed = 0;
          const Int_t parsed = std::stoi(channel.substr(1), &consumed);
          if (consumed == channel.size() - 1 && parsed >= 1 && parsed <= 16) {
            strip = parsed;
            column = channel.front() == 'R' ? 0 : 1;
          }
        } catch (const std::exception &) {
        }
      }
      if (strip < 0)
        errors.push_back("invalid detector channel '" + channel + "' in '" +
                         path + "'");
      else
        next.Noise[ElectrodeIndex(strip, column)] = *value;
    }
  };
  loadNoise("noise_sigma_mev", 1);
  loadNoise("resolution_fwhm_percent", 2);

  if (const auto *physics = tbl["physics"].as_table()) {
    CheckKeys(*physics,
              {"z_effective", "low_energy", "straggling_z_effective",
               "straggling", "stopping"},
              "physics.", errors);
    if (auto value = tbl.at_path("physics.z_effective").value<std::string>()) {
      const Int_t code = ParseZeffName(*value);
      if (code < 0)
        errors.push_back("'physics.z_effective' has unknown value '" + *value +
                         "'");
      else if (NeedsCatimaGlobal(code))
        errors.push_back(CatimaGlobalError("physics.z_effective", *value));
      else
        nextPhysics.mean.z_effective = static_cast<UChar_t>(code);
    } else if (tbl.at_path("physics.z_effective")) {
      errors.push_back("'physics.z_effective' must be a string");
    }
    if (auto value = tbl.at_path("physics.low_energy").value<std::string>()) {
      const Int_t code = ParseLowEnergyName(*value);
      if (code < 0)
        errors.push_back("'physics.low_energy' must be 'srim_85' or 'srim_95'");
      else
        nextPhysics.mean.low_energy = static_cast<UChar_t>(code);
    } else if (tbl.at_path("physics.low_energy")) {
      errors.push_back("'physics.low_energy' must be a string");
    }
    nextPhysics.straggling = nextPhysics.mean;
    nextPhysics.straggling.z_effective = catima::z_eff_type::pierce_blann;
    if (auto value = tbl.at_path("physics.straggling_z_effective")
                         .value<std::string>()) {
      const Int_t code = ParseZeffName(*value);
      if (code < 0)
        errors.push_back(
            "'physics.straggling_z_effective' has unknown value '" + *value +
            "'");
      else if (NeedsCatimaGlobal(code))
        errors.push_back(
            CatimaGlobalError("physics.straggling_z_effective", *value));
      else
        nextPhysics.straggling.z_effective = static_cast<UChar_t>(code);
    } else if (tbl.at_path("physics.straggling_z_effective")) {
      errors.push_back("'physics.straggling_z_effective' must be a string");
    }
    readBool("physics.straggling", nextPhysics.stragglingEnabled);
    if (auto value = tbl.at_path("physics.stopping").value<std::string>()) {
      if (*value == "catima")
        next.stoppingModel = 0;
      else if (*value == "srim")
        next.stoppingModel = 1;
      else if (*value == "mean")
        next.stoppingModel = 2;
      else
        errors.push_back(
            "'physics.stopping' must be 'catima', 'srim', or 'mean'");
    } else if (tbl.at_path("physics.stopping")) {
      errors.push_back("'physics.stopping' must be a string");
    }
  } else if (tbl["physics"] && !tbl["physics"].is_table()) {
    errors.push_back("'physics' must be a table");
  }

  if (const auto *reaction = tbl["reaction"].as_table()) {
    CheckKeys(*reaction,
              {"step", "talys_class", "theta_cm_min_deg", "residue_dedx_scale",
               "evap_dedx_scale"},
              "reaction.", errors);
    const Bool_t classGiven = reaction->contains("talys_class");
    if (const std::optional<std::string> name =
            (*reaction)["talys_class"].value<std::string>())
      next.talysClass = *name;
    else if (classGiven)
      errors.push_back("'reaction.talys_class' must be a string");
    const char *classKeys[3] = {"theta_cm_min_deg", "residue_dedx_scale",
                                "evap_dedx_scale"};
    Double_t classValues[3] = {next.classThetaCmMinDeg, 1.0, 1.0};
    for (Int_t k = 0; k < 3; ++k) {
      if (!reaction->contains(classKeys[k]))
        continue;
      const std::string where = std::string("'reaction.") + classKeys[k] + "'";
      const std::optional<Double_t> value =
          (*reaction)[classKeys[k]].value<Double_t>();
      if (!classGiven)
        errors.push_back(where + " is only used with reaction.talys_class");
      else if (!value || !std::isfinite(*value) || *value <= 0.0)
        errors.push_back(where + " must be a finite number > 0");
      else
        classValues[k] = *value;
    }
    if (classValues[0] >= 180.0)
      errors.push_back("'reaction.theta_cm_min_deg' must be in (0, 180)");
    next.classThetaCmMinDeg = classValues[0];
    next.classResidueDedxScale = static_cast<Float_t>(classValues[1]);
    next.classEvapDedxScale = static_cast<Float_t>(classValues[2]);
    const auto stepsNode = tbl.at_path("reaction.step");
    const auto *steps = stepsNode.as_array();
    if (classGiven) {
      if (stepsNode)
        errors.push_back("'reaction.step' cannot be combined with "
                         "reaction.talys_class, whose chains come from TALYS");
    } else if (!steps) {
      errors.push_back("'reaction.step' must be an array of tables");
    } else if (steps->size() > controlFileParams::MaxNumEvapPart) {
      errors.push_back("reaction has more than 10 steps");
    } else {
      Int_t index = 0;
      for (const auto &element : *steps) {
        const auto *step = element.as_table();
        const std::string base = "reaction.step[" + std::to_string(index) + "]";
        if (!step) {
          errors.push_back(base + " must be a table");
          ++index;
          continue;
        }
        CheckKeys(*step,
                  {"evap", "res", "residue_excitation", "angular_distribution",
                   "theta_cm_min_deg", "level_table", "branch", "weight"},
                  base + ".", errors);
        if (step->contains("weight"))
          errors.push_back("'" + base +
                           ".weight' is only valid inside "
                           "[[reaction.step.branch]]");
        auto loadProduct = [&](const char *side, TString &name, Int_t &color,
                               Float_t &scale) {
          const auto *product = (*step)[side].as_table();
          const std::string productPath = base + "." + side;
          if (!product) {
            errors.push_back("missing required table '" + productPath + "'");
            return;
          }
          CheckKeys(*product, {"name", "color", "dedx_scale"},
                    productPath + ".", errors);
          if (auto value = (*product)["name"].value<std::string>())
            name = *value;
          else
            errors.push_back("missing or invalid string '" + productPath +
                             ".name'");
          if (auto value = (*product)["color"].value<int64_t>())
            color = static_cast<Int_t>(*value);
          else if ((*product)["color"])
            errors.push_back("'" + productPath + ".color' must be an integer");
          if (auto value = (*product)["dedx_scale"].value<Double_t>())
            scale = static_cast<Float_t>(*value);
          else if ((*product)["dedx_scale"])
            errors.push_back("'" + productPath +
                             ".dedx_scale' must be a number");
          if (!std::isfinite(scale) || scale <= 0.0f)
            errors.push_back("'" + productPath +
                             ".dedx_scale' must be finite and > 0");
        };
        loadProduct("evap", next.evap[index], next.colorEvap[index],
                    next.dEdxScaleEvap[index]);
        loadProduct("res", next.res[index], next.colorRes[index],
                    next.dEdxScaleRes[index]);

        if (auto value = (*step)["residue_excitation"].value<std::string>()) {
          next.residueExc[index] = ParseExcitationName(*value);
          if (next.residueExc[index] < 0)
            errors.push_back("'" + base +
                             ".residue_excitation' must be 'forced', "
                             "'ground', 'uniform', or 'talys'");
        } else if ((*step)["residue_excitation"])
          errors.push_back("'" + base +
                           ".residue_excitation' must be a string");

        if (auto value = (*step)["angular_distribution"].value<std::string>()) {
          next.angularDist[index] = ParseAngularName(*value);
          if (next.angularDist[index] < 0)
            errors.push_back(
                "'" + base +
                ".angular_distribution' must be 'isotropic' or 'rutherford'");
        } else if ((*step)["angular_distribution"])
          errors.push_back("'" + base +
                           ".angular_distribution' must be a string");
        if (auto value = (*step)["theta_cm_min_deg"].value<Double_t>())
          next.thetaCmMinDeg[index] = *value;
        else if ((*step)["theta_cm_min_deg"])
          errors.push_back("'" + base + ".theta_cm_min_deg' must be a number");
        if (!std::isfinite(next.thetaCmMinDeg[index]) ||
            next.thetaCmMinDeg[index] <= 0.0 ||
            next.thetaCmMinDeg[index] >= 180.0)
          errors.push_back("'" + base +
                           ".theta_cm_min_deg' must be in (0, 180)");
        ReadLevelTable(*step, base, next.residueExc[index], controlDir,
                       next.levelTable[index], errors);
        if (step->contains("branch")) {
          next.hasBranchTables[index] = kTRUE;
          ParseStepBranches(*step, base, next.thetaCmMinDeg[index], controlDir,
                            next.branches[index], errors);
        }
        ++index;
      }
      next.NumEvapPart = index;
    }
  } else if (tbl["reaction"] && !tbl["reaction"].is_table()) {
    errors.push_back("'reaction' must be a table");
  }

  if (const auto *run = requireTable("run"))
    CheckKeys(*run,
              {"events_per_strip", "threads", "wait", "update", "max_time",
               "sim_step", "method", "output", "file_opt", "print_opt", "seed"},
              "run.", errors);
  readInt("run.events_per_strip", next.NEvents, kTRUE);
  readInt("run.threads", next.Threads);
  readInt("run.wait", next.Wait);
  readInt("run.update", next.Update);
  readDouble("run.max_time", next.MaxTime);
  readDouble("run.sim_step", next.SimStep);
  readInt("run.method", next.Method);
  readString("run.output", next.FileName, kTRUE);
  readString("run.file_opt", next.FileOpt);
  readInt("run.print_opt", next.PrintOpt);
  if (const auto seedNode = tbl.at_path("run.seed")) {
    if (auto value = seedNode.value<int64_t>()) {
      if (*value <= 0 || uint64_t(*value) > std::numeric_limits<UInt_t>::max())
        errors.push_back("'run.seed' must be in [1, 4294967295]");
      else
        next.Seed = static_cast<ULong64_t>(*value);
    } else {
      errors.push_back("'run.seed' must be an integer");
    }
  } else {
    errors.push_back("missing required key 'run.seed'");
  }

  const Int_t unset = controlFileParams::kStripUnset;
  const Bool_t singleSet = next.strip != unset;
  const Bool_t firstSet = next.stripFirst != unset;
  const Bool_t lastSet = next.stripLast != unset;
  if (firstSet != lastSet)
    errors.push_back(
        "detector.strip_first and detector.strip_last must be set together");
  if (singleSet && firstSet)
    errors.push_back("detector.strip cannot be combined with a strip range");
  if (!singleSet && !firstSet)
    errors.push_back(
        "select detector.strip or detector.strip_first/strip_last");
  if (singleSet) {
    if (next.strip < -1 || next.strip > 17)
      errors.push_back("detector.strip must be -1 or in [0, 17]");
    next.stripFirst = next.stripLast = next.strip;
  } else if (firstSet) {
    if (next.stripFirst < 0 || next.stripLast > 17 ||
        next.stripFirst > next.stripLast)
      errors.push_back(
          "detector strip range must satisfy 0 <= first <= last <= 17");
  }

  auto finitePositive = [&](Double_t value, const char *name) {
    if (!std::isfinite(value) || value <= 0.0)
      errors.push_back(std::string("'") + name + "' must be finite and > 0");
  };
  auto finiteNonnegative = [&](Double_t value, const char *name) {
    if (!std::isfinite(value) || value < 0.0)
      errors.push_back(std::string("'") + name + "' must be finite and >= 0");
  };
  finiteNonnegative(next.pressure, "gas.pressure");
  finitePositive(next.temperature, "gas.temperature");
  finitePositive(next.BeamEnergy, "beam.energy");
  finiteNonnegative(next.KbFWHM, "beam.energy_fwhm");
  finitePositive(next.dEdxScaleBeam, "beam.dedx_scale");
  finitePositive(next.MaxELoss, "detector.max_eloss");
  finitePositive(next.MaxTime, "run.max_time");
  finitePositive(next.SimStep, "run.sim_step");
  if (next.ELossBins <= 0)
    errors.push_back("'detector.eloss_bins' must be > 0");
  if (next.NEvents <= 0)
    errors.push_back("'run.events_per_strip' must be > 0");
  if (next.FileName.Length() == 0)
    errors.push_back("'run.output' must not be empty");
  if (next.Threads <= 0 || next.Threads > 256)
    errors.push_back("'run.threads' must be in [1, 256]");
  if ((next.Wait != 0 && next.Wait != 1) ||
      (next.Update != 0 && next.Update != 1))
    errors.push_back("'run.wait' and 'run.update' must be 0 or 1");
  if (next.PrintOpt < 0 || next.PrintOpt > 2)
    errors.push_back("'run.print_opt' must be in [0, 2]");
  if (next.Method != 0)
    errors.push_back(
        "'run.method' must be 0; legacy trace-database mode is unsupported");
  next.FileOpt.ToLower();
  if (next.FileOpt != "recreate" && next.FileOpt != "update")
    errors.push_back("'run.file_opt' must be 'recreate' or 'update'");
  if (next.Threads > 1 && next.FileOpt == "update")
    errors.push_back(
        "run.file_opt='update' is not supported with run.threads > 1");
  if (next.Threads > 1 && next.Update)
    errors.push_back("visual updates are not supported with run.threads > 1");

  const Bool_t reacted = next.stripFirst >= 0;
  if (reacted && next.NumEvapPart == 0 && next.talysClass.IsNull())
    errors.push_back(
        "reacted runs require at least one complete reaction.step");
  for (Int_t i = 0; i < next.NumEvapPart; ++i) {
    const std::string base = "reaction.step[" + std::to_string(i) + "]";
    if (reacted && !next.hasBranchTables[i] && next.residueExc[i] < 0)
      errors.push_back("missing required key '" + base +
                       ".residue_excitation'");
    if (reacted && !next.hasBranchTables[i] && next.angularDist[i] < 0)
      errors.push_back("missing required key '" + base +
                       ".angular_distribution'");
    if (reacted && i == next.NumEvapPart - 1 && next.residueExc[i] == 0)
      errors.push_back(
          "the final reaction step cannot use residue_excitation='forced'");
    // Later steps are decays of a moving residue, not scattering off the beam.
    if (reacted && i > 0 && next.angularDist[i] == 1)
      errors.push_back("'" + base +
                       ".angular_distribution' cannot be 'rutherford' after "
                       "step 0");
    for (size_t b = 0; reacted && i > 0 && next.hasBranchTables[i] &&
                       b < next.branches[i].size();
         ++b)
      if (next.branches[i][b].angularDist == 1)
        errors.push_back("'" + base + ".branch[" + std::to_string(b) +
                         "].angular_distribution' cannot be 'rutherford' "
                         "after step 0");
  }

  // A step without branch tables samples from one implicit branch built from
  // its step-level keys, so kinematics reads every step the same way.
  for (Int_t i = 0; i < next.NumEvapPart; ++i) {
    if (next.hasBranchTables[i])
      continue;
    ExcitationBranch implicit;
    implicit.residueExc = next.residueExc[i];
    implicit.angularDist = next.angularDist[i];
    implicit.thetaCmMinDeg = next.thetaCmMinDeg[i];
    implicit.weight = 1.0;
    implicit.levelTable = next.levelTable[i];
    next.branches[i].assign(1, implicit);
  }
  // A class run's chains come from talys-cache's table at each strip's centre
  // energy, so it needs a reacted run, an alpha target, and no event display.
  const Bool_t talysClass = !next.talysClass.IsNull();
  if (next.hasTalysBlock && !talysClass)
    errors.push_back("[talys] is only used with reaction.talys_class");
  if (talysClass && !reacted)
    errors.push_back("reaction.talys_class needs a reacted run "
                     "(detector.strip >= 0 or a strip range)");
  if (talysClass && next.Update)
    errors.push_back("reaction.talys_class does not support run.update = 1");
  // Each table is read from disk once; workers parse the master's copy, so a
  // table edited mid-run cannot reach part of the run.
  std::map<std::string, std::string> nextTables;
  for (Int_t i = 0; i < next.NumEvapPart; ++i) {
    for (size_t b = 0; b < next.branches[i].size(); ++b) {
      ExcitationBranch &branch = next.branches[i][b];
      if (branch.residueExc != 3 || branch.levelTable.empty())
        continue;
      std::string where = "reaction.step[" + std::to_string(i) + "]";
      if (next.hasBranchTables[i])
        where += ".branch[" + std::to_string(b) + "]";
      const std::string &path = branch.levelTable;
      if (nextTables.count(path) == 0) {
        std::string tableText;
        if (levelTables) {
          const std::map<std::string, std::string>::const_iterator found =
              levelTables->find(path);
          if (found == levelTables->end()) {
            errors.push_back(where + ".level_table: '" + path +
                             "' was not read by the master");
            continue;
          }
          tableText = found->second;
        } else if (!music::ReadFileText(path, tableText)) {
          errors.push_back(where + ".level_table: cannot open level table '" +
                           path + "'");
          continue;
        }
        nextTables[path] = tableText;
      }
      std::string error;
      if (!branch.table.Parse(nextTables[path], path, error))
        errors.push_back(where + ".level_table: " + error);
    }
  }

  const std::array<TString, 7> gases = {"4He", "He",   "helium", "3He",
                                        "Ar",  "40Ar", "argon"};
  const Bool_t simpleGas =
      std::find(gases.begin(), gases.end(), next.gas) != gases.end();
  if (!simpleGas && next.gas != "CF4" && next.gas != "CH4" &&
      next.gas != "methane" && next.gas != "P10" && next.gas != "iC4H10" &&
      next.gas != "isobutane")
    errors.push_back("unsupported gas species '" +
                     std::string(next.gas.Data()) + "'");

  // Validate the reaction graph before any physics objects are allocated.
  auto validateNuclide = [&](const TString &name, const std::string &where) {
    if (name.IsNull() || !NuF->Contains(name.Data()))
      errors.push_back("unknown nuclide '" + std::string(name.Data()) +
                       "' at " + where);
  };
  validateNuclide(next.beamName, "beam.species");
  validateNuclide(next.target, "target.species");
  validateNuclide(next.compound, "target.compound");
  for (Int_t i = 0; i < next.NumEvapPart; ++i) {
    validateNuclide(next.evap[i],
                    "reaction.step[" + std::to_string(i) + "].evap.name");
    validateNuclide(next.res[i],
                    "reaction.step[" + std::to_string(i) + "].res.name");
  }
  if (errors.empty()) {
    const Int_t beamA = NuF->GetA(next.beamName.Data());
    const Int_t beamZ = NuF->GetZ(next.beamName.Data());
    const Int_t targetA = NuF->GetA(next.target.Data());
    const Int_t targetZ = NuF->GetZ(next.target.Data());
    if (NuF->GetA(next.compound.Data()) != beamA + targetA ||
        NuF->GetZ(next.compound.Data()) != beamZ + targetZ)
      errors.push_back(
          "target.compound does not conserve A/Z for beam + target");
    TString parent = next.compound;
    for (Int_t i = 0; i < next.NumEvapPart; ++i) {
      if (NuF->GetA(parent.Data()) !=
              NuF->GetA(next.evap[i].Data()) + NuF->GetA(next.res[i].Data()) ||
          NuF->GetZ(parent.Data()) !=
              NuF->GetZ(next.evap[i].Data()) + NuF->GetZ(next.res[i].Data()))
        errors.push_back("reaction.step[" + std::to_string(i) +
                         "] does not conserve A/Z from parent '" +
                         std::string(parent.Data()) + "'");
      const Int_t evapA = NuF->GetA(next.evap[i].Data());
      const Int_t evapZ = NuF->GetZ(next.evap[i].Data());
      const Int_t resA = NuF->GetA(next.res[i].Data());
      const Int_t resZ = NuF->GetZ(next.res[i].Data());
      const Bool_t elasticAsLabelled = evapA == targetA && evapZ == targetZ &&
                                       resA == beamA && resZ == beamZ;
      const Bool_t elasticSwapped = !elasticAsLabelled && evapA == beamA &&
                                    evapZ == beamZ && resA == targetA &&
                                    resZ == targetZ;
      const Bool_t elastic = i == 0 && reacted && elasticAsLabelled;
      if (elastic && !next.hasBranchTables[i] &&
          (next.residueExc[i] != 1 || next.angularDist[i] != 1))
        errors.push_back("elastic step 0 requires residue_excitation='ground' "
                         "and angular_distribution='rutherford'");
      // REVIEW-GATED relaxation of the audit F-12 rule: a branched elastic
      // step may add inelastic branches but must keep a weighted elastic one.
      if (elastic && next.hasBranchTables[i]) {
        Bool_t elasticBranch = kFALSE;
        for (size_t b = 0; b < next.branches[i].size(); ++b) {
          const ExcitationBranch &branch = next.branches[i][b];
          if (branch.residueExc == 1 && branch.angularDist == 1 &&
              branch.weight > 0.0)
            elasticBranch = kTRUE;
        }
        if (!elasticBranch)
          errors.push_back("elastic step 0 requires a branch with "
                           "residue_excitation='ground', "
                           "angular_distribution='rutherford', and weight > 0");
      }
      // theta_cm is the res direction, so the Rutherford deflection is only
      // the beam's when res is the beam-like product.
      if (i == 0 && reacted && elasticSwapped)
        errors.push_back("elastic step 0 must list the beam-like product as "
                         "'res' and the target-like product as 'evap'");
      parent = next.res[i];
    }
    // TALYS runs alpha on the beam nucleus, and the class is an element its
    // residues can reach; elastic is in the beam's class only.
    if (talysClass) {
      next.classZ = ElementZ(*NuF, next.talysClass.Data());
      if (targetA != 4 || targetZ != 2)
        errors.push_back("reaction.talys_class needs target.species = '4He'");
      if (next.classZ < 1 || next.classZ > beamZ + targetZ)
        errors.push_back("'reaction.talys_class' = '" +
                         std::string(next.talysClass.Data()) +
                         "' is not an element at or below the compound's");
      if (tbl.at_path("reaction.theta_cm_min_deg") && next.classZ != beamZ)
        errors.push_back("'reaction.theta_cm_min_deg' sets the elastic "
                         "cut, but elastic is only in the beam's class");
    }
  }

  if (!errors.empty()) {
    std::cerr << "musicsim ERROR: invalid control file '" << sourcePath
              << "':" << JoinErrors(errors) << std::endl;
    return 0;
  }

  nextPhysics.stoppingModel = next.stoppingModel;
  const std::string output(next.FileName.Data());
  const size_t slash = output.find_last_of('/');
  nextPhysics.srimDir =
      slash == std::string::npos ? "." : output.substr(0, slash);
  char tag[128];
  std::snprintf(tag, sizeof(tag), "%s_%gTorr_%gK", next.gas.Data(),
                Double_t(next.pressure), Double_t(next.temperature));
  nextPhysics.srimGasTag = tag;

  ctf = std::move(next);
  physics_ = std::move(nextPhysics);
  ctrlFilePath_ = sourcePath;
  ctrlFileText_ = text;
  levelTableText_ = std::move(nextTables);
  masterTables_ = levelTables;
  SeedRandom(static_cast<UInt_t>(ctf.Seed));
  return 1;
}

void Simulator::InitCTF() { ctf = controlFileParams(); }

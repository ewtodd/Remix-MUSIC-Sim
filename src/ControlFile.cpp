#include "EnergyLoss.hpp"
#include "Simulator.hpp"

#include <cmath>
#include <cstdint>
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

std::string JoinErrors(const std::vector<std::string> &errors) {
  std::ostringstream out;
  for (const auto &error : errors)
    out << "\n  - " << error;
  return out.str();
}

} // namespace

Int_t Simulator::loadCtrlFile(const char *fileName) {
  if (!fileName || fileName[0] == '\0') {
    std::cerr << "musicsim ERROR: a control-file path is required."
              << std::endl;
    return 0;
  }

  toml::table tbl;
  try {
    tbl = toml::parse_file(fileName);
  } catch (const toml::parse_error &error) {
    std::cerr << "musicsim ERROR: failed to parse TOML control file '"
              << fileName << "': " << error.description() << " at "
              << error.source().begin << std::endl;
    return 0;
  }

  controlFileParams next;
  music::PhysicsConfig nextPhysics;
  nextPhysics.mean.z_effective = catima::z_eff_type::atima14;
  nextPhysics.mean.low_energy = catima::low_energy_types::srim_95;
  nextPhysics.straggling = nextPhysics.mean;
  nextPhysics.straggling.z_effective = catima::z_eff_type::pierce_blann;
  std::vector<std::string> errors;

  CheckKeys(tbl,
            {"gas", "beam", "target", "windows", "detector", "reaction",
             "physics", "run"},
            "", errors);

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
      if (code < 0 || code == catima::z_eff_type::atima14)
        errors.push_back("'physics.straggling_z_effective' must name a "
                         "straggling-capable catima charge model");
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
    CheckKeys(*reaction, {"step"}, "reaction.", errors);
    const auto stepsNode = tbl.at_path("reaction.step");
    const auto *steps = stepsNode.as_array();
    if (!steps) {
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
                   "theta_cm_min_deg"},
                  base + ".", errors);
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
          if (*value == "forced")
            next.residueExc[index] = 0;
          else if (*value == "ground")
            next.residueExc[index] = 1;
          else if (*value == "uniform")
            next.residueExc[index] = 2;
          else
            errors.push_back(
                "'" + base +
                ".residue_excitation' must be 'forced', 'ground', or 'uniform'");
        } else if ((*step)["residue_excitation"])
          errors.push_back("'" + base +
                           ".residue_excitation' must be a string");

        if (auto value = (*step)["angular_distribution"].value<std::string>()) {
          if (*value == "isotropic")
            next.angularDist[index] = 0;
          else if (*value == "rutherford")
            next.angularDist[index] = 1;
          else
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
  if (reacted && next.NumEvapPart == 0)
    errors.push_back(
        "reacted runs require at least one complete reaction.step");
  for (Int_t i = 0; i < next.NumEvapPart; ++i) {
    const std::string base = "reaction.step[" + std::to_string(i) + "]";
    if (reacted && next.residueExc[i] < 0)
      errors.push_back("missing required key '" + base +
                       ".residue_excitation'");
    if (reacted && next.angularDist[i] < 0)
      errors.push_back("missing required key '" + base +
                       ".angular_distribution'");
    if (reacted && i == next.NumEvapPart - 1 && next.residueExc[i] == 0)
      errors.push_back(
          "the final reaction step cannot use residue_excitation='forced'");
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
      if (i == 0 && NuF->GetA(next.evap[i].Data()) == targetA &&
          NuF->GetZ(next.evap[i].Data()) == targetZ &&
          NuF->GetA(next.res[i].Data()) == beamA &&
          NuF->GetZ(next.res[i].Data()) == beamZ && reacted &&
          (next.residueExc[i] != 1 || next.angularDist[i] != 1))
        errors.push_back("elastic step 0 requires residue_excitation='ground' "
                         "and angular_distribution='rutherford'");
      parent = next.res[i];
    }
  }

  if (!errors.empty()) {
    std::cerr << "musicsim ERROR: invalid control file '" << fileName
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
  ctrlFilePath_ = fileName;
  SeedRandom(static_cast<UInt_t>(ctf.Seed));
  return 1;
}

void Simulator::InitCTF() { ctf = controlFileParams(); }

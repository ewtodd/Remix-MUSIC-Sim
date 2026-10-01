// Convert the original whitespace-delimited .msc controls into the strict
// Remix-MUSIC-Sim TOML schema. The converter is dependency-free so it works in
// the project's Nix development shell without Python.

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct Field {
  std::string value;
  std::string comment;
};

using Fields = std::map<std::string, Field>;

struct Step {
  Fields evap;
  Fields residue;
};

struct Conversion {
  Fields gas;
  Fields beam;
  Fields target;
  Fields detector;
  Fields run;
  std::array<Fields, 3> layers;
  std::map<int, Step> steps;
  std::vector<std::string> reviewNotes;
  std::vector<std::string> errors;
};

std::string Trim(std::string text) {
  const size_t first = text.find_first_not_of(" \t\r\n");
  if (first == std::string::npos)
    return {};
  const size_t last = text.find_last_not_of(" \t\r\n");
  return text.substr(first, last - first + 1);
}

std::string TomlString(const std::string &text) {
  std::ostringstream out;
  out << '"';
  for (const unsigned char c : text) {
    switch (c) {
    case '\\':
      out << "\\\\";
      break;
    case '"':
      out << "\\\"";
      break;
    case '\b':
      out << "\\b";
      break;
    case '\t':
      out << "\\t";
      break;
    case '\n':
      out << "\\n";
      break;
    case '\f':
      out << "\\f";
      break;
    case '\r':
      out << "\\r";
      break;
    default:
      if (c < 0x20 || c == 0x7f)
        out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
            << static_cast<unsigned int>(c) << std::dec << std::setfill(' ');
      else
        out << static_cast<char>(c);
    }
  }
  out << '"';
  return out.str();
}

std::string CleanComment(std::string comment) {
  std::replace(comment.begin(), comment.end(), '\t', ' ');
  std::replace(comment.begin(), comment.end(), '\r', ' ');
  std::replace(comment.begin(), comment.end(), '\n', ' ');
  return Trim(comment);
}

std::optional<std::string> IntegerValue(const std::string &raw) {
  errno = 0;
  char *end = nullptr;
  const long long value = std::strtoll(raw.c_str(), &end, 10);
  if (errno == ERANGE || end == raw.c_str() || *end != '\0')
    return std::nullopt;
  return std::to_string(value);
}

std::optional<std::string> FloatValue(const std::string &raw) {
  errno = 0;
  char *end = nullptr;
  const double value = std::strtod(raw.c_str(), &end);
  if (errno == ERANGE || end == raw.c_str() || *end != '\0' ||
      !std::isfinite(value))
    return std::nullopt;
  std::ostringstream out;
  out << std::setprecision(17) << value;
  std::string rendered = out.str();
  if (rendered.find_first_of(".eE") == std::string::npos)
    rendered += ".0";
  return rendered;
}

enum class Kind { Integer, Float, String };

std::optional<std::string> RenderValue(const std::string &raw, Kind kind) {
  if (kind == Kind::String)
    return TomlString(raw);
  return kind == Kind::Integer ? IntegerValue(raw) : FloatValue(raw);
}

void Set(Fields &fields, const std::string &key, const std::string &value,
         const std::string &comment = {}) {
  fields[key] = {value, CleanComment(comment)};
}

bool ParseIndexed(const std::string &key, const std::string &prefix,
                  const std::string &suffix, int &index) {
  if (key.size() <= prefix.size() + suffix.size() ||
      key.compare(0, prefix.size(), prefix) != 0 ||
      key.compare(key.size() - suffix.size(), suffix.size(), suffix) != 0)
    return false;
  const std::string middle =
      key.substr(prefix.size(), key.size() - prefix.size() - suffix.size());
  const auto parsed = IntegerValue(middle);
  if (!parsed)
    return false;
  index = std::stoi(*parsed);
  return index >= 0 && index < 10;
}

Conversion ParseLegacy(const std::filesystem::path &source) {
  Conversion result;
  std::ifstream input(source);
  if (!input) {
    result.errors.push_back("cannot open input '" + source.string() + "'");
    return result;
  }

  struct Mapping {
    const char *section;
    const char *key;
    Kind kind;
  };
  const std::unordered_map<std::string, Mapping> scalar = {
      {"Gas", {"gas", "species", Kind::String}},
      {"Pressure", {"gas", "pressure", Kind::Float}},
      {"pressure", {"gas", "pressure", Kind::Float}},
      {"Temperature", {"gas", "temperature", Kind::Float}},
      {"temperature", {"gas", "temperature", Kind::Float}},
      {"beam", {"beam", "species", Kind::String}},
      {"BeamEnergy", {"beam", "energy", Kind::Float}},
      {"Kb", {"beam", "energy", Kind::Float}},
      {"Ebeam", {"beam", "energy", Kind::Float}},
      {"KbFWHM", {"beam", "energy_fwhm", Kind::Float}},
      {"BeamEnergyFWHM", {"beam", "energy_fwhm", Kind::Float}},
      {"EbeamFWHM", {"beam", "energy_fwhm", Kind::Float}},
      {"dEdxScaleBeam", {"beam", "dedx_scale", Kind::Float}},
      {"target", {"target", "species", Kind::String}},
      {"compound", {"target", "compound", Kind::String}},
      {"ELossBins", {"detector", "eloss_bins", Kind::Integer}},
      {"MaxELoss", {"detector", "max_eloss", Kind::Float}},
      {"strip", {"detector", "strip", Kind::Integer}},
      {"stripFirst", {"detector", "strip_first", Kind::Integer}},
      {"stripLast", {"detector", "strip_last", Kind::Integer}},
      {"Eres", {"detector", "resolution_fwhm_percent", Kind::Float}},
      {"NEvents", {"run", "events_per_strip", Kind::Integer}},
      {"Threads", {"run", "threads", Kind::Integer}},
      {"Wait", {"run", "wait", Kind::Integer}},
      {"Update", {"run", "update", Kind::Integer}},
      {"MaxTime", {"run", "max_time", Kind::Float}},
      {"SimStep", {"run", "sim_step", Kind::Float}},
      {"Method", {"run", "method", Kind::Integer}},
      {"FileName", {"run", "output", Kind::String}},
      {"FileOpt", {"run", "file_opt", Kind::String}},
      {"PrintOpt", {"run", "print_opt", Kind::Integer}},
  };
  struct LayerMapping {
    size_t layer;
    const char *key;
    Kind kind;
  };
  const std::unordered_map<std::string, LayerMapping> layerMap = {
      {"EntranceMaterial", {0, "material", Kind::String}},
      {"EntranceThickness", {0, "thickness_mg_cm2", Kind::Float}},
      {"ExitMaterial", {1, "material", Kind::String}},
      {"ExitThickness", {1, "thickness_mg_cm2", Kind::Float}},
      {"DegraderMaterial", {2, "material", Kind::String}},
      {"DegraderLength", {2, "thickness_um", Kind::Float}},
  };
  const std::set<std::string> dropped = {
      "NumEvapPart", "AnodeGeom", "SRIMdir",   "CSVfile",   "SRIMbeam",
      "SRIMevap0",   "SRIMevap1", "SRIMevap2", "SRIMevap3", "SRIMres0",
      "SRIMres1",    "SRIMres2",  "SRIMres3",  "reacClass", "PARAMETER",
      "Parameter",   "Control"};

  std::string line;
  size_t lineNumber = 0;
  while (std::getline(input, line)) {
    ++lineNumber;
    if (lineNumber <= 2 || Trim(line).empty())
      continue;
    std::istringstream tokens(line);
    std::string key;
    std::string raw;
    if (!(tokens >> key >> raw)) {
      result.errors.push_back(source.string() + ':' +
                              std::to_string(lineNumber) +
                              ": expected PARAMETER VALUE");
      continue;
    }
    std::string comment;
    std::getline(tokens, comment);
    comment = CleanComment(comment);

    const auto scalarIt = scalar.find(key);
    if (scalarIt != scalar.end()) {
      const auto value = RenderValue(raw, scalarIt->second.kind);
      if (!value) {
        result.errors.push_back(source.string() + ':' +
                                std::to_string(lineNumber) + ": invalid " +
                                key);
        continue;
      }
      Fields *fields =
          scalarIt->second.section == std::string("gas")      ? &result.gas
          : scalarIt->second.section == std::string("beam")   ? &result.beam
          : scalarIt->second.section == std::string("target") ? &result.target
          : scalarIt->second.section == std::string("detector")
              ? &result.detector
              : &result.run;
      if ((key == "Kb" || key == "Ebeam") && !comment.empty())
        comment += " -- ";
      if (key == "Kb" || key == "Ebeam")
        comment += "REVIEW: legacy gas-surface KE is now accelerator KE";
      Set(*fields, scalarIt->second.key, *value, comment);
      continue;
    }

    const auto layerIt = layerMap.find(key);
    if (layerIt != layerMap.end()) {
      const auto value = RenderValue(raw, layerIt->second.kind);
      if (!value)
        result.errors.push_back(source.string() + ':' +
                                std::to_string(lineNumber) + ": invalid " +
                                key);
      else
        Set(result.layers[layerIt->second.layer], layerIt->second.key, *value,
            comment);
      continue;
    }

    int index = -1;
    bool matched = false;
    for (const auto &[prefix, side] :
         std::array<std::pair<std::string, bool>, 2>{
             {{"evap", true}, {"res", false}}}) {
      for (const auto &[suffix, field, kind] :
           std::array<std::tuple<std::string, std::string, Kind>, 2>{
               {{"Name", "name", Kind::String},
                {"Color", "color", Kind::Integer}}}) {
        if (!ParseIndexed(key, prefix, suffix, index))
          continue;
        const auto value = RenderValue(raw, kind);
        if (!value)
          result.errors.push_back(source.string() + ':' +
                                  std::to_string(lineNumber) + ": invalid " +
                                  key);
        else
          Set(side ? result.steps[index].evap : result.steps[index].residue,
              field, *value, comment);
        matched = true;
        break;
      }
      if (matched)
        break;
    }
    if (!matched) {
      for (const auto &[prefix, side] :
           std::array<std::pair<std::string, bool>, 2>{
               {{"dEdxScaleEvap", true}, {"dEdxScaleRes", false}}}) {
        if (!ParseIndexed(key, prefix, "", index))
          continue;
        const auto value = FloatValue(raw);
        if (!value)
          result.errors.push_back(source.string() + ':' +
                                  std::to_string(lineNumber) + ": invalid " +
                                  key);
        else
          Set(side ? result.steps[index].evap : result.steps[index].residue,
              "dedx_scale", *value, comment);
        matched = true;
        break;
      }
    }
    if (matched)
      continue;
    if (dropped.count(key)) {
      if (key != "PARAMETER" && key != "Parameter" && key != "Control")
        result.reviewNotes.push_back("dropped legacy key " + key + " = " + raw);
      continue;
    }
    result.reviewNotes.push_back("unrecognized legacy key " + key + " = " +
                                 raw);
  }

  auto defaultField = [](Fields &fields, const std::string &key,
                         const std::string &value,
                         const std::string &comment = {}) {
    if (!fields.count(key))
      Set(fields, key, value, comment);
  };
  defaultField(result.gas, "pressure", "760.0", "legacy default [Torr]");
  defaultField(result.gas, "temperature", "293.0", "legacy default [K]");
  defaultField(result.beam, "energy_fwhm", "0.0", "legacy default [MeV]");
  defaultField(result.beam, "dedx_scale", "1.0");
  defaultField(result.detector, "eloss_bins", "300");
  defaultField(result.detector, "max_eloss", "10.0", "MeV");
  defaultField(result.run, "events_per_strip", "10");
  defaultField(result.run, "threads", "1");
  defaultField(result.run, "wait", "0");
  defaultField(result.run, "update", "0");
  defaultField(result.run, "max_time", "2000.0", "ns");
  defaultField(result.run, "sim_step", "0.001", "cm");
  defaultField(result.run, "method", "0");
  defaultField(result.run, "file_opt", TomlString("recreate"));
  defaultField(result.run, "print_opt", "0");
  defaultField(result.run, "seed", "1",
               "REVIEW: deterministic converter default");
  defaultField(result.layers[0], "material", TomlString("Ti"));
  if (!result.layers[0].count("thickness_mg_cm2") &&
      !result.layers[0].count("thickness_um"))
    Set(result.layers[0], "thickness_mg_cm2", "-1.0", "disabled");
  defaultField(result.layers[1], "material", TomlString("Ti"));
  if (!result.layers[1].count("thickness_mg_cm2") &&
      !result.layers[1].count("thickness_um"))
    Set(result.layers[1], "thickness_mg_cm2", "-1.0", "disabled");
  if (!result.layers[2].empty() &&
      !result.layers[2].count("thickness_mg_cm2") &&
      !result.layers[2].count("thickness_um"))
    Set(result.layers[2], "thickness_um", "-1.0", "disabled");

  if (result.run.at("method").value != "0")
    result.errors.push_back(
        "legacy trace-database Method is unsupported; expected Method 0");

  for (const auto &[section, key] :
       std::array<std::pair<const Fields *, const char *>, 9>{
           {{&result.gas, "species"},
            {&result.beam, "species"},
            {&result.beam, "energy"},
            {&result.target, "species"},
            {&result.target, "compound"},
            {&result.detector, "eloss_bins"},
            {&result.detector, "max_eloss"},
            {&result.run, "events_per_strip"},
            {&result.run, "output"}}}) {
    if (!section->count(key))
      result.errors.push_back("missing required legacy value for '" +
                              std::string(key) + "'");
  }
  const bool hasStrip = result.detector.count("strip") != 0;
  const bool hasFirst = result.detector.count("strip_first") != 0;
  const bool hasLast = result.detector.count("strip_last") != 0;
  if ((hasStrip && (hasFirst || hasLast)) ||
      (!hasStrip && hasFirst != hasLast) || (!hasStrip && !hasFirst))
    result.errors.push_back(
        "legacy control must select one strip or one complete range");

  int expectedIndex = 0;
  for (auto &[index, step] : result.steps) {
    if (index != expectedIndex++)
      result.errors.push_back(
          "reaction step indices must be contiguous from zero");
    if (!step.evap.count("name") || !step.residue.count("name"))
      result.errors.push_back("reaction step " + std::to_string(index) +
                              " is missing an evap/res name");
    defaultField(step.evap, "color", "616");
    defaultField(step.evap, "dedx_scale", "1.0");
    defaultField(step.residue, "color", "416");
    defaultField(step.residue, "dedx_scale", "1.0");
  }
  return result;
}

void EmitFields(std::ostream &out, const std::string &heading,
                const Fields &fields, const std::vector<std::string> &order) {
  out << '[' << heading << "]\n";
  size_t width = 0;
  for (const auto &key : order)
    if (fields.count(key))
      width = std::max(width, key.size());
  for (const auto &key : order) {
    const auto found = fields.find(key);
    if (found == fields.end())
      continue;
    out << std::left << std::setw(static_cast<int>(width)) << key << " = "
        << found->second.value;
    if (!found->second.comment.empty())
      out << " # " << found->second.comment;
    out << '\n';
  }
  out << '\n';
}

std::string Render(const Conversion &conversion) {
  std::ostringstream out;
  EmitFields(out, "gas", conversion.gas,
             {"species", "pressure", "temperature"});
  EmitFields(out, "beam", conversion.beam,
             {"species", "energy", "energy_fwhm", "dedx_scale"});
  EmitFields(out, "target", conversion.target, {"species", "compound"});
  EmitFields(out, "windows.entrance", conversion.layers[0],
             {"material", "thickness_mg_cm2", "thickness_um", "dedx_scale"});
  EmitFields(out, "windows.exit", conversion.layers[1],
             {"material", "thickness_mg_cm2", "thickness_um", "dedx_scale"});
  if (!conversion.layers[2].empty())
    EmitFields(out, "windows.degrader", conversion.layers[2],
               {"material", "thickness_mg_cm2", "thickness_um", "dedx_scale"});
  EmitFields(out, "detector", conversion.detector,
             {"eloss_bins", "max_eloss", "strip", "strip_first", "strip_last",
              "noise_sigma_mev", "resolution_fwhm_percent"});

  const std::string beam = conversion.beam.at("species").value;
  const std::string target = conversion.target.at("species").value;
  for (auto it = conversion.steps.begin(); it != conversion.steps.end(); ++it) {
    const bool last = std::next(it) == conversion.steps.end();
    const bool elastic = it == conversion.steps.begin() &&
                         it->second.evap.at("name").value == target &&
                         it->second.residue.at("name").value == beam;
    out << "[[reaction.step]]\n"
        << "residue_excitation = "
        << TomlString(elastic || last ? "ground" : "forced")
        << " # REVIEW: inferred by converter\n"
        << "angular_distribution = "
        << TomlString(elastic ? "rutherford" : "isotropic")
        << " # REVIEW: inferred by converter\n";
    if (elastic)
      out << "theta_cm_min_deg = 1.0 # REVIEW: Rutherford cutoff\n";
    out << '\n';
    EmitFields(out, "reaction.step.evap", it->second.evap,
               {"name", "color", "dedx_scale"});
    EmitFields(out, "reaction.step.res", it->second.residue,
               {"name", "color", "dedx_scale"});
  }

  out << "[physics]\n"
      << "z_effective = \"atima14\"\n"
      << "low_energy = \"srim_95\"\n"
      << "straggling_z_effective = \"pierce_blann\"\n"
      << "straggling = true\n"
      << "stopping = \"catima\"\n\n";
  EmitFields(out, "run", conversion.run,
             {"events_per_strip", "threads", "wait", "update", "max_time",
              "sim_step", "method", "output", "file_opt", "print_opt", "seed"});
  if (!conversion.reviewNotes.empty()) {
    out << "# Conversion notes; these are comments so the strict loader still accepts the file.\n";
    for (const auto &note : conversion.reviewNotes)
      out << "# REVIEW: " << CleanComment(note) << '\n';
  }
  return out.str();
}

void Usage(std::ostream &out) {
  out << "usage: legacy-msc-to-toml [--force] [--keep|--remove-source] "
         "[--stdout] [-o FILE] INPUT...\n";
}

} // namespace

int main(int argc, char **argv) {
  bool force = false;
  bool keepRequested = false;
  bool removeSource = false;
  bool stdoutMode = false;
  std::optional<std::filesystem::path> explicitOutput;
  std::vector<std::filesystem::path> inputs;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      Usage(std::cout);
      return 0;
    }
    if (arg == "--force") {
      force = true;
    } else if (arg == "--keep") {
      keepRequested = true;
    } else if (arg == "--remove-source") {
      removeSource = true;
    } else if (arg == "--stdout") {
      stdoutMode = true;
    } else if (arg == "--out" || arg == "-o") {
      if (++i >= argc) {
        Usage(std::cerr);
        return 2;
      }
      explicitOutput = argv[i];
    } else if (!arg.empty() && arg[0] == '-') {
      std::cerr << "legacy-msc-to-toml: unknown option '" << arg << "'\n";
      return 2;
    } else {
      inputs.emplace_back(arg);
    }
  }
  if (keepRequested && removeSource) {
    std::cerr << "legacy-msc-to-toml: --keep and --remove-source are mutually "
                 "exclusive\n";
    return 2;
  }
  if (inputs.empty() || (explicitOutput && inputs.size() != 1) ||
      (explicitOutput && stdoutMode)) {
    Usage(std::cerr);
    return 2;
  }

  for (const auto &input : inputs) {
    std::error_code error;
    if (!std::filesystem::is_regular_file(input, error)) {
      std::cerr << "legacy-msc-to-toml: input does not exist or is not a file: "
                << input << '\n';
      return 1;
    }
  }

  for (const auto &input : inputs) {
    Conversion conversion = ParseLegacy(input);
    if (!conversion.errors.empty()) {
      std::cerr << "legacy-msc-to-toml: cannot convert " << input << ":\n";
      for (const auto &error : conversion.errors)
        std::cerr << "  - " << error << '\n';
      return 1;
    }
    const std::string outputText = Render(conversion);
    if (stdoutMode) {
      if (inputs.size() > 1)
        std::cout << "# === " << input.string() << " ===\n";
      std::cout << outputText;
      continue;
    }

    std::filesystem::path output =
        explicitOutput ? *explicitOutput : input.parent_path() / input.stem();
    if (!explicitOutput)
      output.replace_extension(".toml");
    std::error_code error;
    if (std::filesystem::exists(output, error) && !force) {
      std::cerr << "legacy-msc-to-toml: refusing to overwrite " << output
                << "; pass --force to replace it\n";
      return 1;
    }
    std::filesystem::path temporary = output;
    temporary += ".tmp";
    if (std::filesystem::exists(temporary, error)) {
      std::cerr << "legacy-msc-to-toml: refusing to overwrite stale temporary "
                << temporary << '\n';
      return 1;
    }
    {
      std::ofstream stream(temporary, std::ios::binary);
      stream << outputText;
      if (!stream) {
        std::cerr << "legacy-msc-to-toml: failed writing " << temporary << '\n';
        return 1;
      }
    }
    if (force)
      std::filesystem::remove(output, error);
    error.clear();
    std::filesystem::rename(temporary, output, error);
    if (error) {
      std::cerr << "legacy-msc-to-toml: failed publishing " << output << ": "
                << error.message() << '\n';
      return 1;
    }
    std::cout << "wrote " << output << '\n';
    if (removeSource) {
      std::filesystem::remove(input, error);
      if (error) {
        std::cerr << "legacy-msc-to-toml: wrote output but could not remove "
                  << input << ": " << error.message() << '\n';
        return 1;
      }
      std::cout << "removed " << input << '\n';
    }
  }
  return 0;
}

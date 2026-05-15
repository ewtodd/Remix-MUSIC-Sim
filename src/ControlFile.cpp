#include "Simulator.hpp"

// TOML control-file parser. Schema (every section optional; missing keys keep
// the defaults set in InitCTF / controlFileParams):
//
//   [gas]      species, pressure (Torr), temperature (K)
//   [beam]     species, energy (MeV), energy_fwhm (MeV), dedx_scale
//   [target]   species, compound
//   [windows.entrance], [windows.exit], [windows.degrader]
//                 material, and exactly one of thickness_mg_cm2 / thickness_um
//   [detector] eloss_bins, max_eloss, strip OR (strip_first, strip_last), eres
//   [[reaction.step]]  evap = {name, color, dedx_scale},
//                      res  = {name, color, dedx_scale}
//   [physics]  z_effective, low_energy   (catima Config knobs)
//   [run]      n_events, threads, wait, update, max_time, sim_step, method,
//              output, file_opt, print_opt, reac_class

namespace {

// Map string -> catima::z_eff_type. Returns -1 if unrecognised.
Int_t parseZeffName(const TString &s) {
  if (s == "none")
    return catima::z_eff_type::none;
  if (s == "pierce_blann")
    return catima::z_eff_type::pierce_blann;
  if (s == "anthony_landorf")
    return catima::z_eff_type::anthony_landorf;
  if (s == "hubert")
    return catima::z_eff_type::hubert;
  if (s == "winger")
    return catima::z_eff_type::winger;
  if (s == "schiwietz")
    return catima::z_eff_type::schiwietz;
  if (s == "global")
    return catima::z_eff_type::global;
  if (s == "atima14")
    return catima::z_eff_type::atima14;
  return -1;
}

// Map string -> catima::low_energy_types. Returns -1 if unrecognised.
Int_t parseLowEnergyName(const TString &s) {
  if (s == "srim_85")
    return catima::low_energy_types::srim_85;
  if (s == "srim_95")
    return catima::low_energy_types::srim_95;
  return -1;
}

} // namespace

Int_t Simulator::loadCtrlFile(char *fileName) {
  ctrlFilePath_ = fileName ? fileName : "";

  // Our defaults differ from catima's library defaults: ATIMA-14 effective
  // charge and SRIM-95 low-energy tables match what LISE++ reports. The user
  // can override via [physics] or revert to catima's compiled defaults.
  catima::default_config.z_effective = catima::z_eff_type::atima14;
  catima::default_config.low_energy = catima::low_energy_types::srim_95;

  toml::table tbl;
  try {
    tbl = toml::parse_file(fileName);
  } catch (const toml::parse_error &err) {
    std::cerr << "musicsim ERROR: failed to parse TOML control file '"
              << fileName << "': " << err.description() << " at "
              << err.source().begin << std::endl;
    return 0;
  }

  // toml++ silently coerces ints to floats where it makes sense, so a key
  // written as `1` still reads back as 1.0 if the target is a Double_t.
  auto getString = [&](const char *section, const char *key, TString &dest) {
    if (auto v =
            tbl.at_path(std::string(section) + "." + key).value<std::string>())
      dest = *v;
  };
  auto getDouble = [&](const char *section, const char *key, auto &dest) {
    using T = std::remove_reference_t<decltype(dest)>;
    if (auto v =
            tbl.at_path(std::string(section) + "." + key).value<Double_t>())
      dest = T(*v);
  };
  auto getInt = [&](const char *section, const char *key, Int_t &dest) {
    if (auto v = tbl.at_path(std::string(section) + "." + key).value<int64_t>())
      dest = Int_t(*v);
  };
  auto getInlineString = [&](toml::node_view<toml::node> node, const char *key,
                             TString &dest) {
    if (auto v = node[key].value<std::string>())
      dest = *v;
  };
  auto getInlineDouble = [&](toml::node_view<toml::node> node, const char *key,
                             auto &dest) {
    using T = std::remove_reference_t<decltype(dest)>;
    if (auto v = node[key].value<Double_t>())
      dest = T(*v);
  };
  auto getInlineInt = [&](toml::node_view<toml::node> node, const char *key,
                          Int_t &dest) {
    if (auto v = node[key].value<int64_t>())
      dest = Int_t(*v);
  };

  // [gas]
  getString("gas", "species", ctf.gas);
  getDouble("gas", "pressure", ctf.pressure);
  getDouble("gas", "temperature", ctf.temperature);

  // [beam]
  getString("beam", "species", ctf.beamName);
  getDouble("beam", "energy", ctf.BeamEnergy);
  getDouble("beam", "energy_fwhm", ctf.KbFWHM);
  getDouble("beam", "dedx_scale", ctf.dEdxScaleBeam);

  // [target]
  getString("target", "species", ctf.target);
  getString("target", "compound", ctf.compound);

  // [windows.entrance], [windows.exit], [windows.degrader]
  // Each subtable accepts exactly one of thickness_mg_cm2 or thickness_um.
  if (auto w = tbl["windows"]; w.is_table()) {
    auto loadLayer = [&](toml::node_view<toml::node> sub, TString &mat,
                         Double_t &thick, Bool_t &byLen,
                         const char *layerName) {
      getInlineString(sub, "material", mat);
      Bool_t hasMg = false, hasUm = false;
      if (auto v = sub["thickness_mg_cm2"].value<Double_t>()) {
        thick = *v;
        byLen = false;
        hasMg = true;
      }
      if (auto v = sub["thickness_um"].value<Double_t>()) {
        thick = *v;
        byLen = true;
        hasUm = true;
      }
      if (hasMg && hasUm) {
        std::cerr << "musicsim ERROR: [windows." << layerName
                  << "] cannot specify both thickness_mg_cm2 and thickness_um."
                  << std::endl;
        std::exit(EXIT_FAILURE);
      }
    };
    auto entrance = w["entrance"];
    if (entrance.is_table())
      loadLayer(entrance, ctf.entranceMaterial, ctf.entranceThickness,
                ctf.entranceByLength, "entrance");
    auto exitw = w["exit"];
    if (exitw.is_table())
      loadLayer(exitw, ctf.exitMaterial, ctf.exitThickness, ctf.exitByLength,
                "exit");
    auto deg = w["degrader"];
    if (deg.is_table())
      loadLayer(deg, ctf.degraderMaterial, ctf.degraderLength,
                ctf.degraderByLength, "degrader");
  }

  // [detector]
  getInt("detector", "eloss_bins", ctf.ELossBins);
  getDouble("detector", "max_eloss", ctf.MaxELoss);
  getInt("detector", "strip", ctf.strip);
  getInt("detector", "strip_first", ctf.stripFirst);
  getInt("detector", "strip_last", ctf.stripLast);
  getDouble("detector", "eres", ctf.Eres);

  // [physics] — catima Config knobs. These mutate the process-wide
  // default_config, so the values picked here affect every dE/dx and
  // straggling call.
  //
  //   z_effective: effective-charge model as the projectile slows down.
  //     Heavy ions at <~ MeV/u are partially stripped; different formulas
  //     give noticeably different dE/dx.
  //         none, pierce_blann (catima default), anthony_landorf, hubert,
  //         winger, schiwietz, global, atima14 (our default; matches LISE++).
  //
  //   low_energy:  low-velocity table used below the Bethe–Bloch regime.
  //         srim_85 (catima default) or srim_95 (newer; our default).
  if (auto v = tbl.at_path("physics.z_effective").value<std::string>()) {
    Int_t code = parseZeffName(*v);
    if (code < 0) {
      std::cerr << "musicsim ERROR: physics.z_effective='" << *v
                << "' is not a recognised catima z_eff_type." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catima::default_config.z_effective = static_cast<UChar_t>(code);
  }
  if (auto v = tbl.at_path("physics.low_energy").value<std::string>()) {
    Int_t code = parseLowEnergyName(*v);
    if (code < 0) {
      std::cerr << "musicsim ERROR: physics.low_energy='" << *v
                << "' must be \"srim_85\" or \"srim_95\"." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    catima::default_config.low_energy = static_cast<UChar_t>(code);
  }

  // [[reaction.step]]
  if (auto steps = tbl.at_path("reaction.step"); steps.is_array()) {
    Int_t idx = 0;
    for (auto &el : *steps.as_array()) {
      if (idx >= controlFileParams::MaxNumEvapPart) {
        std::cerr << "musicsim warning: more than "
                  << controlFileParams::MaxNumEvapPart
                  << " reaction steps in TOML; truncating." << std::endl;
        break;
      }
      auto step = toml::node_view<toml::node>(&el);
      auto evap = step["evap"];
      if (evap.is_table()) {
        getInlineString(evap, "name", ctf.evap[idx]);
        getInlineInt(evap, "color", ctf.colorEvap[idx]);
        getInlineDouble(evap, "dedx_scale", ctf.dEdxScaleEvap[idx]);
      }
      auto res = step["res"];
      if (res.is_table()) {
        getInlineString(res, "name", ctf.res[idx]);
        getInlineInt(res, "color", ctf.colorRes[idx]);
        getInlineDouble(res, "dedx_scale", ctf.dEdxScaleRes[idx]);
      }
      ++idx;
    }
    ctf.NumEvapPart = idx;
  } else {
    ctf.NumEvapPart = 0;
  }

  // [run]
  getInt("run", "n_events", ctf.NEvents);
  getInt("run", "threads", ctf.Threads);
  getInt("run", "wait", ctf.Wait);
  getInt("run", "update", ctf.Update);
  getDouble("run", "max_time", ctf.MaxTime);
  getDouble("run", "sim_step", ctf.SimStep);
  getInt("run", "method", ctf.Method);
  getString("run", "output", ctf.FileName);
  getString("run", "file_opt", ctf.FileOpt);
  getInt("run", "print_opt", ctf.PrintOpt);
  getInt("run", "reac_class", ctf.reacClass);

  // Resolve and validate the reaction-strip selection. Bail with a clear
  // message instead of silently defaulting; see kStripUnset.
  const Int_t U = controlFileParams::kStripUnset;
  const Bool_t singleSet = (ctf.strip != U);
  const Bool_t firstSet = (ctf.stripFirst != U);
  const Bool_t lastSet = (ctf.stripLast != U);
  if (firstSet != lastSet) {
    std::cerr
        << "musicsim ERROR: 'stripFirst' and 'stripLast' must be set together."
        << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (singleSet && firstSet) {
    std::cerr
        << "musicsim ERROR: cannot mix 'strip' with 'stripFirst'/'stripLast'; pick one."
        << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (!singleSet && !firstSet) {
    std::cerr
        << "musicsim ERROR: must specify either 'strip' or both 'stripFirst' and 'stripLast'."
        << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if (firstSet) {
    if (ctf.stripFirst > ctf.stripLast) {
      std::cerr << "musicsim ERROR: stripFirst (" << ctf.stripFirst
                << ") must be <= stripLast (" << ctf.stripLast << ")."
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    ctf.stripFirst = ctf.strip;
    ctf.stripLast = ctf.strip;
  }
  return 1;
}

void Simulator::InitCTF() {
  ctf.gas = "4He";
  ctf.pressure = 760.0;
  ctf.temperature = 293.0;
  ctf.ELossBins = 300;
  ctf.MaxELoss = 10.0;
  ctf.beamName = "unassigned beam";
  ctf.dEdxScaleBeam = 1.0;
  ctf.target = "unassigned target";
  ctf.compound = "unassigned compound";
  ctf.NumEvapPart = ctf.MaxNumEvapPart;
  for (Int_t i = 0; i < ctf.NumEvapPart; i++) {
    ctf.res[i] = "unassigned res";
    ctf.dEdxScaleRes[i] = 1.0;
    ctf.colorRes[i] = 416;
    ctf.evap[i] = "unassigned evap";
    ctf.dEdxScaleEvap[i] = 1.0;
    ctf.colorEvap[i] = 616;
  }
  ctf.BeamEnergy = 100;
  ctf.KbFWHM = 0.0;
  ctf.entranceMaterial = "Ti";
  ctf.entranceThickness = 0.9;
  ctf.exitMaterial = "Ti";
  ctf.exitThickness = 0.9;
  ctf.degraderMaterial = "";
  ctf.degraderLength = 0.0;
  ctf.strip = controlFileParams::kStripUnset;
  ctf.stripFirst = controlFileParams::kStripUnset;
  ctf.stripLast = controlFileParams::kStripUnset;
  ctf.Eres = -1.0;
  ctf.NEvents = 10;
  ctf.Wait = 1;
  ctf.Update = 1;
  ctf.MaxTime = 2000.0;
  ctf.SimStep = 0.001;
  ctf.Method = 0;
  ctf.FileName = "";
  ctf.FileOpt = "";
  ctf.reacClass = -1;
  ctf.PrintOpt = 0;
}

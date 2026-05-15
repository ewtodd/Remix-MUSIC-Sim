#include "Simulator.hpp"

// Build a catima::Material describing the fill gas at ctf.pressure /
// ctf.temperature. Density follows the ideal gas law:
//   rho = P*M / (R*T) with R = 62363.59 cm³·Torr/(K·mol)
// so P in Torr and T in K give rho in g/cm³.
void Simulator::BuildGasMaterial() {
  gasEnabled_ = (ctf.pressure > 0.0);
  if (!gasEnabled_) {
    if (verbose_)
      std::cout << "musicsim warning: Pressure <= 0; gas energy loss disabled "
                   "(degrader/windows still apply if configured)."
                << std::endl;
    // Stub vacuum-like material so any code still touching gas_ gets a
    // well-formed catima::Material. The simulator gates the actual per-step
    // catima call by gasEnabled_, so this stub is never read in the hot path.
    gas_ = catima::Material();
    gas_.add_element(4.002602, 2, 1);
    gas_.density(1e-12);
    gas_.M(4.002602);
    return;
  }
  struct GasSpec {
    Int_t A;
    Int_t Z;
    Double_t mass_u;
    Int_t stn;
  };
  std::vector<GasSpec> components;
  Double_t molar_mass = 0.0;
  const TString &g = ctf.gas;
  if (g == "4He" || g == "He" || g == "helium") {
    components = {{4, 2, 4.002602, 1}};
    molar_mass = 4.002602;
  } else if (g == "3He") {
    components = {{3, 2, 3.0160293, 1}};
    molar_mass = 3.0160293;
  } else if (g == "Ar" || g == "40Ar" || g == "argon") {
    components = {{40, 18, 39.948, 1}};
    molar_mass = 39.948;
  } else if (g == "CF4") {
    components = {{12, 6, 12.011, 1}, {19, 9, 18.998, 4}};
    molar_mass = 12.011 + 4 * 18.998;
  } else if (g == "CH4" || g == "methane") {
    components = {{12, 6, 12.011, 1}, {1, 1, 1.008, 4}};
    molar_mass = 12.011 + 4 * 1.008;
  } else if (g == "P10") {
    // 90% Ar + 10% CH4 by volume → mole fractions = volume fractions (ideal
    // gas).
    components = {{40, 18, 39.948, 9}, {12, 6, 12.011, 1}, {1, 1, 1.008, 4}};
    molar_mass = 0.9 * 39.948 + 0.1 * (12.011 + 4 * 1.008);
  } else if (g == "iC4H10" || g == "isobutane") {
    components = {{12, 6, 12.011, 4}, {1, 1, 1.008, 10}};
    molar_mass = 4 * 12.011 + 10 * 1.008;
  } else {
    std::cout << "BuildGasMaterial ERROR: unknown gas '" << g << "'. "
              << "Supported: 4He, 3He, Ar, CF4, CH4, P10, iC4H10." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  const Double_t R_gas = 62363.59;
  Double_t rho = ctf.pressure * molar_mass / (R_gas * ctf.temperature);

  gas_ = catima::Material();
  for (const auto &c : components)
    gas_.add_element(c.mass_u, c.Z, Double_t(c.stn));
  gas_.density(rho);
  gas_.M(molar_mass);
}

// Composition + density only; callers set thickness via the helpers below.
catima::Material Simulator::LookupMaterial(const TString &name) {
  catima::Material m;
  Double_t density_g_per_cc = 0.0;
  if (name == "Ti" || name == "titanium") {
    m.add_element(47.867, 22, 1);
    density_g_per_cc = 4.506;
  } else if (name == "Al" || name == "Aluminum" || name == "aluminum") {
    m.add_element(26.982, 13, 1);
    density_g_per_cc = 2.6989;
  } else if (name == "Havar") {
    // Havar: Co 42 / Cr 19.5 / Fe 17.9 / Ni 12.7 / Mo 2.2 / W 2.7 (wt %).
    m.add_element(58.933, 27, 42.0);
    m.add_element(51.996, 24, 19.5);
    m.add_element(55.845, 26, 17.9);
    m.add_element(58.693, 28, 12.7);
    m.add_element(95.95, 42, 2.2);
    m.add_element(183.84, 74, 2.7);
    density_g_per_cc = 8.3;
  } else if (name == "Kapton") {
    // Kapton monomer C22 H10 N2 O5.
    m.add_element(12.011, 6, 22);
    m.add_element(1.008, 1, 10);
    m.add_element(14.007, 7, 2);
    m.add_element(15.999, 8, 5);
    density_g_per_cc = 1.42;
  } else if (name == "Mylar") {
    // Mylar C10 H8 O4.
    m.add_element(12.011, 6, 10);
    m.add_element(1.008, 1, 8);
    m.add_element(15.999, 8, 4);
    density_g_per_cc = 1.40;
  } else {
    std::cout << "LookupMaterial ERROR: unknown material '" << name << "'. "
              << "Supported: Ti, Al, Havar, Kapton, Mylar." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  m.density(density_g_per_cc);
  return m;
}

// Thickness in mg/cm² (areal density). For thin windows whose spec is areal
// density.
catima::Material Simulator::BuildSolidMaterial(const TString &name,
                                               Double_t thickness_mg_per_cm2) {
  catima::Material m = LookupMaterial(name);
  m.thickness(thickness_mg_per_cm2 / 1000.0); // catima thickness is g/cm²
  return m;
}

// Thickness in microns (linear length along beam axis). For bulk degraders.
catima::Material Simulator::BuildBulkMaterial(const TString &name,
                                              Double_t thickness_um) {
  catima::Material m = LookupMaterial(name);
  m.thickness_cm(thickness_um * 1e-4);
  return m;
}

void Simulator::BuildWindows() {
  auto buildLayer = [&](const TString &material, Double_t thick, Bool_t byLen) {
    return byLen ? BuildBulkMaterial(material, thick)
                 : BuildSolidMaterial(material, thick);
  };
  entranceWindowEnabled_ = (ctf.entranceThickness > 0.0);
  if (entranceWindowEnabled_)
    entranceWindow_ = buildLayer(ctf.entranceMaterial, ctf.entranceThickness,
                                 ctf.entranceByLength);
  else if (verbose_)
    std::cout << "musicsim warning: entrance window thickness <= 0; disabled."
              << std::endl;

  exitWindowEnabled_ = (ctf.exitThickness > 0.0);
  if (exitWindowEnabled_)
    exitWindow_ =
        buildLayer(ctf.exitMaterial, ctf.exitThickness, ctf.exitByLength);
  else if (verbose_)
    std::cout << "musicsim warning: exit window thickness <= 0; disabled."
              << std::endl;
}

void Simulator::BuildDegrader() {
  hasDegrader_ = (!ctf.degraderMaterial.IsNull() && ctf.degraderLength > 0);
  if (hasDegrader_)
    degrader_ =
        ctf.degraderByLength
            ? BuildBulkMaterial(ctf.degraderMaterial, ctf.degraderLength)
            : BuildSolidMaterial(ctf.degraderMaterial, ctf.degraderLength);
  else if (!ctf.degraderMaterial.IsNull() && verbose_)
    std::cout << "musicsim warning: degrader material='" << ctf.degraderMaterial
              << "' but thickness <= 0; degrader disabled." << std::endl;
}

// Mean energy out of a material, no straggling. Used for thin entrance/exit
// windows. Returns Ein for neutrals (Z <= 0).
Double_t Simulator::EnergyOutOfMaterial(Int_t A, Int_t Z, Double_t Ein_MeV,
                                        const catima::Material &mat) {
  if (Z <= 0 || A <= 0 || Ein_MeV <= 0.0)
    return Ein_MeV;
  catima::Projectile proj{Double_t(A), Double_t(Z)};
  proj.T = Ein_MeV / A;
  Double_t Eout_per_u = catima::energy_out(proj, mat);
  return Eout_per_u * A;
}

// Per-event propagation with Gaussian straggling sampled from catima's sigma_E.
// Floors Eout at 0 to keep it physical.
Double_t Simulator::EnergyThroughWithStraggling(Int_t A, Int_t Z,
                                                Double_t Ein_MeV,
                                                const catima::Material &mat) {
  if (Z <= 0 || A <= 0 || Ein_MeV <= 0.0)
    return Ein_MeV;
  catima::Projectile proj{Double_t(A), Double_t(Z)};
  proj.T = Ein_MeV / A;
  catima::Result r = catima::calculate(proj, mat);
  Double_t Eout = r.Eout * A; // catima reports MeV/u
  Double_t sigma_E = r.sigma_E * A;
  if (sigma_E > 0)
    Eout += Rdm->Gaus(0.0, sigma_E);
  return std::max(0.0, Eout);
}

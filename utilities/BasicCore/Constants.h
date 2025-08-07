// Constants.h
// Fundamental constants - pure header-only with no external dependencies

#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace rad {

// Mathematical constants
inline constexpr double PI{3.141592653589793238462643383279502884};

// Speed of light in vacuum (m/s)
inline constexpr double C{2.99792458e8};

// Elementary charge (C)
inline constexpr double QE{1.602176634e-19};

// Permittivity of free space (F/m)
inline constexpr double EPSILON0{8.8541878128e-12};

// Permeability of free space (H/m)
inline constexpr double MU0{1.25663706212e-6};

// Electron rest mass in kg
inline constexpr double ME{9.1093837015e-31};

// Electron rest mass in eV
inline constexpr double ME_EV{510998.95000};  // Computed: ME * C^2 / QE

// Classical electron radius in metres
inline constexpr double R_E{2.8179403227e-15};

// Fermi coupling constant (eV^-2)
inline constexpr double G_F{1.1663787e-23};

// Fine structure constant
inline constexpr double ALPHA{7.2973525698e-3};

// Bohr radius in m
inline constexpr double A0{5.29177210903e-11};

// Rydberg energy in eV
inline constexpr double RYDBERG_EV{13.605693122994};

// Ionisation energy of tritium in eV
inline constexpr double TRITIUM_I{13.603};

// Boltzmann constant in Joules/Kelving
inline constexpr double K_B{1.380649e-23};

}  // namespace rad

#endif
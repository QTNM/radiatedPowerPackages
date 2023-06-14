// Constants.h
// Fundamental constants not included in TMath

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "TMath.h"

namespace rad {
// Permittivity of free space
inline constexpr double EPSILON0{8.8541878128e-12};
// Permeability of free space
inline constexpr double MU0{1.25663706212e-6};

// Electron rest mass in kg
inline constexpr double ME{9.1093837015e-31};
// Electron rest mass in eV
inline constexpr double ME_EV{ME * TMath::C() * TMath::C() / TMath::Qe()};

// Classical electron radius in metres
inline constexpr double R_E{2.8179403227e-15};

// Fermi coupling constant
inline constexpr double G_F{1.1663787e-23};  // eV^{-2}

// Fine structure constant
inline constexpr double ALPHA{7.2973525698e-3};

// Bohr radius in m
inline constexpr double A0{5.29177210903e-11};

// Rydberg energy in eV
inline constexpr double RYDBERG_EV{13.605693122994};
}  // namespace rad

#endif

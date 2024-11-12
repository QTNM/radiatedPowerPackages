/*
  InelasticScatter.h
*/

#ifndef INELASTIC_SCATTER_H
#define INELASTIC_SCATTER_H

#include <cmath>

#include "BasicFunctions/Constants.h"
#include "Scattering/BaseScatter.h"
#include "TMath.h"

namespace rad {
enum class CalcType { Rudd1991, Kim1994 };
enum class Species { H, He, H2 };

class InelasticScatter : public BaseScatter {
 private:
  // Private member functions
  // Double differential cross-sections calculated using formulae from either
  // M. E. Rudd 1991 or Kim & Rudd 1994
  double BETA();

  double GAMMA();

  double G_B();

  double n();

  double A1();

  double A2();

  double A3();

  double S();

  /// @brief
  /// @param omega
  /// @param t
  /// @return
  double G2(double omega, double t);

  /// @brief
  /// @param omega
  /// @param t
  /// @return
  double G3(double omega, double t);

  /// @brief
  /// @param omega
  /// @param t
  /// @return
  double G4(double omega, double t);

  double G5();

  /// @brief
  /// @param omega
  /// @param t
  /// @param theta Angle in radians
  /// @return
  double f_BE(double omega, double t, double theta);

  /// @brief
  /// @param theta Angle in radians
  /// @return
  double f_b(double theta);

  /// @brief
  /// @param omega
  /// @param t
  /// @return
  double g_BE(double omega, double t);

  /// @brief
  /// @param t
  /// @return
  double F(double t);

  /// @brief
  /// @param omega
  /// @param t
  /// @return
  double f_1(double omega, double t);

  /// @brief
  /// @param omega
  /// @param t
  /// @return
  double G1(double omega, double t);

  /// @brief Calculate the CDF for the singly-differential cross section (in
  /// omega)
  /// @param omega Normalised outgoing electron energy
  /// @param t Normalised incident electron energy
  /// @return Cumulative distribution function for a given value of omega
  double CDF_SingleDiffXSec_W(double omega, double t);

  /// @brief Calculate the CDF for the double-differential cross section (in
  /// theta and omega)
  /// @param omega Normalised outgoing electron energy
  /// @param t Normalised incident electron energy
  /// @param theta Scattering angle of secondary electron in radians
  /// @return CDF for a given value of theta and omega
  double CDF_DoubleDiffXSec_theta(double omega, double t, double theta);

  /// @brief Return binding energy for the given atom/molecule
  /// @return Binding energy in eV
  double B();

  /// @brief Get the average KE of the electron in a subshell
  /// @return Average KE in eV
  double U();

  double N();

  double Ni();

  /// @brief Integral of generalised oscillator strength
  /// @return Integral of GOS
  double D();

  /// @brief Calculate the differential oscillator strength
  /// @param omega Normalised outgoing electron energy
  /// @return Differential oscillator strength (wrt omega)
  double DiffOscillatorStrength(double omega);

  CalcType calcType;

  Species species;

 public:
  /// @brief Override constructor
  /// @param T Incident kinetic energy in eV
  /// @param calc Calculation type tp use
  InelasticScatter(double T, CalcType calc = CalcType::Rudd1991,
                   Species spec = Species::H);

  /// @brief Calculate total inelastic cross section on atomic hydrogen
  /// @return Cross section in m^2
  double GetTotalXSec() override;

  /// @brief Calculate singly-differential cross section (in W)
  /// @param W Kinetic energy of ejected electron in eV
  double GetSingleDiffXSec_W(double W);

  /// @brief Calculate the double differential cross section
  /// @param W Kinetic energy of ejected electron in eV
  /// @param theta Scattering angle of secondary electron in radians
  /// @return Cross-section in m^2 / eV / rad
  double GetDoubleDiffXSec(double W, double theta);

  /// @brief Get a random secondary electron KE
  /// @return Energy in eV
  double GetRandomW();

  /// @brief Get a random secondary electron scattering angle
  /// @param W Secondary electron kinetic energy in eV
  /// @return Angle in radians
  double GetRandomTheta(double W);

  /// @brief Calculate kinetic energy of scattered primary
  /// @param W Secondary electron kinetic energy in eV
  /// @param theta Scattered angle of secondary in radians
  /// @return Scattered energy of primary in eV
  double GetPrimaryScatteredE(double W, double theta);

  /// @brief Calculate scattering angle of primary
  /// @param W Secondary electron kinetic energy in eV
  /// @param theta Scattered angle of secondary in radians
  /// @return Scattered angle of primary in radians
  double GetPrimaryScatteredAngle(double W, double theta);
};
}  // namespace rad

#endif
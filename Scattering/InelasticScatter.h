/*
  InelasticScatter.h
*/

#ifndef INELASTIC_SCATTER_H
#define INELASTIC_SCATTER_H

#include <cmath>

#include "BasicFunctions/Constants.h"
#include "Scattering/BaseScatter.h"

namespace rad {
enum Species { H, H2, He };

class InelasticScatter : public BaseScatter {
 private:
  // Private member functions
  // Double differential cross-sections calculated using formulae from
  // M. E. Rudd 1991

  /// Gives the binding energy for the species in question
  /// @return Binding energy in eV
  double I();

  /// @brief Average kinetic energy of electron in orbital
  /// @return Average kinetic energy in eV
  double U();

  double N();

  double Ni();

  /// @brief Integral of oscillator strength
  double D();

  double BETA();

  double GAMMA();

  double G_B();

  double n();

  double A1();

  double A2();

  double A3();

  double S();

  /// @brief
  /// @param W
  /// @param T
  /// @return
  double G2(double W, double T);

  /// @brief
  /// @param W
  /// @param T
  /// @return
  double G3(double W, double T);

  /// @brief
  /// @param W
  /// @param T
  /// @return
  double G4(double W, double T);

  double G5();

  /// @brief
  /// @param W
  /// @param T
  /// @param theta Angle in radians
  /// @return
  double f_BE(double W, double T, double theta);

  /// @brief
  /// @param theta Angle in radians
  /// @return
  double f_b(double theta);

  /// @brief
  /// @param W
  /// @param T
  /// @return
  double g_BE(double W, double T);

  /// @brief
  /// @param T
  /// @return
  double F(double T);

  /// @brief
  /// @param W
  /// @param T
  /// @return
  double f_1(double W, double T);

  /// @brief
  /// @param W
  /// @param T
  /// @return
  double G1(double W, double T);

  /// @brief
  /// @param W
  /// @param T
  /// @param theta
  /// @return
  double G4fb(double W, double T, double theta);

  /// @brief
  /// @param T
  /// @return
  double g1(double T);

  Species theSpecies;

 public:
  /// @brief Override constructor
  /// @param T Incident kinetic energy in eV
  /// @param species Species of the incident particle
  InelasticScatter(double T, Species species = H2);

  /// @brief Calculate total inelastic cross section on atomic hydrogen
  /// @return Cross section in m^2
  double GetTotalXSec() override;

  /// @brief Calculate the double differential cross section
  /// @param W Kinetic energy of ejected electron in eV
  /// @param theta Scattering angle of secondary electron in radians
  /// @return Cross-section in m^2 / eV / rad
  double GetDoubleDiffXSec(double W, double theta);

  /// @brief Single differential cross-section (in W)
  /// @param W Kinetic energy of ejected electron in eV
  /// @return Cross-section in m^2 / eV
  double GetSDCS_W(double W);

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
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
class InelasticScatter : public BaseScatter {
 private:
  // Private member functions
  // Double differenial cross-sections calculated using formulae from
  // M. E. Rudd 1991

  double BETA() { return 0.60; }

  double GAMMA() { return 10.0; }

  double G_B() { return 2.9; }

  double n() { return 2.4; }

  double A1() { return 0.74; }

  double A2() { return 0.87; }

  double A3() { return -0.60; }

  double S() {
    return 4 * TMath::Pi() * pow(A0, 2) * pow(RYDBERG_EV / TRITIUM_I, 2);
  }

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

  double G5() { return 0.33; }

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

 public:
  /// @brief Override constructor
  /// @param T Incident kinetic energy in eV
  InelasticScatter(double T) : BaseScatter(T) {}

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
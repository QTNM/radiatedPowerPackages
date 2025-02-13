/*
  ElasticScatter.h
*/

#ifndef ELASTIC_SCATTER_H
#define ELASTIC_SCATTER_H

#include <vector>

#include "Scattering/BaseScatter.h"

namespace rad {
class ElasticScatter : public BaseScatter {
 private:
  // Member data
  // Used for calculating differential cross section
  const std::vector<double> eVec{1e3,  2e3,  4e3,   8e3,  16e3,
                                 32e3, 64e3, 128e3, 256e3};  // eV
  const std::vector<double> A1Vec{-9.007e-1, -6.539e-1, -3.655e-1,
                                  -5.499e-1, -1.96e-2,  4.526e-2,
                                  -6.58e-1,  8.393e-3,  -3.739e-1};
  const std::vector<double> A2Vec{3.975e-1, 3.38e-1,  2.884e-1,
                                  3.151e-1, 2.809e-1, 2.774e-1,
                                  3.126e-1, 2.787e-1, 2.928e-1};
  const std::vector<double> A3Vec{2.344e-03, 3.208e-3, 2.94e-3,
                                  1.429e-3,  9.329e-4, 4.1e-4,
                                  3.017e-5,  1.038e-4, 1.757e-5};
  const std::vector<double> A4Vec{-3.534e-5, -1.59e-5, -5.392e-6,
                                  9.522e-6,  8.538e-7, -4.278e-8,
                                  7.506e-7,  4.492e-9, 3.551e-8};
  const std::vector<double> BVec{8.057e-3, 4.506e-3, 2.592e-3,
                                 1.872e-3, 8.431e-4, 3.444e-4,
                                 3.049e-4, 8.926e-5, 6.648e-5};
  const std::vector<double> C0Vec{1.105,    8.986e-1,  6.487e-1,
                                  8.062e-1, 1.901e-2,  -9.682e-2,
                                  9.669e-1, -1.011e-1, 4.769e-1};
  const std::vector<double> C1Vec{1.172,    1.05,      9.256e-1,
                                  9.955e-1, 2.643e-2,  -1.263e-1,
                                  1.229,    -1.410e-1, 6.287e-1};
  const std::vector<double> C2Vec{7.611e-1, 7.519e-1,  8.045e-1,
                                  7.51e-1,  2.258e-2,  -1.017e-1,
                                  9.513e-1, -1.224e-1, 5.042e-1};
  const std::vector<double> C3Vec{4.001e-1, 4.377e-1,  5.676e-1,
                                  4.597e-1, 1.605e-2,  -6.736e-2,
                                  5.969e-1, -8.834e-2, 3.282e-1};
  const std::vector<double> C4Vec{1.718e-1, 2.092e-1,  3.277e-1,
                                  2.304e-1, 9.861e-3,  -3.748e-2,
                                  3.072e-1, -5.421e-2, 1.760e-1};
  const std::vector<double> C5Vec{5.558e-2, 7.568e-2,  1.426e-1,
                                  8.723e-2, 4.891e-3,  -1.640e-2,
                                  1.202e-1, -2.652e-2, 7.261e-2};
  const std::vector<double> C6Vec{1.031e-2,  1.571e-2,  3.491e-2,
                                  1.878 - 2, 1.71e-3,   -4.697e-3,
                                  2.774e-2,  -8.267e-3, 1.820e-2};

  // Arrays used for calculating the differential cross section
  // Data from Salvat et al. Phys. Rev. A 36, 467 (1987)
  const std::array<std::array<double, 2>, 5> Ai{{{-184.39, 185.39},
                                                 {-0.2259, 1.2259},
                                                 {0.6045, 0.3955},
                                                 {0.3278, 0.6722},
                                                 {0.2327, 0.7673}}};
  const std::array<std::array<double, 2>, 5> alphai{{{2.0027, 1.9973},
                                                     {5.5272, 2.3992},
                                                     {2.8174, 0.6625},
                                                     {4.5430, 0.9852},
                                                     {5.9900, 1.2135}}};

  /// @brief Calculates the rutherford differential cross section
  /// @param theta Scattering angle in radians
  /// @return Differential cross section in units of m^2/sr
  double RutherfordDCS(double theta);

  /// @brief Calculate the momentum transfer for a given scattering angle
  /// @param theta Scattering angle in radians
  /// @return Momentum transfer in kg m/s
  double Calcq(double theta);

  /// @brief Calculate correction factor
  /// @return Correction factor
  double Calct();

  /// @brief Calculate form factor
  /// @return Form factor
  double CalcFe(double q);

  /// @brief Calculate form factor
  /// @param R Radius in m
  /// @param q Momentum transfer in kg m/s
  /// @return Form factor
  double CalcF(double R, double q);

  /// @brief Nuclear radius for a given atomic number
  /// @return Nuclear radius in m
  double R0();

  /// @brief Nuclear form factor
  /// @param q Momentum transfer in kg m/s
  /// @return Form factor
  double NuclearFormFactor(double q);

  /// @brief Calculates the screening factor for a given momentum transfer
  /// @param theta Scattering angle in radians
  /// @return Screening factor
  double ScreeningFactor(double theta);

  /// @brief Calculates the rutherford cross section on atomic hydrogen
  /// @return Cross section in m^2
  double TotalRutherfordXSec();

  /// @brief Calculate differential cross section using parametrization from
  /// Riley, MacCallum, Biggs (1975)
  /// @param A Vector of appropriate A params
  /// @param B Dimensionless B parameter
  /// @param C Vector of C parameters
  /// @param cosTheta
  /// @return Differential cross section in units of m^2/r
  double GetDiffXSec(std::vector<double> A, double B, std::vector<double> C,
                     double cosTheta);

  unsigned int Z;  // Atomic number of target

  unsigned int A;  // Atomic mass of target

 public:
  /// @brief Override constructor
  /// @param T Incident kinetic energy in eV
  /// @param aNum Atomic number of target (default is 1)
  /// @param aMass Atomic mass of target (default is 3)
  ElasticScatter(double T, unsigned int aNum = 1, unsigned int aMass = 3);

  /// @brief Total Mott cross section
  /// @return Cross section in m^2
  double GetTotalXSec() override;

  /// @brief Calculate the elastic differential cross section on hydrogen
  /// @param cosTheta Cosine of scattering angle of the electron
  /// @return Differential cross section in units of m^2 / sr
  double GetDiffXSec(double cosTheta);

  /// @brief Generates a random scattering angle according to the differential
  /// cross-section
  /// @return The scattering angle in radians
  double GetRandomScatteringAngle();
};
}  // namespace rad

#endif
/*
  QTNMFields.h

  Derived QTNM field classes

  Transcription of T. Goffrey's original Python code
*/

#ifndef QTNM_FIELDS_H
#define QTNM_FIELDS_H

#include "BasicFunctions/Constants.h"
#include "ElectronDynamics/BaseField.h"
#include "TGraph.h"
#include "TSpline.h"
#include "TVector3.h"

namespace rad {
/// Class describing a uniform magnetic field pointing in the +z direction
class UniformField : public BaseField {
 private:
  double fieldStrength;

 public:
  /// @param field The magnitude of the field in Tesla
  UniformField(double field) : fieldStrength(field) {}

  /// @param vec Position vector of charge
  /// @returns The magnetic field vector at the point
  TVector3 evaluate_field_at_point(const TVector3 vec) override;

  /// @brief Calculate electric field at a point
  /// @param v Position at which to calculate field
  /// @return Electric field = 0
  TVector3 evaluate_e_field_at_point(TVector3 v) override {
    return TVector3(0, 0, 0);
  }
};

/// Class describing the magnetic field of an infinitesimally thin current coil
/// in the x-y plane
class CoilField : public BaseField {
 private:
  double coilRadius;
  double coilCurrent;
  double coilZ;
  double coilMu;

  double central_field();

  double on_axis_field(const double z);

 public:
  /// \param radius The radius of the coil (in metres)
  /// \param current The current passing through the coil (in Amps)
  /// \param z The z position of the coil (in metres)
  /// \param Magnetic permeability
  CoilField(const double radius = 0.005, const double current = 40,
            const double z = 0.0, const double mu = MU0);

  /// \param vec Position vector of charge
  /// \returns The magnetic field vector at the point
  TVector3 evaluate_field_at_point(const TVector3 vec) override;

  /// @brief Calculate electric field at a point
  /// @param v Position at which to calculate field
  /// @return Electric field = 0
  TVector3 evaluate_e_field_at_point(TVector3 v) override {
    return TVector3(0, 0, 0);
  }
};

/// Class describing the bathtub field
/// Made up of two superposed current loops and a background magnetic field
/// The coils must be located along the z axis
class BathtubField : public BaseField {
 private:
  CoilField coil1;
  CoilField coil2;
  TVector3 btBkg;

 public:
  /// \param radius The radius of both coils
  /// \param The current (in amps) flowing through each trap coil
  /// \param Z1 The z position of the centre of the first trap coil
  /// \param Z2 The z position of the centre of the second trap coil
  /// \param background A vector describing the background field in Tesla
  BathtubField(const double radius, const double current, const double Z1,
               const double Z2, TVector3 background);

  /// \param vec Position vector of charge
  /// \returns The magnetic field vector at the point
  TVector3 evaluate_field_at_point(const TVector3 vec) override;

  /// @brief Calculate electric field at a point
  /// @param v Position at which to calculate field
  /// @return Electric field = 0
  TVector3 evaluate_e_field_at_point(TVector3 v) override {
    return TVector3(0, 0, 0);
  }
};

/// Class describing the field of a finite solenoid
/// Axis of the solenoid coincides with z axis
class SolenoidField : public BaseField {
 private:
  double r;
  double l;
  double i;
  double n;
  double mu;
  double zOff;

  /// Gives the on axis field for the solenoid
  /// \param z The position along the z axis
  /// \Returns The magnetic field in Tesla
  double on_axis_field(const double z);

 public:
  /// \param radius The internal radius of the solenoid (in metres)
  /// \param length The length of the bore of the solenoid (in metres)
  /// \param current The current (in amps) passing through each filament
  /// \param turnsPerMetre Number of coil terms per metre
  /// \param perm Permeability (default is permeability of free space)
  /// \param zOffset Offset (in metres) of the solenoid centre from z = 0
  /// (default is 0)
  SolenoidField(double radius = 0.2, double length = 0.05, double current = 1e3,
                double turnsPerMetre = 400, double perm = MU0,
                double zOffset = 0.0)
      : r(radius),
        l(length),
        i(current),
        n(turnsPerMetre),
        mu(perm),
        zOff(zOffset) {}

  /// Gives the field at a positon vector
  /// \param vec Position vector of charge
  /// \Returns The magnetic field vector at the point (in Tesla)
  TVector3 evaluate_field_at_point(const TVector3 vec) override;

  /// @brief Calculate electric field at a point
  /// @param v Position at which to calculate field
  /// @return Electric field = 0
  TVector3 evaluate_e_field_at_point(TVector3 v) override {
    return TVector3(0, 0, 0);
  }
};

/// Class describing something like a non-ideal background field, i.e. from a
/// finite solenoid The field variation is expressed as quadratic functions of
/// the axial and radial distances The axis of the field coincides with the z
/// axis
class InhomogeneousBackgroundField : public BaseField {
 private:
  double BCent;
  double inhomAx;
  double inhomAxPosition;
  double inhomRad;
  double inhomRadRadius;

 public:
  /// Parametrised constructor for the inhomogeneous background field
  /// \param centralField The magnetic field (in Tesla) at the centre of the
  /// trap \param fractionalInhomAx The maximum axial inhomogeneity as a
  /// fraction of the maximum field \param inhomZPos Z position at which the
  /// inhomogeneity reaches fractionalInhomAx \param fractionalInhomRad The
  /// maximum radial inhomogeneity as a fraction of the maximum field \param
  /// inhomRadPos Radial position at which the inhomogeneity reaches
  /// fractionalInhomRad
  InhomogeneousBackgroundField(const double centralField = 1.0,
                               const double fractionalInhomAx = 1e-6,
                               const double inhomZPos = 0.15,
                               const double fractionalInhomRad = 0.0,
                               const double inhomRadPos = 0.05)
      : BCent(centralField),
        inhomAx(fractionalInhomAx),
        inhomAxPosition(inhomZPos),
        inhomRad(fractionalInhomRad),
        inhomRadRadius(inhomRadPos) {}

  /// Gives the field at a positon vector
  /// \param vec Position vector of charge
  /// \Returns The magnetic field vector at the point (in Tesla)
  TVector3 evaluate_field_at_point(const TVector3 vec) override;

  /// @brief Calculate electric field at a point
  /// @param v Position at which to calculate field
  /// @return Electric field = 0
  TVector3 evaluate_e_field_at_point(TVector3 v) override {
    return TVector3(0, 0, 0);
  }
};

/// Class describing a bathtub trap with some inhomogeneity added
/// Coils are located along the z axis which is the same as the background field
class InhomogeneousBathtubField : public BaseField {
 private:
  CoilField coil1;
  CoilField coil2;
  InhomogeneousBackgroundField bkgField;

 public:
  /// \param radius The radius of both coils in metres
  /// The radial inhomogeneity is also measured to here
  /// \param The current (in amps) flowing through each trap coil
  /// \param Z The distance from the trap centre of the trap coils (in metres)
  /// The axial inhomogeneity is also measured to here
  /// \param centralField The central background field (in Tesla)
  /// \param fractionalInhomZ The maximum axial inhomogeneity as a fraction of
  /// the maximum field \param fractionalInhomR The maximum radial inhomogeneity
  /// as a fraction of the maximum field
  InhomogeneousBathtubField(const double radius, const double current,
                            const double Z, const double centralField,
                            const double fractionalInhomZ,
                            const double fractionalInhomR = 0.0);

  /// Get the magnetic field at a point
  /// \param vec The position vector of the charge
  /// \Returns The magnetic field vector at the point (in Tesla)
  TVector3 evaluate_field_at_point(const TVector3 vec) override;

  /// @brief Calculate electric field at a point
  /// @param v Position at which to calculate field
  /// @return Electric field = 0
  TVector3 evaluate_e_field_at_point(TVector3 v) override {
    return TVector3(0, 0, 0);
  }
};

/// Class describing a harmonic trap
/// Formed using a single coil in a background field
/// Trap coil is located at z = 0 with field anti parallel to z axis
class HarmonicField : public BaseField {
 private:
  CoilField coil;
  TVector3 btBkg;

 public:
  /// \param radius of the trap coil (in metres)
  /// \param current The current (in amps) flowing through each trap coil
  /// \param background The background field in Tesla
  HarmonicField(const double radius, const double current,
                const double background);

  /// \param vec Position vector of charge
  /// \returns The magnetic field at the point
  TVector3 evaluate_field_at_point(const TVector3 vec) override;

  /// @brief Calculate electric field at a point
  /// @param v Position at which to calculate field
  /// @return Electric field = 0
  TVector3 evaluate_e_field_at_point(TVector3 v) override {
    return TVector3(0, 0, 0);
  }
};

/// Class describing the custom HTS magnet designed for the UCL AMOPP group
class HTSMagnetUCL : public BaseField {
 private:
  TGraph* grFieldZ;
  TGraph* grFieldR;

  TSpline3* spFieldZ;
  TSpline3* spFieldR;

 public:
  /// Default constructor
  HTSMagnetUCL();

  /// Default destructor
  ~HTSMagnetUCL() override;

  /// Calculates the magnetic field at a point
  /// \param vec The position vector (in metres)
  /// \Returns The magnetic field vector at the point (in Tesla)
  TVector3 evaluate_field_at_point(const TVector3 vec) override;

  /// @brief Calculate electric field at a point
  /// @param v Position at which to calculate field
  /// @return Electric field = 0
  TVector3 evaluate_e_field_at_point(TVector3 v) override {
    return TVector3(0, 0, 0);
  }
};

class HTSMagnetTrap : public BaseField {
 private:
  CoilField coil;
  HTSMagnetUCL bkg;

 public:
  HTSMagnetTrap(double radius, double current);

  TVector3 evaluate_field_at_point(const TVector3 vec) override;

  /// @brief Calculate electric field at a point
  /// @param v Position at which to calculate field
  /// @return Electric field = 0
  TVector3 evaluate_e_field_at_point(TVector3 v) override {
    return TVector3(0, 0, 0);
  }
};

class HelmholtzField : public BaseField {
 private:
  SolenoidField coil1;
  SolenoidField coil2;

 public:
  /// @brief Generate field from Helmholtz coils
  /// @param coilRadius Radius of both coils in metres
  /// @param coilLength Length of both coils in metres
  /// @param coilCurrent Current through both coils in Amps
  /// @param turnsPerMetre Turns per metre for coils
  HelmholtzField(double coilRadius = 0.2, double coilLength = 0.05,
                 double coilCurrent = 1e3, double turnsPerMetre = 400);

  /// @brief Get magnetic field at point
  /// @param vec 3-vector at which to evaluate B field
  /// @return Magnetic field vector in Tesla
  TVector3 evaluate_field_at_point(const TVector3 vec) override;

  /// @brief Calculate electric field at a point
  /// @param v Position at which to calculate field
  /// @return Electric field = 0
  TVector3 evaluate_e_field_at_point(TVector3 v) override {
    return TVector3(0, 0, 0);
  }
};

// Class representing 4 current loops
// Useful for representing a shaped bathtub
// Two pair of coils symmetric about z = 0
class FourCoilField : public BaseField {
 private:
  CoilField coil1a;
  CoilField coil1b;
  CoilField coil2a;
  CoilField coil2b;

 public:
  /// @brief Parametrised constructor
  /// @param rPair1 Radius in metres of first coil pair
  /// @param iPair1 Current in amps of first coil pair
  /// @param zOffPair1 z offset of first coil pair
  /// @param rPair2 Radius in metres of second coil pair
  /// @param iPair2 Current in amps of second coil pair
  /// @param zOffPair1 z offset of second coil pair
  FourCoilField(double rPair1 = 15e-3, double iPair1 = 100,
                double zOffPair1 = 25e-3, double rPair2 = 15e-3,
                double iPair2 = -50, double zOffPair2 = 19e-3);

  /// @brief Get magnetic field at point
  /// @param vec 3-vector at which to evaluate B field
  /// @return Magnetic field vector in Tesla
  TVector3 evaluate_field_at_point(const TVector3 vec) override;

  /// @brief Calculate electric field at a point
  /// @param v Position at which to calculate field
  /// @return Electric field = 0
  TVector3 evaluate_e_field_at_point(TVector3 v) override {
    return TVector3(0, 0, 0);
  }
};

// Bathtub field with the above 4-coil trap
class FourCoilBathtub : public BaseField {
 private:
  FourCoilField fourCoils;
  TVector3 bkg;

 public:
  /// @brief Parametrised constructor
  FourCoilBathtub(double rPair1, double iPair1, double zOffPair1, double rPair2,
                  double iPair2, double zOffPair2, double bkgField);

  /// @brief Get magnetic field at point
  /// @param vec 3-vector at which to evaluate B field
  /// @return Magnetic field vector in Tesla
  TVector3 evaluate_field_at_point(const TVector3 vec) override;

  /// @brief Calculate electric field at a point
  /// @param v Position at which to calculate field
  /// @return Electric field = 0
  TVector3 evaluate_e_field_at_point(TVector3 v) override {
    return TVector3(0, 0, 0);
  }
};

// Bathtub field with a Helmholtz background field
class FourCoilHelmholtz : public BaseField {
 private:
  FourCoilField trapField;
  HelmholtzField bkgField;

 public:
  /// @brief Parametrised constructor
  /// @param tField Trap field
  /// @param bgField Background field
  FourCoilHelmholtz(FourCoilField tField, HelmholtzField bgField);

  /// @brief Get magnetic field at point
  /// @param vec 3-vector at which to evaluate B field
  /// @return Magnetic field vector in Tesla
  TVector3 evaluate_field_at_point(const TVector3 vec) override;

  /// @brief Calculate electric field at a point
  /// @param v Position at which to calculate field
  /// @return Electric field = 0
  TVector3 evaluate_e_field_at_point(TVector3 v) override {
    return TVector3(0, 0, 0);
  }
};

/// @brief Class representing an ideal harmonic trap
class IdealHarmonicTrap : public BaseField {
 private:
  double L0;
  double B0;

 public:
  /// @brief Parametrised constructor
  /// @param l0 Controls the steepness of the trap
  /// @param b0 The background magnetic field in tesla
  IdealHarmonicTrap(double l0, double b0) : L0(l0), B0(b0) {}

  /// @brief Get magnetic field at point
  /// @param vec 3-vector at which to evaluate B field
  /// @return Magnetic field vector in Tesla
  TVector3 evaluate_field_at_point(const TVector3 vec) override;

  /// @brief Calculate electric field at a point
  /// @param v Position at which to calculate field
  /// @return Electric field = 0
  TVector3 evaluate_e_field_at_point(TVector3 v) override {
    return TVector3(0, 0, 0);
  }
};
}  // namespace rad

#endif
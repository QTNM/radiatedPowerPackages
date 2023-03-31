/*
    ComsolFields.h

    Class used for realistic fields from COMSOL simulations
    Values come from M. Fleck's simulations

    S. Jones 26-10-2022
*/

#ifndef COMSOL_FIELDS_H
#define COMSOL_FIELDS_H

#include <string>

#include "ElectronDynamics/BaseField.h"
#include "ElectronDynamics/QTNMFields.h"
#include "TGraph2D.h"
#include "TVector3.h"

namespace rad {
class ComsolField : public BaseField {
 private:
  TGraph2D *fieldValues = 0;
  double scaleFactor;

 public:
  /// Parametrised constructor
  /// \param fieldFile CSV file path
  /// \param centralField The desired central magnetic field
  ComsolField(std::string fieldFile, double centralField = 0.0);

  /// Destructor
  ~ComsolField() override;

  /// Calculates the magnetic field at a point in space
  /// \param vec The position vector (units of metres)
  /// \return The magnetic field vector (units of tesla)
  TVector3 evaluate_field_at_point(const TVector3 vec) override;

  /// @brief Calculate electric field at a point
  /// @param v Position at which to calculate field
  /// @return Electric field = 0
  TVector3 evaluate_e_field_at_point(TVector3 v) override {
    return TVector3(0, 0, 0);
  }
};

class ComsolHarmonicField : public BaseField {
 private:
  CoilField coil;
  TGraph2D *fieldValues = 0;
  double scaleFactor;

 public:
  /// Parametrised constructor
  /// \param fieldFile CSV file path
  /// \param centralField The desired central magnetic field
  ComsolHarmonicField(double radius, double current, std::string fieldFile,
                      double centralField = 0.0);

  /// Destructor
  ~ComsolHarmonicField() override;

  /// Calculates the magnetic field at a point in space
  /// \param vec The position vector (units of metres)
  /// \return The magnetic field vector (units of tesla)
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
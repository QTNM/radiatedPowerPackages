// BaseField.cxx

#include "ElectronDynamics/BaseField.h"

double rad::BaseField::evaluate_field_magnitude(TVector3 v) {
  return evaluate_field_at_point(v).Mag();
}

double rad::BaseField::evaluate_e_field_magnitude(TVector3 v) {
  return evaluate_e_field_at_point(v).Mag();
}
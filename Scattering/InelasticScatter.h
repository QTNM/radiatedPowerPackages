/*
  InelasticScatter.h
*/

#ifndef INELASTIC_SCATTER_H
#define INELASTIC_SCATTER_H

#include "Scattering/BaseScatter.h"

namespace rad {
class InelasticScatter : public BaseScatter {
 public:
  /// @brief Calculate total inelastic cross section on atomic hydrogen
  /// @return Cross section in m^2
  double GetTotalXSec() override;
};
}  // namespace rad

#endif
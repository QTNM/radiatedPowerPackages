/*
  ElasticScatter.h
*/

#ifndef ELASTIC_SCATTER_H
#define ELASTIC_SCATTER_H

#include "Scattering/BaseScatter.h"

namespace rad {
class ElasticScatter : public BaseScatter {
 private:
  double TotalRutherfordXSec();

 public:
  double GetTotalXSec() override;
};
}  // namespace rad

#endif
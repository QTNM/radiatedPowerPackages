/*
  NuFitValues.h
  
  Values from NuFit 5.1 analysis of oscillation data
  JHEP 09 (2020) 178 [arXiv:2007.14792]
*/

#ifndef NUFIT_VALUES_H
#define NUFIT_VALUES_H

#include "Constants.h"

namespace rad
{
  inline constexpr double kNuFitDmsq21 { 7.42e-5 };
  
  // Assuming normal hierarchy
  inline constexpr double kNuFitTh12NH { 33.44 * PI / 180.0 };
  inline constexpr double kNuFitTh23NH { 49.2 * PI / 180.0 };
  inline constexpr double kNuFitTh13NH { 8.57 * PI / 180.0 };
  inline constexpr double kNuFitdCPNH  { 194 * PI / 180.0 };
  inline constexpr double kNuFitDmsq32NH { 2.515e-3 - kNuFitDmsq21 };
      
  // Assuming inverted hierarchy
  inline constexpr double kNuFitTh12IH { 33.45 * PI / 180.0 };
  inline constexpr double kNuFitTh23IH { 49.5 * PI / 180.0 };
  inline constexpr double kNuFitTh13IH { 8.60 * PI / 180.0 };
  inline constexpr double kNuFitdCPIH  { 287 * PI / 180.0 };
  inline constexpr double kNuFitDmsq32IH { -2.498e-3 };
}

#endif
/// IWaveguide.cxx

#include "Waveguides/IWaveguide.h"

#include <iostream>

#include "BasicFunctions/Constants.h"

double rad::IWaveguide::GetModeImpedance(WaveguideMode mode, double omega) {
  double Z{-1.0};
  double k_c{GetCutoffWavenumber(mode)};
  double k{omega / TMath::C()};
  double beta{sqrt(k * k - k_c * k_c)};

  if (mode.GetModeType() == ModeType::kTE) {
    Z = k * sqrt(MU0 / EPSILON0) / beta;
  } else if (mode.GetModeType() == ModeType::kTM) {
    Z = beta * sqrt(MU0 / EPSILON0) / k;
  } else {
    std::cout << "Currently we are not set up to deal with TEM modes. You will "
                 "not get a result that makes sense."
              << std::endl;
  }
  return Z;
}

bool rad::IWaveguide::ModePropagates(WaveguideMode mode, double f) {
  const double fc{GetCutoffFrequency(mode)};
  if (fc <= 0) {
    return false;
  } else {
    return f > fc;
  }
}

double rad::IWaveguide::GetPhaseVelocity(WaveguideMode mode, double f) {
  const double f_c{GetCutoffFrequency(mode)};
  return TMath::C() / sqrt(1 - pow(f_c / f, 2));
}

double rad::IWaveguide::GetGroupVelocity(WaveguideMode mode, double f) {
  const double f_c{GetCutoffFrequency(mode)};
  return TMath::C() * sqrt(1 - pow(f_c / f, 2));
}
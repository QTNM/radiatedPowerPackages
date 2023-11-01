/// IWaveguide.cxx

#include "Waveguides/IWaveguide.h"

#include <iostream>

#include "BasicFunctions/Constants.h"

double rad::IWaveguide::GetModeImpedance(Mode_t modeType, unsigned int n,
                                         unsigned int m, double omega) {
  double Z{-1.0};
  double k_c{GetCutoffWavenumber(modeType, n, m)};
  double k{omega / TMath::C()};
  double beta{sqrt(k * k - k_c * k_c)};

  if (modeType == kTE) {
    Z = k * sqrt(MU0 / EPSILON0) / beta;
  } else if (modeType == kTM) {
    Z = beta * sqrt(MU0 / EPSILON0) / k;
  } else {
    std::cout << "Currently we are not set up to deal with TEM modes. You will "
                 "not get a result that makes sense."
              << std::endl;
  }
  return Z;
}

double rad::IWaveguide::GetModeImpedance(WaveguideMode mode, double omega) {
  if (mode.GetModeType() == ModeType::kTE) {
    return GetModeImpedance(IWaveguide::kTE, mode.GetModeIndex1(),
                            mode.GetModeIndex2(), omega);
  } else if (mode.GetModeType() == ModeType::kTM) {
    return GetModeImpedance(IWaveguide::kTM, mode.GetModeIndex1(),
                            mode.GetModeIndex2(), omega);
  } else {
    return GetModeImpedance(IWaveguide::kTEM, mode.GetModeIndex1(),
                            mode.GetModeIndex2(), omega);
  }
}

TVector3 rad::IWaveguide::GetModeEField(TVector3 pos, WaveguideMode mode,
                                        double A, double omega, bool state) {
  if (mode.GetModeType() == ModeType::kTE) {
    return GetModeEField(pos, IWaveguide::kTE, A, mode.GetModeIndex1(),
                         mode.GetModeIndex2(), omega, state);
  } else if (mode.GetModeType() == ModeType::kTM) {
    return GetModeEField(pos, IWaveguide::kTM, A, mode.GetModeIndex1(),
                         mode.GetModeIndex2(), omega, state);
  } else {
    return GetModeEField(pos, IWaveguide::kTEM, A, mode.GetModeIndex1(),
                         mode.GetModeIndex2(), omega, state);
  }
}

double rad::IWaveguide::GetCutoffFrequency(WaveguideMode mode) {
  if (mode.GetModeType() == ModeType::kTE) {
    return GetCutoffFrequency(IWaveguide::kTE, mode.GetModeIndex1(),
                              mode.GetModeIndex2());
  } else if (mode.GetModeType() == ModeType::kTM) {
    return GetCutoffFrequency(IWaveguide::kTM, mode.GetModeIndex1(),
                              mode.GetModeIndex2());
  } else {
    return GetCutoffFrequency(IWaveguide::kTEM, mode.GetModeIndex1(),
                              mode.GetModeIndex2());
  }
}

double rad::IWaveguide::GetCutoffWavenumber(WaveguideMode mode) {
  if (mode.GetModeType() == ModeType::kTE) {
    return GetCutoffWavenumber(IWaveguide::kTE, mode.GetModeIndex1(),
                               mode.GetModeIndex2());
  } else if (mode.GetModeType() == ModeType::kTM) {
    return GetCutoffWavenumber(IWaveguide::kTM, mode.GetModeIndex1(),
                               mode.GetModeIndex2());
  } else {
    return GetCutoffWavenumber(IWaveguide::kTEM, mode.GetModeIndex1(),
                               mode.GetModeIndex2());
  }
}

bool rad::IWaveguide::ModePropagates(WaveguideMode mode, double f) {
  const double fc{GetCutoffFrequency(mode)};
  if (fc <= 0) {
    return false;
  } else {
    return f > fc;
  }
}

double rad::IWaveguide::GetEFieldIntegral(WaveguideMode mode, double omega,
                                          double A, int nSurfPnts, bool state) {
  if (mode.GetModeType() == ModeType::kTE) {
    return GetEFieldIntegral(IWaveguide::kTE, mode.GetModeIndex1(),
                             mode.GetModeIndex2(), omega, A, nSurfPnts, state);
  } else if (mode.GetModeType() == ModeType::kTM) {
    return GetEFieldIntegral(IWaveguide::kTM, mode.GetModeIndex1(),
                             mode.GetModeIndex2(), omega, A, nSurfPnts, state);
  } else {
    return GetEFieldIntegral(IWaveguide::kTEM, mode.GetModeIndex1(),
                             mode.GetModeIndex2(), omega, A, nSurfPnts, state);
  }
}
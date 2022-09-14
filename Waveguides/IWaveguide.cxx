/// IWaveguide.cxx

#include "Waveguides/IWaveguide.h"
#include "BasicFunctions/Constants.h"

#include <iostream>

double rad::IWaveguide::GetModeImpedance(Mode_t modeType, unsigned int n, unsigned int m, double omega)
{
  double Z{ -1.0 };
  double k_c{ GetCutoffWavenumber(modeType, n, m) };
  double k{ omega/TMath::C() };
  double beta{ sqrt(k*k - k_c*k_c) };
  
  if (modeType == kTE) {
    Z = k * sqrt(MU0/EPSILON0) / beta;
  }
  else if (modeType == kTM) {
    Z = beta * sqrt(MU0/EPSILON0) / k;
  }
  else {
    std::cout<<"Currently we are not set up to deal with TEM modes. You will not get a result that makes sense."<<std::endl;
  }
  return Z;
}

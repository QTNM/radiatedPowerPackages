/// BasicTritiumSpectrum.cxx

#include "EventGeneration/BasicTritiumSpectrum.h"
#include "BasicFunctions/Constants.h"

#include "TMath.h"
#include "TRandom3.h"

#include <climits>

rad::BasicTritiumSpectrum::BasicTritiumSpectrum(double mass1, double mass2, double mass3, double theta12, double theta13, double theta23, double endpointE)
{
  m1 = mass1;
  m2 = mass2;
  m3 = mass3;
  th12 = theta12;
  th13 = theta13;
  th23 = theta23;
  endpoint = endpointE;
  Ue1Sq = CalculateUe1Sq();
  Ue2Sq = CalculateUe2Sq();
  Ue3Sq = CalculateUe3Sq();

  maxRate = -DBL_MAX;
  int sampledPoints = 50000;
  for (int i = 0; i < 50000; i++) {
    double thisRate = GetDecayRate( double(i)/double(sampledPoints) * endpoint );
    if (thisRate > maxRate) maxRate = thisRate;
  }
}

double rad::BasicTritiumSpectrum::GetDecayRate(double electronEnergy)
{
  double beta = sqrt( 1.0 - pow(ME_EV/(electronEnergy + ME_EV), 2) );
  double p = sqrt( electronEnergy*electronEnergy + 2*ME_EV*electronEnergy );
  double eta = ALPHA * 2 / beta;
  double fermiFunction = 2*TMath::Pi()*eta / ( 1 - exp(-2*TMath::Pi()*eta) );
  double nuE = endpoint - electronEnergy;
  double rate = G_F*G_F * pow(0.97425, 2) * fermiFunction * (1 + 3.0 * pow(-1.2646, 2)) * p*(electronEnergy+ME_EV) / (2*pow(TMath::Pi(), 3));

  double neutrinoPhaseSpc1 = (m1 <= nuE) ? Ue1Sq * nuE * sqrt(nuE*nuE - m1*m1) : 0.0;
  double neutrinoPhaseSpc2 = (m2 <= nuE) ? Ue2Sq * nuE * sqrt(nuE*nuE - m2*m2) : 0.0;
  double neutrinoPhaseSpc3 = (m3 <= nuE) ? Ue3Sq * nuE * sqrt(nuE*nuE - m3*m3) : 0.0;
  double neutrinoPhaseSpc = neutrinoPhaseSpc1 + neutrinoPhaseSpc2 + neutrinoPhaseSpc3;  
  rate *= neutrinoPhaseSpc;
  rate *= 1.0 / 6.58e-16; // Account for natural units
  
  return rate;
}

double rad::BasicTritiumSpectrum::DrawRandomEnergy()
{
  // Use brute force accept-reject method to get the point
  TRandom3* numGen = new TRandom3(0);
  double eng = 0.0;
  for (int iPnt = 0; iPnt < INT_MAX; iPnt++) {
    eng = numGen->Uniform(0, endpoint);
    double randomRate = numGen->Uniform(0, maxRate);
    if (randomRate < GetDecayRate(eng)) break;
  }
  delete numGen;
  
  return eng;
}

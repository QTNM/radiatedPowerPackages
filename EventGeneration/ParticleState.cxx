/// ParticleState.cxx

#include "EventGeneration/ParticleState.h"

#include "TVector3.h"
#include "TMath.h"

rad::ParticleState::ParticleState(double tStart, double mass, double charge, TVector3 initialPos, TVector3 initialVel)
{
  startTime = tStart;
  currentTime = tStart;
  particleMass = mass;
  particleCharge = charge;
  positionVector = initialPos;
  velocityVector = initialVel;
}

double rad::ParticleState::GetE()
{
  double betaSq{pow(velocityVector.Mag() / TMath::C(), 2)};
  return {(1.0 / sqrt(1.0 - betaSq)) * particleMass * pow(TMath::C(), 2)};
}

double rad::ParticleState::GetKE()
{
  return GetE() - particleMass * pow(TMath::C(), 2);
}
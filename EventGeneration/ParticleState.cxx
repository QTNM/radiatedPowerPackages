/// ParticleState.cxx

#include "EventGeneration/ParticleState.h"

#include "TVector3.h"

rad::ParticleState::ParticleState(double tStart, double mass, double charge, TVector3 initialPos, TVector3 initialVel)
{
  startTime   = tStart;
  currentTime = tStart;
  particleMass   = mass;
  particleCharge = charge;
  positionVector = initialPos;
  velocityVector = initialVel;
}

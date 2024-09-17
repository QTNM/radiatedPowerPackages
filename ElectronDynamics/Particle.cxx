#include "ElectronDynamics/Particle.h"

rad::Particle::Particle(TVector3 pos, TVector3 vel, double m, double q) {
  position = pos;
  velocity = vel;
  mass = m;
  charge = q;
}

void rad::Particle::UpdateState(TVector3 newpos, TVector3 newvel,
                                TVector3 newacc) {
  position = newpos;
  velocity = newvel;
  acceleration = newacc;
}
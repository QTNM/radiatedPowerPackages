#include "ElectronDynamics/Particle.h"

unsigned int rad::Particle::counter = 0;

rad::Particle::Particle(TVector3 pos, TVector3 vel, TVector3 acc, double m,
                        double q) {
  position = pos;
  velocity = vel;
  acceleration = acc;
  mass = m;
  charge = q;
  id = ++counter;
}

void rad::Particle::UpdateState(TVector3 newpos, TVector3 newvel,
                                TVector3 newacc) {
  position = newpos;
  velocity = newvel;
  acceleration = newacc;
}
/// Event.cxx

#include "EventGeneration/Event.h"

#include "EventGeneration/ParticleState.h"
#include "ElectronDynamics/BaseField.h"
#include "ElectronDynamics/BorisSolver.h"
#include "Antennas/IAntenna.h"

#include "TMath.h"
#include "TVector3.h"

#include <cassert>
#include <tuple>

rad::Event::Event(std::vector<ParticleState> particles, BaseField *field,
                  double simStepSize, double simTime)
{
  assert(particles.size() > 0);

  particleList = particles;
  magneticField = field;
  simulationStepSize = simStepSize;
  maximumSimulationTime = simTime;
  antennaList = {};

  clockTime = 0.0; // Just for simplicity

  // Create a separate Boris solver for each particle
  // keep them in the same ordering as the particles
  solverList.clear();
  for (auto &part : particleList)
  {
    BorisSolver solver(field, part.GetParticleCharge(), part.GetParticleMass(),
                       2 * R_E / (3 * TMath::C()));
    solverList.push_back(solver);
  }
}

rad::Event::Event(std::vector<ParticleState> particles, std::vector<IAntenna *> antennas,
                  BaseField *field, double simStepSize, double simTime)
{
  assert(particles.size() > 0);

  particleList = particles;
  magneticField = field;
  simulationStepSize = simStepSize;
  maximumSimulationTime = simTime;
  antennaList = antennas;

  clockTime = 0.0; // Just for simplicity

  // Create a separate Boris solver for each particle
  // keep them in the same ordering as the particles
  solverList.clear();
  for (auto &part : particleList)
  {
    BorisSolver solver(field, part.GetParticleCharge(), part.GetParticleMass(),
                       2 * R_E / (3 * TMath::C()));
    solverList.push_back(solver);
  }
}

bool rad::Event::ParticleStartCheck(ParticleState part)
{
  return (part.currentTime <= clockTime);
}

void rad::Event::AdvanceParticleStep(int particleNumber)
{
  double timeStep{clockTime - particleList[particleNumber].currentTime};
  std::tuple<TVector3, TVector3> outputVectors =
      solverList[particleNumber].advance_step(timeStep,
                                              particleList[particleNumber].GetPositionVector(),
                                              particleList[particleNumber].GetVelocityVector());

  // Update the relevant stats
  particleList[particleNumber].positionVector = std::get<0>(outputVectors);
  particleList[particleNumber].velocityVector = std::get<1>(outputVectors);
  particleList[particleNumber].currentTime += simulationStepSize;
}

void rad::Event::PropagateParticles(const char *outputFile, std::vector<OutputVar> vars)
{
  std::cout << "Simulation time is " << maximumSimulationTime << " s" << std::endl;
  std::cout << "Simulation step size is " << simulationStepSize << " s" << std::endl;
  int nTimeSteps = maximumSimulationTime / simulationStepSize;
  std::cout << "Time steps = " << nTimeSteps << std::endl;

  // Move forward through time
  for (int iStep = 0; iStep < nTimeSteps; iStep++)
  {
    clockTime += simulationStepSize;
    for (int iPart = 0; iPart < particleList.size(); iPart++)
    {
      // Check if this particle is meant to be active yet
      if (ParticleStartCheck(particleList[iPart]))
      {
        AdvanceParticleStep(iPart);
      }
    }
  }
}

rad::ParticleState rad::Event::GetParticle(int particleIndex)
{
  if (particleIndex >= 0 && particleIndex < particleList.size())
  {
    return particleList[particleIndex];
  }
  else
  {
    // Return a placeholder state so the application doesn't just crash
    std::cout << "Requested a particle that does not exist!" << std::endl;
    ParticleState placeholder(clockTime, 0.0, 0.0,
                              TVector3(9999.9, 9999.9, 9999.9),
                              TVector3(0.0, 0.0, 0.0));
    return placeholder;
  }
}

double rad::Event::GetPropagationTime(ParticleState particle, IAntenna *antenna)
{
  TVector3 disp = particle.GetPositionVector() - antenna->GetAntennaPosition();
  return (disp.Mag() / TMath::C());
}

unsigned int rad::Event::GetNParticles()
{
  return particleList.size();
}
/// Event.cxx

#include "EventGeneration/Event.h"

#include "EventGeneration/ParticleState.h"
#include "ElectronDynamics/BaseField.h"
#include "ElectronDynamics/BorisSolver.h"
#include "Antennas/IAntenna.h"

#include "TFile.h"
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

  // Check if we have chosen to write output to file
  bool writeToFile{outputFile != NULL};
  TFile *fout = 0;
  TTree *tree = 0;

  // Create the output file
  if (writeToFile)
  {
    std::cout << "We are writing to file.\n";
    fout = new TFile(outputFile, "RECREATE");
    // Create the TTree with the correct variables
    tree = CreateOutputTree(vars);
    AddParticleData(tree, vars);
  }

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
        AddParticleData(tree, vars);
      }
    }
  }

  // If we opened the file, then close it
  if (writeToFile)
  {
    fout->cd();
    tree->Write();
    delete tree;
    fout->Close();
    delete fout;
  }
  else
  {
    delete tree;
    delete fout;
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

TTree *rad::Event::CreateOutputTree(std::vector<OutputVar> vars)
{
  auto *tree = new TTree("output", "output");

  // It's hard to imagine why we would want to create a tree without the time info
  tree->Branch("time", &particleTime, "time/D");

  // Loop through the variables and create the appropriate branches
  for (auto &quantity : vars)
  {
    if (quantity == OutputVar::kPos)
    {
      tree->Branch("pos", pos, "pos[3]/D");
    }
    else if (quantity == OutputVar::kVel)
    {
      tree->Branch("vel", vel, "vel[3]/D");
    }
    else if (quantity == OutputVar::kAcc)
    {
      tree->Branch("acc", acc, "acc[3]/D");
    }
    else if (quantity == OutputVar::kBField)
    {
      tree->Branch("bField", bField, "bField[3]/D");
    }
    else
    {
      std::cout << "Unsupported option. This will not do what you think.\n";
    }
  }
  return tree;
}

void rad::Event::AddParticleData(TTree *outputTree, std::vector<OutputVar> vars)
{
  // Get the particle state so we can extract kinematic info
  particleTime = clockTime;
  ParticleState part{GetParticle(0)};
  TVector3 particlePos{part.GetPositionVector()};
  TVector3 particleVel{part.GetVelocityVector()};

  // Branch addresses should be set already so just add the appropriate data
  for (auto &quantity : vars)
  {
    if (quantity == OutputVar::kPos)
    {
      pos[0] = particlePos.X();
      pos[1] = particlePos.Y();
      pos[2] = particlePos.Z();
    }
    else if (quantity == OutputVar::kVel)
    {
      vel[0] = particleVel.X();
      vel[1] = particleVel.Y();
      vel[2] = particleVel.Z();
    }
    else if (quantity == OutputVar::kAcc)
    {
      TVector3 particleAcc{solverList.at(0).acc(particlePos, particleVel)};
      acc[0] = particleAcc.X();
      acc[1] = particleAcc.Y();
      acc[2] = particleAcc.Z();
    }
    else if (quantity == OutputVar::kBField)
    {
      TVector3 fieldAtPoint{magneticField->evaluate_field_at_point(particlePos)};
      bField[0] = fieldAtPoint.X();
      bField[1] = fieldAtPoint.Y();
      bField[2] = fieldAtPoint.Z();
    }
    else
    {
      std::cout << "Unsupported option. This will not do what you think.\n";
    }
  }

  // Actually fill the tree
  outputTree->Fill();
}
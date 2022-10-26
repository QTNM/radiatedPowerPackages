/// Event.cxx

#include "EventGeneration/Event.h"
#include "EventGeneration/ParticleState.h"
#include "EventGeneration/FieldStorer.h"

#include "ElectronDynamics/BaseField.h"
#include "ElectronDynamics/BorisSolver.h"

#include "Antennas/IAntenna.h"

#include "TFile.h"
#include "TMath.h"
#include "TVector3.h"

#include <cassert>
#include <tuple>
#include <algorithm>

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

  nAntennas = 0;
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

  nAntennas = antennaList.size();
  nArrays = antennaList.size();

  // All these antennas are considered to be separate elements
  antennaArrayMap.resize(nAntennas, 0);
  for (int iAnt{0}; iAnt < nAntennas; iAnt++)
  {
    antennaArrayMap.at(iAnt) = iAnt;
  }
}

rad::Event::Event(std::vector<ParticleState> particles, std::vector<AntennaArray> arrays,
                  BaseField *field, double simStepSize, double simTime)
{
  assert(particles.size() > 0);

  particleList = particles;
  magneticField = field;
  simulationStepSize = simStepSize;
  maximumSimulationTime = simTime;
  arrayList = arrays;

  clockTime = 0.0; // For simplicity

  // Create a separate Boris solver for each particle
  // keep them in the same ordering as the particles
  solverList.clear();
  for (auto &part : particleList)
  {
    BorisSolver solver(field, part.GetParticleCharge(), part.GetParticleMass(),
                       2 * R_E / (3 * TMath::C()));
    solverList.push_back(solver);
  }

  nArrays = arrays.size();
  nAntennas = 0;

  int startAntennaIndex{0};
  for (int iArr{0}; iArr < nArrays; iArr++)
  {
    nAntennas += arrayList.at(iArr).GetNElements();
    // Add all the component antennas to the antenna list
    for (int iElement{0}; iElement < arrayList.at(iArr).GetNElements(); iElement++)
    {
      antennaList.push_back(arrayList.at(iArr).GetAntenna(iElement));
    }

    // Now fill in the map so we know which array these belong to
    for (int ii{startAntennaIndex}; ii < startAntennaIndex + nAntennas; ii++)
    {
      antennaArrayMap.push_back(iArr);
    }
  }

  std::cout << "We have a total of " << nArrays << " arrays, made up of " << nAntennas << " elements\n";
  std::cout << "Antenna list has size of " << antennaList.size() << std::endl;
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
  for (int i{0}; i < antennaArrayMap.size(); i++)
  {
    std::cout << "Antenna " << i << " is part of array " << antennaArrayMap.at(i) << std::endl;
  }

  std::cout << "Simulation time is " << maximumSimulationTime << " s" << std::endl;
  std::cout << "Simulation step size is " << simulationStepSize << " s" << std::endl;
  const int nTimeSteps = maximumSimulationTime / simulationStepSize;
  std::cout << "Time steps = " << nTimeSteps << std::endl;

  // Set 2D arrays to defined nonsense values
  ResetArrays();

  // Do we need to compute fields or voltages?
  bool computeFields{(VectorContainsVar(vars, kAntEField) ||
                      VectorContainsVar(vars, kAntBField) ||
                      VectorContainsVar(vars, kAntVoltage))};
  if (computeFields)
  {
    std::cout << "We are computing EM fields\n";
  }

  // Check if we have chosen to write output to file
  bool writeToFile{outputFile != NULL};
  TFile *fout = 0;
  TTree *tree = 0;

  TVector3 eFieldInitial[nAntennas];
  TVector3 bFieldInitial[nAntennas];
  double tAInitial[nAntennas];
  // Create the output file
  if (writeToFile)
  {
    std::cout << "We are writing to file.\n";
    fout = new TFile(outputFile, "RECREATE");
    // Create the TTree with the correct variables
    tree = CreateOutputTree(vars);
    AddLocalParticleData(tree, vars);

    if (antennaList.size() > 0 && computeFields)
    {
      // If we are writing fields, get fields for each antenna
      for (int iAnt{0}; iAnt < nAntennas; iAnt++)
      {
        tAInitial[iAnt] = GetPropagationTime(particleList.at(0),
                                             antennaList.at(iAnt));
        bFieldInitial[iAnt] = GetBFieldAtAntenna(0, iAnt);
        antBField[iAnt][0] = 0.0;
        antBField[iAnt][1] = 0.0;
        antBField[iAnt][2] = 0.0;
        eFieldInitial[iAnt] = GetEFieldAtAntenna(0, iAnt);
        antEField[iAnt][0] = 0.0;
        antEField[iAnt][1] = 0.0;
        antEField[iAnt][2] = 0.0;
      }

      for (int iArr{0}; iArr < nArrays; iArr++)
      {
        antVoltage[iArr] = 0.0;
      }
    }

    tree->Fill();
  }

  // Create an instance of FieldStorer for every antenna
  std::vector<FieldStorer> fsVec;
  if (computeFields)
  {
    // Reserve the appropriate amount of memory
    fsVec.reserve(nAntennas);

    TVector3 posInitial{particleList.at(0).GetPositionVector()};
    // Create the FieldStorers
    for (int iAnt{0}; iAnt < nAntennas; iAnt++)
    {
      FieldStorer fieldStorage(eFieldInitial[iAnt], bFieldInitial[iAnt], posInitial,
                               tAInitial[iAnt], antennaList.at(iAnt));
      fsVec.push_back(fieldStorage);
    }
  }

  // Move forward through time
  for (int iStep = 1; iStep < nTimeSteps; iStep++)
  {
    clockTime += simulationStepSize;

    // Reset to 0 so we can start adding voltages
    for (int iArr{0}; iArr < nArrays; iArr++)
    {
      antVoltage[iArr] = 0.0;
    }

    // Create a sum for each of the arrays, initially all elements zero
    // This will store the element sums of each of the arrays
    int whichArray{0};
    // std::vector<double> arrayVoltageSums(nArrays, 0.0);

    for (int iPart = 0; iPart < particleList.size(); iPart++)
    {
      // Check if this particle is meant to be active yet
      if (ParticleStartCheck(particleList[iPart]))
      {
        AdvanceParticleStep(iPart);
        AddLocalParticleData(tree, vars);

        // Check if we are saving antenna data (and have antenna data to save)
        if (antennaList.size() > 0 && computeFields)
        {
          // Loop over the antennas again and calculate the fields
          for (int iAnt{0}; iAnt < nAntennas; iAnt++)
          {
            double tA{clockTime + GetPropagationTime(particleList.at(0),
                                                     antennaList.at(iAnt))};
            TVector3 eField{GetEFieldAtAntenna(0, iAnt)};
            TVector3 bField{GetBFieldAtAntenna(0, iAnt)};
            TVector3 pos{particleList[iPart].GetPositionVector()};
            fsVec.at(iAnt).AddNewFields(eField, bField, pos, tA);

            // Write to file
            if (VectorContainsVar(vars, kAntEField))
            {
              antEField[iAnt][0] = fsVec.at(iAnt).GetInterpolatedEField(clockTime).X();
              antEField[iAnt][1] = fsVec.at(iAnt).GetInterpolatedEField(clockTime).Y();
              antEField[iAnt][2] = fsVec.at(iAnt).GetInterpolatedEField(clockTime).Z();
            }

            if (VectorContainsVar(vars, kAntBField))
            {
              antBField[iAnt][0] = fsVec.at(iAnt).GetInterpolatedBField(clockTime).X();
              antBField[iAnt][1] = fsVec.at(iAnt).GetInterpolatedBField(clockTime).Y();
              antBField[iAnt][2] = fsVec.at(iAnt).GetInterpolatedBField(clockTime).Z();
            }

            if (VectorContainsVar(vars, kAntVoltage))
            {
              antVoltage[antennaArrayMap.at(iAnt)] += fsVec.at(iAnt).GetAntennaLoadVoltage(clockTime);
            }
          }
        }
      }
    }
    tree->Fill();
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

  // If necessary, write the number of antennas
  if (VectorContainsVar(vars, kAntBField) || VectorContainsVar(vars, kAntEField) ||
      VectorContainsVar(vars, kAntVoltage))
  {
    tree->Branch("nAntennas", &nAntennas, "nAntennas/I");
    tree->Branch("nArrays", &nArrays, "nArrays/I");
  }

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
    else if (quantity == OutputVar::kAntEField)
    {
      tree->Branch("antEField", antEField, "antEField[nAntennas][3]/D");
    }
    else if (quantity == OutputVar::kAntBField)
    {
      tree->Branch("antBField", antBField, "antBField[nAntennas][3]/D");
    }
    else if (quantity == OutputVar::kAntVoltage)
    {
      tree->Branch("antVoltage", antVoltage, "antVoltage[nArrays]/D");
    }
    else
    {
      std::cout << "Unsupported option. This will not do what you think.\n";
    }
  }
  return tree;
}

void rad::Event::AddLocalParticleData(TTree *outputTree, std::vector<OutputVar> vars)
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
  }
}

TVector3 rad::Event::GetEFieldAtAntenna(unsigned int particleIndex, IAntenna *ant)
{
  TVector3 antennaPos{ant->GetAntennaPosition()};
  TVector3 particlePos{particleList.at(particleIndex).GetPositionVector()};
  TVector3 particleVel{particleList.at(particleIndex).GetVelocityVector()};
  double q{particleList.at(particleIndex).GetParticleCharge()};
  TVector3 particleAcc{solverList.at(0).acc(particlePos, particleVel)};
  TVector3 beta{particleVel * (1.0 / TMath::C())};
  TVector3 betaDot{particleAcc * (1.0 / TMath::C())};
  double premult{q / (4 * TMath::Pi() * EPSILON0)};
  double r{(antennaPos - particlePos).Mag()};
  TVector3 rHat{(antennaPos - particlePos).Unit()};

  TVector3 term1{(rHat - beta) * (1 - beta.Dot(beta)) *
                 (1.0 / (pow(1 - beta.Dot(rHat), 3) * r * r))};
  TVector3 term2{rHat.Cross((rHat - beta).Cross(betaDot)) *
                 (1.0 / (TMath::C() * r * pow(1 - beta.Dot(rHat), 3)))};
  return (term1 + term2) * premult;
}

TVector3 rad::Event::GetBFieldAtAntenna(unsigned int particleIndex, IAntenna *ant)
{
  TVector3 eField{GetEFieldAtAntenna(particleIndex, ant)};
  TVector3 antennaPos{ant->GetAntennaPosition()};
  TVector3 particlePos{particleList.at(particleIndex).GetPositionVector()};
  TVector3 rHat{(antennaPos - particlePos).Unit()};
  return rHat.Cross(eField) * (1.0 / TMath::C());
}

TVector3 rad::Event::GetEFieldAtAntenna(unsigned int particleIndex,
                                        unsigned int antennaIndex)
{
  return GetEFieldAtAntenna(particleIndex, antennaList.at(antennaIndex));
}

TVector3 rad::Event::GetBFieldAtAntenna(unsigned int particleIndex,
                                        unsigned int antennaIndex)
{
  return GetBFieldAtAntenna(particleIndex, antennaList.at(antennaIndex));
}

bool rad::Event::VectorContainsVar(std::vector<OutputVar> vars, OutputVar testVar)
{
  return std::find(vars.begin(), vars.end(), OutputVar::kAntEField) != vars.end();
}

void rad::Event::ResetArrays()
{
  for (int iAnt = 0; iAnt < nMaxAntennas; iAnt++)
  {
    for (int iDir = 0; iDir < 3; iDir++)
    {
      antBField[iAnt][iDir] = -9999;
      antEField[iAnt][iDir] = -9999;
    }
    antVoltage[iAnt] = -9999;
  }
}
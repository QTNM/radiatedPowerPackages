/*
    DemonstrateEventGeneration.cxx
    Seb Jones 12-10-22

    Demonstrate the functionality of the event generation branch
*/

#include "BasicFunctions/Constants.h"
#include "BasicFunctions/BasicFunctions.h"

#include "EventGeneration/Event.h"
#include "EventGeneration/ParticleState.h"

#include "ElectronDynamics/QTNMFields.h"

#include "Antennas/HalfWaveDipole.h"
#include "Antennas/AntennaArray.h"

#include "TMath.h"
#include "TVector3.h"

#include <iostream>

using namespace rad;
using std::cout;

int main(int argc, char *argv[])
{
    TVector3 initialVel{0.26 * TMath::C(), 0, 0};
    TVector3 initialPos{0, 0, 0};

    ParticleState testParticle(0.0, ME, -TMath::Qe(), initialPos, initialVel);
    const double initialKE{testParticle.GetKE() / TMath::Qe()}; // eV
    cout << "Initial particle energy is " << testParticle.GetE() / TMath::Qe() << " eV\n";
    cout << "Initial particle KE is " << initialKE / 1e3 << " keV\n";

    const double selectedField{1.0};
    UniformField *field = new UniformField(selectedField);

    TVector3 antennaPos1(0.05, 0, 0);
    TVector3 antennaXAxis1(1, 0, 0);
    TVector3 antennaZAxis1(0, 1, 0);
    const double centralFreq{27.01e9};
    HalfWaveDipole *testPoint1 = new HalfWaveDipole(antennaPos1,
                                                    antennaXAxis1, antennaZAxis1,
                                                    centralFreq);
    HalfWaveDipole *testPoint1Delay = new HalfWaveDipole(antennaPos1,
                                                         antennaXAxis1, antennaZAxis1,
                                                         centralFreq, 1e-10);

    std::cout << "Delay 1, 2 = " << testPoint1->GetTimeDelay() << ", " << testPoint1Delay->GetTimeDelay() << std::endl;

    TVector3 antennaPos2(-0.05, 0, 0);
    TVector3 antennaXAxis2(-1, 0, 0);
    TVector3 antennaZAxis2(0, -1, 0);
    HalfWaveDipole *testPoint2 = new HalfWaveDipole(antennaPos2,
                                                    antennaXAxis2, antennaZAxis2,
                                                    centralFreq);
    std::vector<IAntenna *> antennaVec{testPoint1, testPoint1Delay};
    AntennaArray twoElementArray(antennaVec);
    std::cout<<"Our array has "<<twoElementArray.GetNElements()<<" elements.\n";
    std::vector<AntennaArray> twoElementArrayVec{twoElementArray};
    // Create an event containing just one particle
    const double simStepSize{1e-12};
    const double simTime{4e-10};
    Event ev({testParticle}, twoElementArrayVec, field, simStepSize, simTime);
    const double initialClockTime{ev.GetClockTime()};

    // Now propagate the particles in time
    const char *filePath{"~/work/qtnm/outputs/DemonstrateEventGeneration/output.root"};
    std::vector<OutputVar> vars{
        OutputVar::kPos,
        OutputVar::kVel,
        OutputVar::kAcc,
        OutputVar::kBField,
        OutputVar::kAntEField,
        OutputVar::kAntBField,
        OutputVar::kAntVoltage};

    ev.PropagateParticles(filePath, vars);
    const double finalKE{ev.GetParticle(0).GetKE() / TMath::Qe()}; // eV
    const double finalClockTime{ev.GetClockTime()};

    // And get the final particle energy
    cout << "There are " << ev.GetNParticles() << " particles in the event.\n";
    cout << "Final particle KE is " << finalKE / 1e3 << " keV\n";

    cout << "Final position vector is (" << ev.GetParticle(0).GetPositionVector().X() << ", " << ev.GetParticle(0).GetPositionVector().Y() << ", " << ev.GetParticle(0).GetPositionVector().Z() << ")\n";

    const double radiatedPower{((initialKE - finalKE) * TMath::Qe()) / simTime};
    cout << "Radiated power is " << radiatedPower * 1e15 << " fW\n";

    return 0;
}
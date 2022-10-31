/// SignalStorer.cxx

#include "EventGeneration/SignalStorer.h"

rad::SignalStorer::SignalStorer(double vInitial, double tInitial,
                                double sRate, LocalOscillator osc)
{
    sampleRate = sRate;
    maxTimeLength = 200 / (sRate); // 200 samples
    lo = osc;

    times.push_back(tInitial);
    vDownmixedI.push_back(vInitial * lo.GetInPhaseComponent(tInitial));
    vDownmixedQ.push_back(vInitial * lo.GetQuadratureComponent(tInitial));
}


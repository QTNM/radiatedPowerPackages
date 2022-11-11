/// SignalStorer.cxx

#include "EventGeneration/SignalStorer.h"
#include "BasicFunctions/BasicFunctions.h"

rad::SignalStorer::SignalStorer(double vInitial, double tInitial,
                                double sRate, LocalOscillator osc)
{
  sampleRate = sRate;
  maxTimeLength = 100 / (sRate); // 200 samples
  // maxTimeLength = 50 / (sRate); // 200 samples
  lo = osc;
  nextSampleTime = tInitial + 1.0 / sRate;
  nextSample10Time = tInitial + 1.0 / (10 * sRate);

  times.push_back(tInitial);
  vDownmixedI.push_back(vInitial * lo.GetInPhaseComponent(tInitial));
  vDownmixedQ.push_back(vInitial * lo.GetQuadratureComponent(tInitial));
  vSample10I.push_back(vInitial * lo.GetInPhaseComponent(tInitial));
  vSample10Q.push_back(vInitial * lo.GetQuadratureComponent(tInitial));
  timesSample10.push_back(tInitial);
}

void rad::SignalStorer::AddNewPoint(double v, double t)
{
  vDownmixedI.push_back(v * lo.GetInPhaseComponent(t));
  vDownmixedQ.push_back(v * lo.GetQuadratureComponent(t));
  times.push_back(t);
  // Delete points from the beginning if we're over the max time
  if (times.at(times.size() - 1) - times.at(0) > maxTimeLength)
  {
    vDownmixedI.erase(vDownmixedI.begin());
    vDownmixedQ.erase(vDownmixedQ.begin());
    times.erase(times.begin());
  }

  // If we are at or beyond a sample time, then sample the values
  if (t >= nextSample10Time)
  {
    // Use a cubic interpolation. Get the relevant values first
    unsigned int nVals = times.size();
    std::vector<double> timeVals{times.at(nVals - 4), times.at(nVals - 3),
                                 times.at(nVals - 2), times.at(nVals - 1)};
    std::vector<double> viVals{vDownmixedI.at(nVals - 4), vDownmixedI.at(nVals - 3),
                               vDownmixedI.at(nVals - 2), vDownmixedI.at(nVals - 1)};
    std::vector<double> vqVals{vDownmixedQ.at(nVals - 4), vDownmixedQ.at(nVals - 3),
                               vDownmixedQ.at(nVals - 2), vDownmixedQ.at(nVals - 1)};
    double vIInterp{CubicInterpolation(timeVals, viVals, nextSample10Time)};
    double vQInterp{CubicInterpolation(timeVals, vqVals, nextSample10Time)};
    vSample10I.push_back(vIInterp);
    vSample10Q.push_back(vQInterp);

    timesSample10.push_back(nextSample10Time);
    nextSample10Time += 1.0 / (10 * sampleRate);

    if (timesSample10.at(timesSample10.size() - 1) - timesSample10.at(0) > maxTimeLength)
    {
      timesSample10.erase(timesSample10.begin());
      vSample10I.erase(vSample10I.begin());
      vSample10Q.erase(vSample10Q.begin());
    }
  }

  // Now check if we are past the next coarse sample time
  // If we are, we want to filter the downmixed and sampled data
  if (t >= nextSampleTime)
  {
    // Clear these ready to have new data added
    vFilterI.clear();
    vFilterQ.clear();
    // Apply a band pass filter
    vFilterI = BandPassFilter(timesSample10, vSample10I, 0.0, sampleRate / 2.0);
    vFilterQ = BandPassFilter(timesSample10, vSample10Q, 0.0, sampleRate / 2.0);

    // Advance time to the next sample rate
    nextSampleTime += 1.0 / sampleRate;
  }
}

double rad::SignalStorer::GetVI(double time)
{
  // We may not actually need to do any interpolation
  if (time == timesSample10.at(timesSample10.size() - 1))
  {
    return vFilterI.at(vFilterI.size() - 1);
  }
  else
  {
    unsigned int nVals =timesSample10.size();
    std::vector<double> timeVals{timesSample10.at(nVals - 4), timesSample10.at(nVals - 3),
                                 timesSample10.at(nVals - 2), timesSample10.at(nVals - 1)};
    std::vector<double> viVals{vFilterI.at(nVals - 4), vFilterI.at(nVals - 3),
                               vFilterI.at(nVals - 2), vFilterI.at(nVals - 1)};
    return CubicInterpolation(timeVals, viVals, time);
  }
}

double rad::SignalStorer::GetVQ(double time)
{
  // We may not actually need to do any interpolation
  if (time == timesSample10.at(timesSample10.size() - 1))
  {
    return vFilterQ.at(vFilterQ.size() - 1);
  }
  else
  {
    unsigned int nVals = timesSample10.size();
    std::vector<double> timeVals{timesSample10.at(nVals - 4), timesSample10.at(nVals - 3),
                                 timesSample10.at(nVals - 2), timesSample10.at(nVals - 1)};
    std::vector<double> vqVals{vFilterQ.at(nVals - 4), vFilterQ.at(nVals - 3),
                               vFilterQ.at(nVals - 2), vFilterQ.at(nVals - 1)};
    return CubicInterpolation(timeVals, vqVals, time);
  }
}
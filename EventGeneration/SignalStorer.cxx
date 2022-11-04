/// SignalStorer.cxx

#include "EventGeneration/SignalStorer.h"
#include "BasicFunctions/BasicFunctions.h"

rad::SignalStorer::SignalStorer(double vInitial, double tInitial,
                                double sRate, LocalOscillator osc)
{
  sampleRate = sRate;
  maxTimeLength = 100 / (sRate); // 200 samples
  lo = osc;
  nextSampleTime = tInitial + 1.0 / sRate;
  nextSample10Time = tInitial + 1.0 / (10 * sRate);

  times.push_back(tInitial);
  vDownmixedI.push_back(vInitial * lo.GetInPhaseComponent(tInitial));
  vDownmixedQ.push_back(vInitial * lo.GetQuadratureComponent(tInitial));
  timesSample10.push_back(tInitial);
}

void rad::SignalStorer::AddNewPoint(double v, double t)
{
  vDownmixedI.push_back(v * lo.GetInPhaseComponent(t));
  vDownmixedQ.push_back(v * lo.GetQuadratureComponent(t));
  times.push_back(t);
  // Delete points from the beginning if we're over the max time
  if (vDownmixedI.at(vDownmixedI.size() - 1) - vDownmixedI.at(0) > maxTimeLength)
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
    double vIInterp{DoCubicInterpolation(timeVals, viVals, nextSample10Time)};
    double vQInterp{DoCubicInterpolation(timeVals, vqVals, nextSample10Time)};
    vSample10I.push_back(vIInterp);
    vSample10Q.push_back(vQInterp);

    timesSample10.push_back(nextSample10Time);
    nextSample10Time += 1.0 / (10 * sampleRate);

    if (vSample10I.at(vSample10I.size() - 1) - vSample10I.at(0) > maxTimeLength)
    {
      timesSample10.erase(timesSample10.begin());
      vSample10I.erase(vSample10I.begin());
      vSample10Q.erase(vSample10Q.begin());
    }
  }

  // Now check if we are past the next coarse sample time
  // If we are, we want to filter the downmixed data
  if (t >= nextSampleTime)
  {
    // Clear these ready to have new data inserted
    vFilterI.clear();
    vFilterQ.clear();

    nextSampleTime += 1.0 / sampleRate;
  }
}

double rad::SignalStorer::DoCubicInterpolation(std::vector<double> xVals,
                                               std::vector<double> yVals, double interp)
{
  // Use the 3rd integrating Lagrange polynomial
  // First calculate the Lagrange interpolating basis functions
  double l0{(interp - xVals[1]) * (interp - xVals[2]) * (interp - xVals[3]) /
            ((xVals[0] - xVals[1]) * (xVals[0] - xVals[2]) * (xVals[0] - xVals[3]))};
  double l1{(interp - xVals[0]) * (interp - xVals[2]) * (interp - xVals[3]) /
            ((xVals[1] - xVals[0]) * (xVals[1] - xVals[2]) * (xVals[1] - xVals[3]))};
  double l2{(interp - xVals[0]) * (interp - xVals[1]) * (interp - xVals[3]) /
            ((xVals[2] - xVals[0]) * (xVals[2] - xVals[1]) * (xVals[2] - xVals[3]))};
  double l3{(interp - xVals[0]) * (interp - xVals[1]) * (interp - xVals[2]) /
            ((xVals[3] - xVals[0]) * (xVals[3] - xVals[1]) * (xVals[3] - xVals[2]))};

  return l0 * yVals[0] + l1 * yVals[1] + l2 * yVals[2] + l3 * yVals[3];
}

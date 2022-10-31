/*
    SignalStorer.h

    Class for storing time series of signal processing quantities

    S. Jones 31/10/2022
*/

#ifndef SIGNAL_STORER_H
#define SIGNAL_STORER_H

#include "SignalProcessing/LocalOscillator.h"

#include <vector>

namespace rad
{
    class SignalStorer
    {
    private:
        // Times
        std::vector<double> times;
        // Downmixed voltage values
        std::vector<double> vDownmixedI;
        std::vector<double> vDownmixedQ;
        // Values after first downsampling
        std::vector<double> vSample1;
        // Values after filtering
        std::vector<double> vFilter;
        // Values after second downsampling
        std::vector<double> vSample2;

        // The desired final sample rate
        double sampleRate;
        // The maximum time length we wish to store
        double maxTimeLength;

        // The local oscillator used in the signal processing
        LocalOscillator lo;

    public:
        SignalStorer(double vInitial, double tInitial,
                     double sRate, LocalOscillator osc);
    };
}

#endif
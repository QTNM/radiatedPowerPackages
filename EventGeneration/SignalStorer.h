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
		std::vector<double> timesSample10;
		std::vector<double> vSample10I;
		std::vector<double> vSample10Q;
		// Values after filtering
		std::vector<double> vFilterI;
		std::vector<double> vFilterQ;
		// Values after second downsampling
		// std::vector<double> vSample2I;
		// std::vector<double> vSample2Q;

		// The desired final sample rate
		double sampleRate;
		// The maximum time length we wish to store
		double maxTimeLength;

		double nextSampleTime;
		double nextSample10Time;

		// The local oscillator used in the signal processing
		LocalOscillator lo;

	public:
		SignalStorer(double vInitial, double tInitial,
								 double sRate, LocalOscillator osc);

		/// Adds a new time point to the system
		/// \param v The new voltage to be added (in volts)
		/// \param t The time at which the voltage is measured (in seconds)
		void AddNewPoint(double v, double t);
	};
}

#endif
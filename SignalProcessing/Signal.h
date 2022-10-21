/*
  Signal.h
  Class for a signal detected on an antenna or antennas
*/

#ifndef SIGNAL_H
#define SIGNAL_H

#include "SignalProcessing/InducedVoltage.h"
#include "SignalProcessing/LocalOscillator.h"
#include "SignalProcessing/NoiseFunc.h"
#include "FieldClasses/FieldClasses.h"
#include "BasicFunctions/BasicFunctions.h"

#include <vector>

#include "TGraph.h"
#include "TH2.h"

namespace rad
{
  /// Signal base class
  /// Class members are:
  /// grVITime In phase processed voltage component
  /// grVQTime Quadrature processed voltage component
  /// sampleRate Hypothetical sample rate of a digitiser
  
  class Signal
  {
  private:
    /// Function to actually do the signal processing for a given time chunk
    void ProcessTimeChunk(InducedVoltage iv, LocalOscillator lo, double thisChunk, double lastChunk,
			  std::vector<GaussianNoise> noiseTerms,
			  double &firstSampleTime, double &firstSample10Time, bool firstVoltage=true);
    
  protected:
    TGraph* grVITime; // In phase component
    TGraph* grVQTime; // Quadrature component
    double sampleRate;

    /// Function for performing downsampling on a time domain wave using the class sample rate
    /// \param grInput The input time domain waveform
    /// Returns the downsampled waveform    
    TGraph* SampleWaveform(TGraph* grInput);

    /// Function for performing downsampling on a time domain wave
    /// \param grInput The input time domain waveform
    /// \param sRate The designated sample rate
    /// Returns the downsampled waveform
    TGraph* SampleWaveform(TGraph* grInput, const double sRate);

    /// Function for performing downsampling on a time domain wave
    /// Allows for the specification of the first sample time
    /// \param grInput The input time domain waveform
    /// \param sRate The designated sample rate
    /// \param firstSampleTime The first time to sample
    /// Returns the downsampled waveform
    TGraph* SampleWaveform(TGraph* grInput, const double sRate, const double firstSampleTime);

    /// Adds Gaussian white noise to a time series of voltage values
    /// \param grInput The input TGraph
    /// \param noiseTerms Vector containing the various noise contributions
    /// \param IsComponent Are we adding noise to in phase and voltage components?
    void AddGaussianNoise(TGraph* grInput, std::vector<GaussianNoise> noiseTerms,
			  bool IsComponent=true);

    /// Default constructor
    Signal();
    
  public:
    ~Signal();
    
    /// \param fp The field points where the antenna is located
    /// \param lo The local oscillator used to do the downmixing
    /// \param srate The sample rate of a hypothetical ADC
    /// \param noiseTerms A vector of noise terms to be added to the signal
    /// \param kUseRetardedTime Boolean on whether or not to use the retarded time
    Signal(std::vector<FieldPoint> fp, LocalOscillator lo, double srate, std::vector<GaussianNoise> noiseTerms={}, const bool kUseRetardedTime=false);

    /// Constructor for a single field point
    /// \param fp The field point where the antenna is located
    /// \param lo The local oscillator used to do the downmixing
    /// \param srate The sample rate of a hypothetical ADC
    /// \param noiseTerms A vector of noise terms to be added to the signal
    /// \param kUseRetardedTime Boolean on whether or not to use the retarded time
    Signal(FieldPoint fp, LocalOscillator lo, double srate, std::vector<GaussianNoise> noiseTerms={}, const bool kUseRetardedTime=false);

    /// Constructor for a single induced voltage
    /// \param iv The induced voltage at a particular field point
    /// \param lo The local oscillator used to do the downmixing
    /// \param srate The sample rate of a hypothetical ADC
    /// \param noiseTerms A vector of noise terms to be added to the signal
    Signal(InducedVoltage iv, LocalOscillator lo, double srate,
	   std::vector<GaussianNoise> noiseTerms={}, double maxTime=-1);

    /// Constructor for multiple induced voltages
    /// \param iv The vector of induced voltages
    /// \param lo The local oscillator used to do the downmixing
    /// \param srate The sample rate of a hypothetical ADC
    /// \param noiseTerms A vector of noise terms to be added to the signal
    Signal(std::vector<InducedVoltage> iv, LocalOscillator lo, double srate,
	   std::vector<GaussianNoise> noiseTerms={}, double maxTime=-1);
    
    /// Copy constructor
    Signal(const Signal &s1);

    /// Returns the summed in-phase voltage component after all signal processing
    TGraph* GetVITimeDomain(int firstPoint=-1, int lastPoint=-1);
    /// Returns the summed quadrature voltage component after all signal processing
    TGraph* GetVQTimeDomain(int firstPoint=-1, int lastPoint=-1);

    /// Returns the power spectrum of the in-phase voltage component after all signal processing
    /// \param loadResistance The load resistance used for the power calculation
    /// \param firstPoint The first point in the time series to use
    /// \param lastPoint The last point in the time series to use
    TGraph* GetVIPowerNorm(const double loadResistance, int firstPoint=-1, int lastPoint=-1);

    /// Returns the power spectrum of the quadrature voltage component after all signal processing
    /// \param loadResistance The load resistance used for the power calculation
    /// \param firstPoint The first point in the time series to use
    /// \param lastPoint The last point in the time series to use
    TGraph* GetVQPowerNorm(const double loadResistance, int firstPoint=-1, int lastPoint=-1);

    /// Returns the power spectrum of the in-phase voltage component after all signal processing
    /// Power spectrum in this case is the periodogram
    /// \param loadResistance The load resistance used for the power calculation
    /// \param firstPoint The first point in the time series to use
    /// \param lastPoint The last point in the time series to use
    TGraph* GetVIPowerPeriodogram(const double loadResistance, int firstPoint=-1, int lastPoint=-1);

    /// Returns the power spectrum of the quadrature voltage component after all signal processing
    /// Power spectrum in this case is the periodogram
    /// \param loadResistance The load resistance used for the power calculation
    /// \param firstPoint The first point in the time series to use
    /// \param lastPoint The last point in the time series to use
    TGraph* GetVQPowerPeriodogram(const double loadResistance, int firstPoint=-1, int lastPoint=-1);

    /// Returns the 2D spectrogram made using the in-phase voltage component
    /// \param loadResistance The load resistance used for the power calculation
    /// \param NSamplesPerTimeBin The number of time samples used to make each time bin
    TH2D* GetVISpectrogram(const double loadResistance, const int NSamplesPerTimeBin);

    /// Returns the 2D spectrogram made using the quadrature voltage component
    /// \param loadResistance The load resistance used for the power calculation
    /// \param NSamplesPerTimeBin The number of time samples used to make each time bin
    TH2D* GetVQSpectrogram(const double loadResistance, const int NSamplesPerTimeBin);

    /// Returns a 2D sparse spectrogram made using the in-phase voltage component
    /// along with a threshold for considering a pixel as a hit
    /// \param loadResistance The load resistance used for power calculation
    /// \param NSamplesPerTimeBin The number of time samples used to make each time bin
    /// \param ThresholdPower The power (in Watts) above which a pixel is considered a hit
    TH2D* GetVISparseSpectrogram(const double loadResistance, const int NSamplesPerTimeBin, const double ThresholdPower);

    /// Returns a 2D sparse spectrogram made using the quadrature voltage component
    /// along with a threshold for considering a pixel as a hit
    /// \param loadResistance The load resistance used for power calculation
    /// \param NSamplesPerTimeBin The number of time samples used to make each time bin
    /// \param ThresholdPower The power (in Watts) above which a pixel is considered a hit
    TH2D* GetVQSparseSpectrogram(const double loadResistance, const int NSamplesPerTimeBin, const double ThresholdPower);
    
    /// Returns the 2D spectrogram made using the in-phase voltage component
    /// \param loadResistance The load resistance used for the power calculation
    /// \param NSamplesPerTimeBin The number of time samples used to make each time bin
    TH2D* GetVISpectrogramNorm(const double loadResistance, const int NSamplesPerTimeBin);

    /// Returns the 2D spectrogram made using the quadrature voltage component
    /// \param loadResistance The load resistance used for the power calculation
    /// \param NSamplesPerTimeBin The number of time samples used to make each time bin
    TH2D* GetVQSpectrogramNorm(const double loadResistance, const int NSamplesPerTimeBin);

    /// Downmixes the inputted signal with given local oscillator
    /// \param grInput The input time domain signal
    /// \param lo The local oscillator used to do the downmixing
    /// \return The in-phase voltage component in the time domain
    TGraph* DownmixInPhase(TGraph* grInput, LocalOscillator lo);

    /// Downmixes the inputted signal with given local oscillator
    /// \param grInput The input time domain signal
    /// \param lo The local oscillator used to do the downmixing
    /// \return The quadrature voltage component in the time domain
    TGraph* DownmixQuadrature(TGraph* grInput, LocalOscillator lo);

    /// Returns a signal with the dechirping operator applied
    /// \param alpha Parameter defining the strength of the chirp (units of s^-2)
    /// \param firstPoint The first point in the time series to return
    /// \param lastPoint The last point in the time series to return
    /// \return The voltage signal after the operator has been applied
    TGraph* GetDechirpedSignalTimeDomain(const double alpha, int firstPoint=-1, int lastPoint=-1);

    /// Returns the 2D spectrogram made using the dechirped voltage signal
    /// \param loadResistance The load resistance used for the power calculation
    /// \param NSamplesPerTimeBin The number of time samples used to make each time bin
    /// \param alpha Parameter defining the strength of the dechirp (units of s^-2)
    /// \return The spectrogram with the dechirping operator applied
    TH2D* GetDechirpedSpectrogram(const double loadResistance, const int NSamplesPerTimeBin, const double alpha);
  };
}

#endif

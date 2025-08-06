/*
  InducedVoltage.h

  Lighter-weight representation of the voltage induced on an antenna compared to FieldPoint
  Just contains the actual output voltage across the terminals of the chosen antenna
*/

#ifndef INDUCED_VOLTAGE_H
#define INDUCED_VOLTAGE_H

#include "physics/FieldClasses/FieldClasses.h"
#include "physics/Antennas/IAntenna.h"

#include "TGraph.h"
#include "TString.h"

#include <vector>

namespace rad
{
  class InducedVoltage
  {
  private:
    // Open circuit voltage with no signal processing performed
    TString theFile;
    TGraph* grVoltage;
    std::vector<IAntenna*> theAntennas;
    bool UseRetardedTime;
    double chunkSize;

    /// Produce a shifted version of a voltage graph according to the delay of the antenna
    /// \param grIn The graph to be delayed
    /// \param ant The antenna which is delayed a certain amount
    /// \return The time series voltage graph of the delayed signal
    TGraph* DelayVoltage(TGraph* grIn, IAntenna* ant);

    /// Function for processing time chunks of the voltage
    /// \param fp The FieldPoint for which to generate the voltage
    /// \param firstTime The first time from which to generate from
    /// \param lastTime The last time from which to generate from
    /// \param minTime The first time in the whole file to generate from
    /// \param latestStartTime Variable keeping track of the points to cut following processing
    /// \param kFirstAntenna Boolean for keeping track of whether this is the first antenna in a series
    void ProcessTimeChunk(FieldPoint fp, double firstTime, double lastTime, double minTime, double &latestStartTime, bool kFirstAntenna);
    
  public:
    /// Constructor for a voltage
    /// \param trajectoryFilePath Path to the file containing the electron trajectory
    /// \param myAntenna Pointer to the chosen antenna
    /// \param kUseRetardedTime Boolean to select the use of the retarded time
    InducedVoltage(TString trajectoryFilePath, IAntenna* myAntenna,
		   const bool kUseRetardedTime=false);

    /// Constructor taking a vector of antennas as an argument
    /// \param trajectoryFilePath Path to the file containing the electron trajectory
    /// \param antennaVec Vector of antenna points to generate the voltages at
    /// \param kUseRetardedTime Boolean to select the use of the retarded time
    InducedVoltage(TString trajectoryFilePath, std::vector<IAntenna*> antennaVec,
		   const bool kUseRetardedTime=false);

    
    /// Generates the voltage between the specified times
    /// \param minTime First time point to generate voltage from
    /// \param maxTime Maximum time period to generate voltage for
    void GenerateVoltage(double minTime=-1, double maxTime=-1);
    
    /// Destructor
    ~InducedVoltage();

    /// Copy constructor
    InducedVoltage(const InducedVoltage &iv);

    /// Returns the voltage graph
    TGraph* GetVoltageGraph();

    /// Deletes the TGraph pointer to free up memory
    void ResetVoltage();

    // Returns the input file name
    TString GetFilename() { return theFile; }

    // Returns the final time in the file
    double GetFinalTime();

    // Returns the upper antenna bandwidth
    double GetUpperAntennaBandwidth();

    // Returns the lower antenna bandwidth limit
    double GetLowerAntennaBandwidth();

    // Applies the antenna bandwidth to the signal
    void ApplyAntennaBandwidth();

    /// Returns periodogram
    /// \param loadResistance The resistance of the load circuit
    TGraph* GetPowerPeriodogram(const double loadResistance);

    /// Returns the chunk size used in generating the voltage
    double GetChunkSize() { return chunkSize; }
  };
}


#endif

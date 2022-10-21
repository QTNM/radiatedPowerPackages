// FieldClasses.h
#ifndef FIELD_CLASSES_H
#define FIELD_CLASSES_H

#include <vector>

#include "TFile.h"
#include "TGraph.h"
#include "TVector3.h"
#include "TString.h"
#include "Math/Vector3D.h"
#include "Math/Point3D.h"

#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/EMFunctions.h"
#include "Antennas/IAntenna.h"

namespace rad
{
  class FieldPoint
  {
  private:
    // Input file name
    TString inputFile;
    
    // Time series for field components
    TGraph* EField[3];
    TGraph* BField[3];

    // Time series of electron position
    TGraph* pos[3];

    TGraph* tPrime; // Relationship between time and retarded time

    void ResetFields();

    /// Produces a TGraph of the same variable on the y axis of the input graph but using the retarded time
    /// \param grOriginal The input TGraph which should be a time series
    /// \return A time series graph plotted using the retarded time
    TGraph* MakeRetardedTimeGraph(const TGraph* grOriginal);

    IAntenna* myAntenna; // The chosen antenna, contains position and direction info

    double fileStartTime; // First time in the file (in seconds)

    /// Variables to record the times to which graphs should be trimmed down to
    /// This is necessary when making multiple retarded time graphs
    double minCutTime;
    double maxCutTime;

    /// Function to cut time series graphs down to a specified time
    /// \param grIn The input time series graph
    /// \param minCutTime The minimum time at which the outputted graph should start
    /// \param maxCutTime The maximum time at which the outputted graph should end
    /// \returns The trimmed graph
    TGraph* TrimGraphToTime(const TGraph* grIn);
    
  public:
    enum Coord_t{
      kX, kY, kZ
    };

    /// Parametrised constructor
    /// \param trajectoryFilePath TString to the input ROOT file containing electron trajectory
    /// \param myAntenna Pointer to the chosen antenna
    FieldPoint(TString trajectoryFilePath, IAntenna* myAntenna);

    /// Destructor
    ~FieldPoint();

    /// Copy constructor
    FieldPoint(const FieldPoint &fp);

    /// Returns a pointer to the antenna used for the calculations
    IAntenna* GetAntenna() { return myAntenna; }
    
    /// Fills the class members between two given times
    /// \param minTime The initial time in the input file
    /// \param maxTime The final time to use in the input file
    void GenerateFields(const double minTime, const double maxTime);

    /// Electron position as a function of time
    /// \param coord The chosen cartesian coordinate
    /// \param kUseRetardedTime Boolean to use retarded time
    /// \param firstPoint First point of the class members to return
    /// \param firstPoint Last point of the class members to return
    /// \return The electron coordinate as a function of time
    TGraph* GetPositionTimeDomain(Coord_t coord, const bool kUseRetardedTime=false,
				  int firstPoint=-1, int lastPoint=-1);
    
    /// Electric field component at the antenna as a function of time
    /// \param coord The chosen cartesian coordinate
    /// \param kUseRetardedTime Boolean to use retarded time
    /// \param firstPoint First point of the class members to return
    /// \param firstPoint Last point of the class members to return
    /// \return The electric field as a function of time    
    TGraph* GetEFieldTimeDomain(Coord_t coord, const bool kUseRetardedTime=false,
				int firstPoint=-1, int lastPoint=-1);

    /// Magnetic field component at the antenna as a function of time
    /// \param coord The chosen cartesian coordinate
    /// \param kUseRetardedTime Boolean to use retarded time
    /// \param firstPoint First point of the class members to return
    /// \param firstPoint Last point of the class members to return
    /// \return The magnetic field as a function of time    
    TGraph* GetBFieldTimeDomain(Coord_t coord, const bool kUseRetardedTime=false);

    /// Poynting vector component at the antenna as a function of time
    /// \param coord The chosen cartesian coordinate
    /// \param kUseRetardedTime Boolean to use retarded time
    /// \param firstPoint First point of the class members to return
    /// \param firstPoint Last point of the class members to return
    /// \return The poynting vector as a function of time
    TGraph* GetPoyntingVecTimeDomain(Coord_t coord, const bool kUseRetardedTime=false);
    
    /// Electric field magnitude at the antenna as a function of time
    /// \param kUseRetardedTime Boolean to use retarded time
    /// \return The electric field magnitude as a function of time
    TGraph* GetEFieldMagTimeDomain(const bool kUseRetardedTime=false);

    /// Magnetic field magnitude at the antenna as a function of time
    /// \param kUseRetardedTime Boolean to use retarded time
    /// \return The magnetic field magnitude as a function of time
    TGraph* GetBFieldMagTimeDomain(const bool kUseRetardedTime=false);

    /// Poynting vector magnitude at the antenna as a function of time
    /// \param kUseRetardedTime Boolean to use retarded time
    /// \return The poynting vector magnitude as a function of time
    TGraph* GetPoyntingMagTimeDomain(const bool kUseRetardedTime=false);

    /// Load voltage from the antenna as a function of time (assumes impedance matching)
    /// \param kUseRetardedTime Boolean to use retarded time
    /// \param firstPoint First point of the class members to return
    /// \param lastPoint Last point of the class members to return
    /// \return The load voltage as a function of time 
    TGraph* GetAntennaLoadVoltageTimeDomain(const bool kUseRetardedTime=false,
					    int firstPoint=-1, int lastPoint=-1);

    /// Load power from the antenna as a function of time
    /// \param loadResistance The load resistance
    /// \param kUseRetardedTime Boolean to use retarded time
    /// \param firstPoint First point of the class members to return
    /// \param lastPoint Last point of the class members to return
    /// \return The load voltage as a function of time 
    TGraph* GetAntennaLoadPowerTimeDomain(const double loadResistance,
					  const bool kUseRetardedTime=false,
					  int firstPoint=-1, int lastPoint=-1);

    /// Frequency domain functions

    /// Get the electric field periodogram for a given component
    /// \param coord The chosen component
    /// \param kUseRetardedTime Boolean to use retarded time
    TGraph* GetEFieldPeriodogram(Coord_t coord, const bool kUseRetardedTime=false);

    /// Get the total electric field periodogram for a given component
    /// \param kUseRetardedTime Boolean to use retarded time
    TGraph* GetTotalEFieldPeriodogram(const bool kUseRetardedTime=false);

    /// These have a normalisation such that they shouldn't be altered by varying
    /// things such as zero padding and the sample rate
    /// \param coord The chosen component
    /// \param kUseRetardedTime Boolean to use retarded time
    TGraph* GetEFieldPowerSpectrumNorm(Coord_t coord, const bool kUseRetardedTime=false);

    /// These have a normalisation such that they shouldn't be altered by varying
    /// things such as zero padding and the sample rate
    /// \param coord The chosen component
    /// \param kUseRetardedTime Boolean to use retarded time
    TGraph* GetTotalEFieldPowerSpectrumNorm(const bool kUseRetardedTime=false);

    /// These have a normalisation such that they shouldn't be altered by varying
    /// things such as zero padding and the sample rate
    /// \param resistance The load resistance (assumes matched impedances)
    /// \param kUseRetardedTime Boolean to use retarded time
    TGraph* GetAntennaLoadPowerSpectrumNorm(const double resistance,
					    const bool kUseRetardedTime=false,
					    int firstPoint=-1, int lastPoint=-1);
    
    // Functions for various useful things such as the final time in a file
    double GetFinalTime();
    double GetSampleRate();
  };

}

#endif

/*
  IAntenna.h
  Contains the abstract base class for a generalised antenna
*/

#ifndef IANTENNA_H
#define IANTENNA_H

#include "TVector3.h"

namespace rad {

  class IAntenna {
    
  public:
    ~IAntenna(){}
    
    // Get the radiation patterns
    // Overidden in each antenna derived class
    virtual TVector3 GetETheta(const TVector3 electronPosition) = 0;
    virtual TVector3 GetEPhi(const TVector3 electronPosition) = 0;

    // Return the effective antenna height/length  
    virtual double GetHEff() = 0;
    
    double GetCentralFrequency(){ return centralFreq; }
    
    double GetCentralWavelength();

    TVector3 GetAntennaPosition(){ return antennaPosition; }

    double GetBandwidthUpperLimit() { return upperBandwidth; }

    double GetBandwidthLowerLimit() { return lowerBandwidth; }

    void SetBandwidth(const double lowerLimit=-DBL_MAX, const double upperLimit=DBL_MAX);

    double GetTimeDelay() { return timeDelay; }
    
  protected:
    TVector3 antennaPosition;
    TVector3 antennaXAxis;
    TVector3 antennaYAxis;
    TVector3 antennaZAxis;
    
    /// The upper and lower frequency bandwidths 
    double lowerBandwidth;
    double upperBandwidth;
    
    /// Central frequency of the antenna bandwidth
    double centralFreq;

    /// Time delay for this antenna
    double timeDelay; 
    
    /// \param electronPosition The position of the electron in global coordinates
    /// \returns The unit vector in the theta direction (relative to antenna axis) in global coords
    TVector3 GetThetaHat(const TVector3 electronPosition);

    /// \param electronPosition The position of the electron in global coordinates
    /// \returns The unit vector in the phi direction (relative to antenna axis) in global coords
    TVector3 GetPhiHat(const TVector3 electronPosition);

    /// \param electronPosition The position of the electron in global coordinates
    /// \returns The polar angle theta w.r.t. the antenna axis
    double GetTheta(const TVector3 electronPosition);

    /// \param electronPosition The position of the electron in global coordinates
    /// \returns The azimuthal angle phi w.r.t. the antenna axis
    double GetPhi(const TVector3 electronPosition);
  };
}

#endif

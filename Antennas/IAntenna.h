/*
  IAntenna.h
  Contains the abstract base class for a generalised antenna
*/

#ifndef IANTENNA_H
#define IANTENNA_H

#include "TVector3.h"

namespace rad {

class IAntenna {
 private:
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

 public:
  virtual ~IAntenna() {}

  /////// Getters ////////
  // Get the radiation patterns
  // Overidden in each antenna derived class
  virtual TVector3 GetETheta(const TVector3 electronPosition) = 0;
  virtual TVector3 GetEPhi(const TVector3 electronPosition) = 0;

  /// \param electronPosition The position of the electron in global coordinates
  /// \returns The polar angle theta w.r.t. the antenna axis
  double GetTheta(const TVector3 electronPosition);

  /// \param electronPosition The position of the electron in global coordinates
  /// \returns The azimuthal angle phi w.r.t. the antenna axis
  double GetPhi(const TVector3 electronPosition);

  // Radiation patterns taking angles as arguments
  virtual double GetETheta(double theta, double phi) = 0;
  virtual double GetEPhi(double theta, double phi) = 0;

  // Returns the effective antenna height/length
  virtual double GetHEff() = 0;

  virtual double GetHEff(TVector3 ePos) = 0;

  /// Virtual function for returning antenna effective area
  /// \param electronPosition The position vector of the electron
  /// \return The effective area of the antenna (in metres)
  virtual double GetAEff(TVector3 electronPosition) = 0;

  double GetCentralFrequency() { return centralFreq; }

  double GetCentralWavelength();

  TVector3 GetAntennaPosition() { return antennaPosition; }

  double GetBandwidthUpperLimit() { return upperBandwidth; }

  double GetBandwidthLowerLimit() { return lowerBandwidth; }

  void SetBandwidth(const double lowerLimit = -DBL_MAX,
                    const double upperLimit = DBL_MAX);

  double GetTimeDelay() { return timeDelay; }

  virtual double GetAEffTheta(TVector3 ePos) = 0;

  virtual double GetAEffPhi(TVector3 ePos) = 0;

 protected:
  /// \param electronPosition The position of the electron in global coordinates
  /// \returns The unit vector in the theta direction (relative to antenna axis)
  /// in global coords
  TVector3 GetThetaHat(const TVector3 electronPosition);

  /// \param electronPosition The position of the electron in global coordinates
  /// \returns The unit vector in the phi direction (relative to antenna axis)
  /// in global coords
  TVector3 GetPhiHat(const TVector3 electronPosition);

  /// Gets the integral of the radiation pattern. Used for normalisation
  double GetPatternIntegral();

  ////////// Setters /////////
  /// @brief Sets the antenna position vector
  /// @param pos Position vector to set to
  void SetAntennaPosition(TVector3 pos) { antennaPosition = pos; }

  /// @brief X axis setter
  /// @param ax
  void SetAntennaXAx(TVector3 ax) { antennaXAxis = ax; }

  /// @brief Y axis setter
  /// @param ax
  void SetAntennaYAx(TVector3 ax) { antennaYAxis = ax; }

  /// @brief Z axis setter
  /// @param ax
  void SetAntennaZAx(TVector3 ax) { antennaZAxis = ax; }

  /// @brief Central frequency setter
  /// @param f Frequency in Hz
  void SetCentralFreq(double f) { centralFreq = f; }

  /// @brief Time delay setter
  /// @param t Time delay in seconds
  void SetTimeDelay(double t) { timeDelay = t; }
};

}  // namespace rad

#endif

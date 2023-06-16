/*
  PatchAntenna.h
  Implementation of a microstrip patch antenna of variable width and length
*/

#ifndef PATCH_ANTENNA_H
#define PATCH_ANTENNA_H

#include "Antennas/IAntenna.h"
#include "TVector3.h"

namespace rad
{
  class PatchAntenna : public IAntenna
  {

  private:
    double L;            // Length of the patch (in metres)
    double W;            // Width of the patch antenna (in metres)
    double H;            // Height of the substrate (in metres)
    double relativePerm; // Relative permettivity
    double LEff;         // Effective length (in metres)

    double PRad;         // Integral of radiation pattern over surface
    double directivity;  // Directivity of antenna

  public:
    PatchAntenna(TVector3 antPos, TVector3 antXAx, TVector3 antYAx,
                 double width, double height, double f0,
                 double permittivity = 1.0, double delay = 0.0);

    /// Calculates the electric field of the antenna in the theta direction
    /// \param electronPosition Electron position (in metres)
    /// \return The electric field (in V/m)
    TVector3 GetETheta(const TVector3 electronPosition);

    /// Calculates the electric field of the antenna in the phi direction
    /// \param electronPosition Electron position (in metres)
    /// \return The electric field (in V/m)
    TVector3 GetEPhi(const TVector3 electronPosition);

    /// Get the theta component of the electric field
    double GetETheta(double theta, double phi);

    /// Get the phi component of the electric field
    double GetEPhi(double theta, double phi);

    /// Get the effective height of the antenna
    /// \return The effective height of the patch antenna (in metres) 
    double GetHEff() { return LEff; }

    double GetHEff(TVector3 ePos);

    /// Calculates the effective area for a given electron position
    /// \param ePos Electron position vector
    /// \return The effective area (in metres squared)
    double GetAEff(TVector3 ePos);

    /// Get the length of the patch
    /// \return The length of the patch in metres 
    double GetL() { return L; }

    /// Get the height of the substrate
    /// \return The height of the substrate in metres
    double GetH() { return H; }

    /// Get the width of the patch
    /// \return The width of the patch in metres
    double GetW() { return W; }

    /// Calculate the input impedance of the patch antenna
    /// \return The impedance (in Ohms)
    double GetImpedance();

    /// Calculates patch antenna bandwidth
    /// \return The antenna bandwidth in Hertz 
    double GetBandwidth();

    /// Get the partial polarisation effective area in the theta direction
    /// \param ePos The electron position vector
    /// \return The partial effective area in metres squared
    double GetAEffTheta(TVector3 ePos);

    /// Get the partial polarisation effective area in the phi direction
    /// \param ePos The electron position vector
    /// \return The partial effective area in metres squared
    double GetAEffPhi(TVector3 ePos);
  };
}

#endif

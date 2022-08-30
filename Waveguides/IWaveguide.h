/*
  IWaveguide.h

  Abstract base class for waveguides
*/

#ifndef IWAVEGUIDE_H
#define IWAVEGUIDE_H

#include "TVector3.h"

namespace rad {
  
  class IWaveguide {

  public:
    // Different kinds of modes we may have
    enum Mode_t {
      kTE, kTM, kTEM
    };

    /// Gets the electric field vector for a given mode at a point                               
    /// \param modeType The mode type to get (either TE or TM)                                     
    /// \param n The angular number of the mode                                               
    /// \param m The radial number of the mode                                                
    /// \param pos The position vector (in metres)                                              
    /// \param omega Angular frequency of the chosen wave                                         
    /// \param A Arbitrary amplitude for part of solution (default = 1)                          
    /// \param B Arbitrary amplitude for part of solution (default = 0)                    
    /// \Returns The mode electric field vector at the supplied point                                
    virtual TVector3 GetModeEField(Mode_t modeType, int n, int m, TVector3 pos, double omega, double A=1, double B=0) = 0;

    /// Gets the H field vector for a given mode at a point                                
    /// \param modeType The mode type to get (either TE or TM)                                   
    /// \param n The angular number of the mode                                                   
    /// \param m The radial number of the mode                                                 
    /// \param pos The position vector (in metres)                                     
    /// \param omega Angular frequency of the chosen wave                              
    /// \param A Arbitrary amplitude for part of solution (default = 1)                                 
    /// \param B Arbitrary amplitude for part of solution (default = 0)                        
    /// \Returns The mode H field vector at the supplied point                                        
    virtual TVector3 GetModeHField(Mode_t modeType, int n, int m, TVector3 pos, double omega, double A=1, double B=0) = 0;

    /// Gets the characteristic mode impedance for a given mode                             
    /// \param modeType The mode type to get (either TE or TM)                                   
    /// \param n The angular number of the mode                                                   
    /// \param m The radial number of the mode                                                   
    /// \param omega Angular frequency of the chosen wave                                      
    /// \Returns The impedance of the mode (in Ohms)                                                  
    virtual double GetModeImpedance(Mode_t modeType, int n, int m, double omega) = 0;
  };
}

#endif
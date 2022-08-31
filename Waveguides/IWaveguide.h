/*
  IWaveguide.h

  Abstract base class for waveguides
*/

#ifndef IWAVEGUIDE_H
#define IWAVEGUIDE_H

#include "BasicFunctions/ComplexVector3.h"

#include "TVector3.h"

namespace rad {
  
  class IWaveguide {

  public:
    // Different kinds of modes we may have
    enum Mode_t {
      kTE, kTM, kTEM
    };

    /// Gets the complex electric field vector for a given mode at a point                
    /// \param modeType The mode type to get (TE, TM, TEM)                                     
    /// \param pos The position vector (in metres)                                              
    /// \param omega Angular frequency of the chosen wave                                         
    /// \param A Arbitrary amplitude for part of solution (default = 1)                          
    /// \param B Arbitrary amplitude for part of solution (default = 0)                    
    /// \Returns The mode electric field vector at the supplied point    
    virtual ComplexVector3 GetModeEFieldComplex(Mode_t modeType, int n, int m, TVector3 pos, double omega, double A=1, double B=0) = 0;
    
    /// Gets the electric field vector for a given mode at a point                               
    /// \param modeType The mode type to get (TE, TM, TEM)                                     
    /// \param n The angular number of the mode                                               
    /// \param m The radial number of the mode                                                
    /// \param pos The position vector (in metres)                                              
    /// \param omega Angular frequency of the chosen wave                                         
    /// \param A Arbitrary amplitude for part of solution (default = 1)                          
    /// \param B Arbitrary amplitude for part of solution (default = 0)                    
    /// \Returns The mode electric field vector at the supplied point                                
    virtual TVector3 GetModeEField(Mode_t modeType, int n, int m, TVector3 pos, double omega, double A=1, double B=0) = 0;

    /// Gets the complex magnetic field strength vector for a given mode at a point       
    /// \param modeType The mode type to get (TE, TM, TEM)                                   
    /// \param n The angular number of the mode                                                   
    /// \param m The radial number of the mode                                                 
    /// \param pos The position vector (in metres)                                     
    /// \param omega Angular frequency of the chosen wave                              
    /// \param A Arbitrary amplitude for part of solution (default = 1)                                 
    /// \param B Arbitrary amplitude for part of solution (default = 0)                        
    /// \Returns The mode H field vector at the supplied point
    virtual ComplexVector3 GetModeHFieldComplex(Mode_t modeType, int n, int m, TVector3 pos, double omega, double A=1, double B=0) = 0;
    
    /// Gets the H field vector for a given mode at a point                                
    /// \param modeType The mode type to get (TE, TM, TEM)                                   
    /// \param n The angular number of the mode                                                   
    /// \param m The radial number of the mode                                                 
    /// \param pos The position vector (in metres)                                     
    /// \param omega Angular frequency of the chosen wave                              
    /// \param A Arbitrary amplitude for part of solution (default = 1)                                 
    /// \param B Arbitrary amplitude for part of solution (default = 0)                        
    /// \Returns The mode H field vector at the supplied point                                        
    virtual TVector3 GetModeHField(Mode_t modeType, int n, int m, TVector3 pos, double omega, double A=1, double B=0) = 0;

    /// Gets the characteristic mode impedance for a given mode                             
    /// \param modeType The mode type to get (TE, TM, TEM)                                   
    /// \param n The angular number of the mode                                                   
    /// \param m The radial number of the mode                                                   
    /// \param omega Angular frequency of the chosen wave                                      
    /// \Returns The impedance of the mode (in Ohms)                                                  
    virtual double GetModeImpedance(Mode_t modeType, int n, int m, double omega) = 0;

    /// Gets the cutoff frequency for a particular waveguide mode
    /// \param modeType The mode type to get (TE, TM, TEM)
    /// \param n The first mode index to get
    /// \param m The second mode index to get
    /// \Returns The cutoff frequency in Hertz
    virtual double GetCutoffFrequency(Mode_t modeType, int n, int m) = 0;
  };
}

#endif

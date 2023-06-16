/*
  FieldPointNR.h

  Non-relativistic version of field point
*/

#include "FieldClasses/FieldClasses.h"

namespace rad
{
  class FieldPointNR : public FieldPoint
  {
  public:
    FieldPointNR(TString trajectoryFilePath, IAntenna *myAntenna) : FieldPoint(trajectoryFilePath, myAntenna) {}

    /// Fills the class members between two given times
    /// \param minTime The initial time in the input file
    /// \param maxTime The final time to use in the input file
    void GenerateFields(const double minTime, const double maxTime);
  };
}
// BasicFunctions.h
#ifndef BASIC_FUNCTIONS_H
#define BASIC_FUNCTIONS_H

#include "TVector3.h"
#include "TGraph.h"

namespace rad
{

  TVector3 CalcEField(const TVector3 fieldPoint, const TVector3 ePosition,
		      const TVector3 eVelocity, const TVector3 eAcceleration);

  TVector3 CalcBField(const TVector3 fieldPoint, const TVector3 ePosition,
		      const TVector3 eVelocity, const TVector3 eAcceleration);

  TVector3 CalcPoyntingVec(const TVector3 fieldPoint, const TVector3 ePosition,
			   const TVector3 eVelocity, const TVector3 eAcceleration);

  TVector3 CalcPoyntingVec(const TVector3 EField, const TVector3 BField);

  void setGraphAttr(TGraph *gr);

  double CalcAeHertzianDipole(const double wavelength, const TVector3 dipoleDir,
			      const TVector3 ePosition, const TVector3 position);

  double CalcRetardedTime(const TVector3 fieldPoint, const TVector3 ePosition, const double labTime);
 
}
 
#endif
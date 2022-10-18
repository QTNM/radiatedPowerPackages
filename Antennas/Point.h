/*
    Point.h

    An class which just represents a point in space
    rather than an actual antenna

    Seb Jones 18-10-2022
*/

#ifndef POINT_H
#define POINT_H

#include "Antennas/IAntenna.h"

#include "TVector3.h"

namespace rad
{
    class Point : public IAntenna
    {
    public:
        /// Parametrised constructor
        /// \param antennaPos The position of the antenna 
        Point(TVector3 antennaPos);

        // Used for radiation pattern usually - returns 0 here 
        TVector3 GetETheta(const TVector3 electronPosition);
        TVector3 GetEPhi(const TVector3 electronPosition);
    };
}

#endif
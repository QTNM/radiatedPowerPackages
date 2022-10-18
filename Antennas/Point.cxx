/// Point.cxx

#include "Antennas/Point.h"

rad::Point::Point(TVector3 antennaPos)
{
    antennaPosition = antennaPos;
    antennaXAxis = TVector3(1, 0, 0);
    antennaYAxis = TVector3(0, 1, 0);
    antennaXAxis = TVector3(0, 0, 1);
    centralFreq = 0.0;
    timeDelay = 0.0;
    SetBandwidth();
}

TVector3 rad::Point::GetETheta(TVector3 electronPosition)
{
    return TVector3(0, 0, 0);
}

TVector3 rad::Point::GetEPhi(TVector3 electronPosition)
{
    return TVector3(0, 0, 0);
}
/*
    FieldStorer.h

    Class for storing the advance time fields.
    Designed to be used by the Event class to keep a moving track of fields

    Seb Jones 18-10-2022
*/

#include "BasicFunctions/Constants.h"
#include "EventGeneration/ParticleState.h"

#include "TVector3.h"

#include <vector>

namespace rad
{
    class FieldStorer
    {
    private:
        std::vector<double> tA;
        std::vector<TVector3> eField;
        std::vector<TVector3> bField;

        double DoCubicInterpolation(std::vector<double> xVals,
                                    std::vector<double> yVals, double interp);

    public:
        /// Parametrised constructor
        /// \param eField0 Initial electric field at field point
        /// \param bField0 Initial magnetic field at field point
        /// \param tA0 Initial advanced time at field point
        FieldStorer(TVector3 eField0, TVector3 bField0, double tA0);

        /// Adds new fields and time to the object and (if necessary) removes old ones
        /// \param newEField Electric field vector to be added
        /// \param newBField Magnetic field vector to be added
        /// \param tANew Advanced time to be added (in seconds)
        void AddNewFields(TVector3 newEField, TVector3 newBField, double tANew);

        /// Gets an interpolated electric field at the desired time 
        /// \param timeInterp The time at which to interpolate the value
        /// \return The interpolated field value (in units of V/m)
        TVector3 GetInterpolatedEField(double timeInterp);

        /// Gets an interpolated magnetic field at the desired time 
        /// \param timeInterp The time at which to interpolate the value
        /// \return The interpolated field value (in units of T
        TVector3 GetInterpolatedBField(double timeInterp);
    };
}
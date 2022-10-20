/// FieldStorer.cxx

#include "EventGeneration/FieldStorer.h"

#include "TSpline.h"

#include <iostream>

rad::FieldStorer::FieldStorer(TVector3 eField0, TVector3 bField0,
                              TVector3 pos0, double tA0)
{
    tA.push_back(tA0);
    eField.push_back(eField0);
    bField.push_back(bField0);
    pos.push_back(pos0);
}

void rad::FieldStorer::AddNewFields(TVector3 newEField, TVector3 newBField,
                                    TVector3 newPos, double tANew)
{
    tA.push_back(tANew);
    eField.push_back(newEField);
    bField.push_back(newBField);
    pos.push_back(newPos);

    const double maxTimeLength{5e-8}; // seconds
    // If the vector is too long, remove the point at the start
    if (tA.at(tA.size() - 1) - tA.at(0) > maxTimeLength)
    {
        tA.erase(tA.begin());
        eField.erase(eField.begin());
        bField.erase(bField.begin());
        pos.erase(pos.begin());
    }
}

double rad::FieldStorer::DoCubicInterpolation(std::vector<double> xVals,
                                              std::vector<double> yVals, double interp)
{
    // Use the 3rd integrating Lagrange polynomial
    // First calculate the Lagrange interpolating basis functions
    double l0{(interp - xVals[1]) * (interp - xVals[2]) * (interp - xVals[3]) /
              ((xVals[0] - xVals[1]) * (xVals[0] - xVals[2]) * (xVals[0] - xVals[3]))};
    double l1{(interp - xVals[0]) * (interp - xVals[2]) * (interp - xVals[3]) /
              ((xVals[1] - xVals[0]) * (xVals[1] - xVals[2]) * (xVals[1] - xVals[3]))};
    double l2{(interp - xVals[0]) * (interp - xVals[1]) * (interp - xVals[3]) /
              ((xVals[2] - xVals[0]) * (xVals[2] - xVals[1]) * (xVals[2] - xVals[3]))};
    double l3{(interp - xVals[0]) * (interp - xVals[1]) * (interp - xVals[2]) /
              ((xVals[3] - xVals[0]) * (xVals[3] - xVals[1]) * (xVals[3] - xVals[2]))};

    return l0 * yVals[0] + l1 * yVals[1] + l2 * yVals[2] + l3 * yVals[3];
}

TVector3 rad::FieldStorer::GetInterpolatedEField(double timeInterp)
{
    std::vector<double> timeVals(4);
    std::vector<double> eXVals(4);
    std::vector<double> eYVals(4);
    std::vector<double> eZVals(4);

    // Signal has not yet reached the antenna
    if (timeInterp < tA.at(0))
    {
        return TVector3(0, 0, 0);
    }

    // Loop through the times and find the points which our
    // data point lies between
    for (int i = 0; i < tA.size(); i++)
    {
        if (timeInterp > tA.at(i) && timeInterp < tA.at(i + 1))
        {
            if (i == 0)
            {
                timeVals.at(0) = 0;
                eXVals.at(0) = 0;
                eYVals.at(0) = 0;
                eZVals.at(0) = 0;
            }
            else
            {
                timeVals.at(0) = tA.at(i - 1);
                eXVals.at(0) = eField.at(i - 1).X();
                eYVals.at(0) = eField.at(i - 1).Y();
                eZVals.at(0) = eField.at(i - 1).Z();
            }
            timeVals.at(1) = tA.at(i);
            timeVals.at(2) = tA.at(i + 1);
            timeVals.at(3) = tA.at(i + 2);
            eXVals.at(1) = eField.at(i).X();
            eYVals.at(1) = eField.at(i).Y();
            eZVals.at(1) = eField.at(i).Z();
            eXVals.at(2) = eField.at(i + 1).X();
            eYVals.at(2) = eField.at(i + 1).Y();
            eZVals.at(2) = eField.at(i + 1).Z();
            eXVals.at(3) = eField.at(i + 2).X();
            eYVals.at(3) = eField.at(i + 2).Y();
            eZVals.at(3) = eField.at(i + 2).Z();
            break;
        }
    }

    // Now interpolate each field component
    // double eFieldXInterp{InterpolateCubicSpline(timeVals, eXVals, timeInterp)};
    double eFieldXInterp{DoCubicInterpolation(timeVals, eYVals, timeInterp)};
    double eFieldYInterp{DoCubicInterpolation(timeVals, eYVals, timeInterp)};
    double eFieldZInterp{DoCubicInterpolation(timeVals, eZVals, timeInterp)};

    return TVector3(eFieldXInterp, eFieldYInterp, eFieldZInterp);
}

TVector3 rad::FieldStorer::GetInterpolatedBField(double timeInterp)
{
    std::vector<double> timeVals(4);
    std::vector<double> bXVals(4);
    std::vector<double> bYVals(4);
    std::vector<double> bZVals(4);

    if (timeInterp < tA.at(0))
    {
        return TVector3(0, 0, 0);
    }

    // Loop through the times and find the points which our
    // data point lies between
    for (int i = 0; i < tA.size(); i++)
    {
        if (timeInterp > tA.at(i))
        {
            timeVals.at(0) = tA.at(i - 1);
            timeVals.at(1) = tA.at(i);
            timeVals.at(2) = tA.at(i + 1);
            timeVals.at(3) = tA.at(i + 2);
            bXVals.at(0) = bField.at(i - 1).X();
            bYVals.at(0) = bField.at(i - 1).Y();
            bZVals.at(0) = bField.at(i - 1).Z();
            bXVals.at(1) = bField.at(i).X();
            bYVals.at(1) = bField.at(i).Y();
            bZVals.at(1) = bField.at(i).Z();
            bXVals.at(2) = bField.at(i + 1).X();
            bYVals.at(2) = bField.at(i + 1).Y();
            bZVals.at(2) = bField.at(i + 1).Z();
            bXVals.at(3) = bField.at(i + 2).X();
            bYVals.at(3) = bField.at(i + 2).Y();
            bZVals.at(3) = bField.at(i + 2).Z();
            break;
        }
    }

    // Now interpolate each field component
    double bFieldXInterp{DoCubicInterpolation(timeVals, bXVals, timeInterp)};
    double bFieldYInterp{DoCubicInterpolation(timeVals, bYVals, timeInterp)};
    double bFieldZInterp{DoCubicInterpolation(timeVals, bZVals, timeInterp)};
    return TVector3(bFieldXInterp, bFieldYInterp, bFieldZInterp);
}
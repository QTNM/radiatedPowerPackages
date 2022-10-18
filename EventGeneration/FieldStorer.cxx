/// FieldStorer.cxx

#include "EventGeneration/FieldStorer.h"

#include <iostream>

rad::FieldStorer::FieldStorer(TVector3 eField0, TVector3 bField0, double tA0)
{
    tA.push_back(tA0);
    eField.push_back(eField0);
    bField.push_back(bField0);
}

void rad::FieldStorer::AddNewFields(TVector3 newEField, TVector3 newBField,
                                    double tANew)
{
    tA.push_back(tANew);
    eField.push_back(newEField);
    bField.push_back(newBField);

    const double maxTimeLength{5e-8}; // seconds
    // If the vector is too long, remove the point at the start
    if (tA.at(tA.size() - 1) - tA.at(0) > maxTimeLength)
    {
        tA.erase(tA.begin());
        eField.erase(eField.begin());
        bField.erase(bField.begin());
    }
}

double rad::FieldStorer::InterpolateCubicSpline(std::vector<double> xVals,
                                                std::vector<double> yVals, double interp)
{
    if (xVals.size() != 4 || yVals.size() != 4)
    {
        std::cout << "Vectors of input values are the wrong size!!!\n";
        return 0;
    }
    else
    {
        double pm1{yVals[0]};
        double p0{yVals[1]};
        double p1{yVals[2]};
        double p2{yVals[3]};

        double xm1{xVals[0]};
        double x0{xVals[1]};
        double x1{xVals[2]};
        double x2{xVals[3]};

        double m0{0.5 * ((p1 - p0) / (x1 - x0) + (p0 - pm1) / (x0 - xm1))};
        double m1{0.5 * ((p2 - p1) / (x2 - x1) + (p1 - p0) / (x1 - x0))};
        double value{(2 * p0 + m1 - 2 * p1 + m1) * pow(interp, 2) +
                     (-3 * p0 + 3 * p1 - 2 * m0 - m1) * pow(interp, 2) +
                     m0 * interp + p0};
        return value;
    }
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
        if (timeInterp > tA.at(i))
        {
            timeVals.at(0) = tA.at(i - 1);
            timeVals.at(1) = tA.at(i);
            timeVals.at(2) = tA.at(i + 1);
            timeVals.at(3) = tA.at(i + 2);
            eXVals.at(0) = eField.at(i - 1).X();
            eYVals.at(0) = eField.at(i - 1).Y();
            eZVals.at(0) = eField.at(i - 1).Z();
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
    double eFieldXInterp{InterpolateCubicSpline(timeVals, eXVals, timeInterp)};
    double eFieldYInterp{InterpolateCubicSpline(timeVals, eYVals, timeInterp)};
    double eFieldZInterp{InterpolateCubicSpline(timeVals, eZVals, timeInterp)};
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
    double bFieldXInterp{InterpolateCubicSpline(timeVals, bXVals, timeInterp)};
    double bFieldYInterp{InterpolateCubicSpline(timeVals, bYVals, timeInterp)};
    double bFieldZInterp{InterpolateCubicSpline(timeVals, bZVals, timeInterp)};
    return TVector3(bFieldXInterp, bFieldYInterp, bFieldZInterp);
}
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
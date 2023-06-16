/// ComsolFields.cxx

#include "ElectronDynamics/ComsolFields.h"

#include "TGraph2D.h"
#include "TVector3.h"

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

rad::ComsolField::ComsolField(std::string fieldFile, double centralField)
{
    fieldValues = new TGraph2D();

    // Need to read in field data first of all
    std::ifstream file110mm(fieldFile);
    if (file110mm.is_open())
    {
        std::string line;
        int nLines{0};
        int ix{1};
        int iy{1};

        // Read in line by line
        while (std::getline(file110mm, line))
        {
            nLines++;
            // Exclude the lines containing the headers and stuff
            if (nLines >= 10)
            {
                std::stringstream ss{line};
                std::string r;
                std::string z;
                std::string magB;
                std::getline(ss, r, ',');
                std::getline(ss, z, ',');
                std::getline(ss, magB, ',');
                double magBDouble{std::stod(magB)};
                double zDouble{std::stod(z)};
                double rDouble{std::stod(r)};
                // Convert from mm to m
                zDouble *= 1e-3;
                rDouble *= 1e-3;
                fieldValues->SetPoint(fieldValues->GetN(),
                                      zDouble, rDouble, magBDouble);
            }
        }

        if (centralField == 0.0)
        {
            scaleFactor = 1.0;
        }
        else
        {
            std::cout<<centralField<<std::endl;
            //std::cout<<evaluate_field_at_point(TVector3(0.01, 0.01, 0.01)).Mag()<<std::endl;
            //scaleFactor = centralField / evaluate_field_at_point(TVector3(0, 0.0001, 0.0001)).Mag();
            scaleFactor = centralField;
        }
        std::cout<<"Scale factor is "<<scaleFactor<<std::endl;
    }
    else
    {
        std::cout << "Unable to open field map file. Exiting." << std::endl;
        exit(1);
    }
}

TVector3 rad::ComsolField::evaluate_field_at_point(const TVector3 vec)
{
    double r{vec.Perp()};
    double z{vec.Z()};
    // Consider B field to only be in z direction

    return TVector3(0, 0, scaleFactor * fieldValues->Interpolate(z, r));
}

rad::ComsolField::~ComsolField()
{
    delete fieldValues;
}

rad::ComsolHarmonicField::ComsolHarmonicField(double radius, double current,
                                              std::string fieldFile, double centralField)
{
    fieldValues = new TGraph2D();

    // Need to read in field data first of all
    std::ifstream file110mm(fieldFile);
    if (file110mm.is_open())
    {
        std::string line;
        int nLines{0};
        int ix{1};
        int iy{1};

        // Read in line by line
        while (std::getline(file110mm, line))
        {
            nLines++;
            // Exclude the lines containing the headers and stuff
            if (nLines >= 10)
            {
                std::stringstream ss{line};
                std::string r;
                std::string z;
                std::string magB;
                std::getline(ss, r, ',');
                std::getline(ss, z, ',');
                std::getline(ss, magB, ',');
                double magBDouble{std::stod(magB)};
                double zDouble{std::stod(z)};
                double rDouble{std::stod(r)};
                // Convert from mm to m
                zDouble *= 1e-3;
                rDouble *= 1e-3;
                fieldValues->SetPoint(fieldValues->GetN(),
                                      zDouble, rDouble, magBDouble);
            }
        }

        if (centralField == 0.0)
        {
            scaleFactor = 1.0;
        }
        else
        {
            // scaleFactor = centralField / evaluate_field_at_point(TVector3(0, 0, 0)).Mag();
            scaleFactor = centralField;
        }
    }
    else
    {
        std::cout << "Unable to open field map file. Exiting." << std::endl;
        exit(1);
    }

    coil = CoilField(radius, current, 0.0, MU0);
}

rad::ComsolHarmonicField::~ComsolHarmonicField()
{
    delete fieldValues;
}

TVector3 rad::ComsolHarmonicField::evaluate_field_at_point(const TVector3 vec)
{
    double r{vec.Perp()};
    double z{vec.Z()};
    // Consider B field to only be in z direction
    TVector3 coilBField = coil.evaluate_field_at_point(vec);
    return TVector3(0, 0, scaleFactor * fieldValues->Interpolate(z, r)) - coilBField;
}
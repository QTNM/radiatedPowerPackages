#include "Antennas/AntennaArray.h"

void rad::AntennaArray::AddElement(IAntenna *ant)
{
    elements.push_back(ant);
}

rad::AntennaArray::AntennaArray(std::vector<IAntenna*> antennas)
{
    elements = antennas;
}

rad::IAntenna *rad::AntennaArray::GetAntenna(unsigned int nAntenna)
{
    return elements.at(nAntenna);
}
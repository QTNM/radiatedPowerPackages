/*
    AntennaArray.h

    Class designed to represent an array (phased or otherwise) of
    multiple individual antenna elements

    S. Jones 24-10-2022
*/

#ifndef ANTENNA_ARRAY_H
#define ANTENNA_ARRAY_H

#include "Antennas/IAntenna.h"

#include <vector>

namespace rad
{
    class AntennaArray
    {
    private:
        std::vector<IAntenna *> elements; // Vector containing elements

    public:
        /// Default constructor with no elements
        AntennaArray(){};

        /// Parametrised constructor
        /// \param antennas Vector of antenna elements
        AntennaArray(std::vector<IAntenna *> antennas);

        /// Adds an antenna element to the array
        /// \param ant The element to be added to the array
        void AddElement(IAntenna *ant);

        /// Gets the number of elements in the array
        /// \return Integer number of elements in the array
        unsigned int GetNElements() { return elements.size(); }
    };
}

#endif
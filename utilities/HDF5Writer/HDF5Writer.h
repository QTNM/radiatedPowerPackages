#ifndef HDF5_WRITER_H
#define HDF5_WRITER_H

#include <string>
#include <vector>

#include "H5Cpp.h"
#include "TGraph.h"

namespace rad {
namespace hdf5 {

/// @brief Write a 1D double dataset into an HDF5 group
/// @param group The HDF5 group to write into
/// @param name Dataset name
/// @param data Pointer to the data buffer
/// @param nPoints Number of elements to write
void WriteDataset(H5::Group& group, const std::string& name,
                  const double* data, hsize_t nPoints);

/// @brief Write a TGraph's Y-values as an HDF5 dataset
///
/// Extracts the Y array from the TGraph via GetY() and writes it as a 1D
/// dataset. Optionally computes and attaches a "Time step [seconds]"
/// attribute from the X spacing of the first two points.
///
/// @param group The HDF5 group to write into
/// @param name Dataset name
/// @param graph The TGraph whose Y-values will be written
/// @param writeTimeStep If true, attach a time-step attribute (default: true)
void WriteGraphDataset(H5::Group& group, const std::string& name,
                       const TGraph* graph, bool writeTimeStep = true);

/// @brief Write a scalar double attribute on an HDF5 object
/// @param obj The HDF5 object (group, dataset, or file) to attach to
/// @param name Attribute name
/// @param value The scalar value to write
void WriteScalarAttribute(H5::H5Object& obj, const std::string& name,
                          double value);

/// @brief Write a scalar string attribute on an HDF5 object
/// @param obj The HDF5 object (group, dataset, or file) to attach to
/// @param name Attribute name
/// @param value The string value to write
void WriteStringAttribute(H5::H5Object& obj, const std::string& name,
                          const std::string& value);

}  // namespace hdf5
}  // namespace rad

#endif

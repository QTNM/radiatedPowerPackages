#include "utilities/HDF5Writer/HDF5Writer.h"

namespace rad {
namespace hdf5 {

void WriteDataset(H5::Group& group, const std::string& name,
                  const double* data, hsize_t nPoints) {
  hsize_t dim[] = {nPoints};
  H5::DataSpace dspace(1, dim);

  double fillValue{0};
  H5::DSetCreatPropList plist;
  plist.setFillValue(H5::PredType::NATIVE_DOUBLE, &fillValue);

  H5::DataSet dataset = group.createDataSet(
      name, H5::PredType::NATIVE_DOUBLE, dspace, plist);
  dataset.write(data, H5::PredType::NATIVE_DOUBLE, dspace);
}

void WriteGraphDataset(H5::Group& group, const std::string& name,
                       const TGraph* graph, bool writeTimeStep) {
  const hsize_t nPoints = graph->GetN();
  hsize_t dim[] = {nPoints};
  H5::DataSpace dspace(1, dim);

  double fillValue{0};
  H5::DSetCreatPropList plist;
  plist.setFillValue(H5::PredType::NATIVE_DOUBLE, &fillValue);

  H5::DataSet dataset = group.createDataSet(
      name, H5::PredType::NATIVE_DOUBLE, dspace, plist);
  dataset.write(graph->GetY(), H5::PredType::NATIVE_DOUBLE, dspace);

  if (writeTimeStep && nPoints >= 2) {
    double timeStep = graph->GetPointX(1) - graph->GetPointX(0);
    WriteScalarAttribute(dataset, "Time step [seconds]", timeStep);
  }
}

void WriteScalarAttribute(H5::H5Object& obj, const std::string& name,
                          double value) {
  H5::DataSpace scalarSpace(H5S_SCALAR);
  H5::Attribute attr = obj.createAttribute(
      name, H5::PredType::NATIVE_DOUBLE, scalarSpace);
  attr.write(H5::PredType::NATIVE_DOUBLE, &value);
}

void WriteStringAttribute(H5::H5Object& obj, const std::string& name,
                          const std::string& value) {
  H5::StrType strType(H5::PredType::C_S1, H5T_VARIABLE);
  H5::DataSpace scalarSpace(H5S_SCALAR);
  H5::Attribute attr = obj.createAttribute(name, strType, scalarSpace);
  attr.write(strType, value);
}

}  // namespace hdf5
}  // namespace rad

#ifndef MESHGEOMETRY_H
#define MESHGEOMETRY_H

#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/BulkData.hpp"

#include <Realm.h>
#include <ComputeGeometryAlgorithmDriver.h>

namespace tioga_nalu {

typedef stk::mesh::Field<double, stk::mesh::Cartesian> VectorFieldType;
typedef stk::mesh::Field<double> ScalarFieldType;
    
class MeshGeometry
{

public:
    
    MeshGeometry(
        sierra::nalu::Realm *realm);

    void setup();
    
    void initialize();

    void execute();

private:
    MeshGeometry() = delete;
    MeshGeometry(const MeshGeometry&) = delete;
    
    sierra::nalu::Realm *realm_;
    
    stk::mesh::MetaData& meta_;
    
    stk::mesh::BulkData& bulk_;
    
    sierra::nalu::ComputeGeometryAlgorithmDriver *computeGeometryAlgDriver_;

    std::vector<std::string> interiorPartNameVec_;
    std::vector<std::string> boundaryPartNameVec_;
    
};

}

#endif /* MESHGEOMETRY_H */

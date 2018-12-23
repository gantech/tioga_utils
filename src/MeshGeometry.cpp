#include <MeshGeometry.h>
#include <ComputeGeometryBoundaryAlgorithm.h>
#include <ComputeGeometryInteriorAlgorithm.h>
#include <Enums.h>

// basic c++
#include <map>

namespace tioga_nalu {

MeshGeometry::MeshGeometry(
    sierra::nalu::Realm * realm
) : realm_(realm),
    meta_(realm_->meta_data()),
    bulk_(realm_->bulk_data()),
    interiorPartNameVec_(0),
    boundaryPartNameVec_(0)
{
    computeGeometryAlgDriver_ = new sierra::nalu::ComputeGeometryAlgorithmDriver(*realm_);

    interiorPartNameVec_ = {"blade1Rot-HEX","blade2Rot-HEX","blade3Rot-HEX","hubRot-HEX","hubRot-WEDGE","towerNacelle-HEX","towerNacelle-WEDGE"};
    
    boundaryPartNameVec_ = {"blade_1","blade_2","blade_3","hub_rot","nacelle","tower"};
    
}

void MeshGeometry::setup()
{

    ScalarFieldType& dual_nodal_volume = meta_.declare_field<ScalarFieldType>(
        stk::topology::NODE_RANK, "dual_nodal_volume");
    stk::mesh::put_field(dual_nodal_volume, meta_.universal_part());
    
}

void MeshGeometry::initialize()
{


  sierra::nalu::AlgorithmType algType = sierra::nalu::INTERIOR;  
  for (const auto& targetName : interiorPartNameVec_) {
      auto* part = meta_.get_part(targetName);
      std::map<sierra::nalu::AlgorithmType, sierra::nalu::Algorithm *>::iterator it
          = computeGeometryAlgDriver_->algMap_.find(algType);
      if ( it == computeGeometryAlgDriver_->algMap_.end() ) {
          sierra::nalu::ComputeGeometryInteriorAlgorithm *theAlg
              = new sierra::nalu::ComputeGeometryInteriorAlgorithm(*realm_, part);
          computeGeometryAlgDriver_->algMap_[algType] = theAlg;
      }
      else {
          it->second->partVec_.push_back(part);
      }
  }

  algType = sierra::nalu::WALL;  //All walls for now
  for (const auto& targetName : boundaryPartNameVec_) {
      auto* part = meta_.get_part(targetName);
      std::map<sierra::nalu::AlgorithmType, sierra::nalu::Algorithm *>::iterator it
          = computeGeometryAlgDriver_->algMap_.find(algType);
      if ( it == computeGeometryAlgDriver_->algMap_.end() ) {
          sierra::nalu::ComputeGeometryBoundaryAlgorithm *theAlg
              = new sierra::nalu::ComputeGeometryBoundaryAlgorithm(*realm_, part);
          computeGeometryAlgDriver_->algMap_[algType] = theAlg;
      }
      else {
          it->second->partVec_.push_back(part);
      }
  }
  
  
}

void MeshGeometry::execute()
{

    computeGeometryAlgDriver_->execute();
    
}

}

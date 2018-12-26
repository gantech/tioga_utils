#include "stk_mesh/base/Field.hpp"
#include "master_element/MasterElement.h"

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
    stk::mesh::put_field_on_mesh(dual_nodal_volume, meta_.universal_part(),nullptr);

    ScalarFieldType& ref_dual_nodal_volume = meta_.declare_field<ScalarFieldType>(
        stk::topology::NODE_RANK, "ref_dual_nodal_volume");
    stk::mesh::put_field_on_mesh(ref_dual_nodal_volume, meta_.universal_part(),nullptr);
    
    if ( realm_->realmUsesEdges_ ) {
        const int nDim = meta_.spatial_dimension();
        VectorFieldType& edgeAreaVec_ = meta_.declare_field<VectorFieldType>(stk::topology::EDGE_RANK, "edge_area_vector");
        stk::mesh::put_field_on_mesh(edgeAreaVec_, meta_.universal_part(), nDim, nullptr);
    }

    for (const auto& targetName : boundaryPartNameVec_) {
        auto* targetPart = meta_.get_part(targetName);
        const std::vector<stk::mesh::Part*> & sub_parts = targetPart->subsets();
        for( std::vector<stk::mesh::Part*>::const_iterator i = sub_parts.begin();
             i != sub_parts.end(); ++i )
        {
            stk::mesh::Part * const part = *i ;
            const stk::topology theTopo = part->topology();
    
            // push back the part for book keeping and, later, skin mesh
            realm_->bcPartVec_.push_back(part);
            const int nDim = meta_.spatial_dimension();
            
            // register fields
            sierra::nalu::MasterElement *meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(theTopo);
            const int numScsIp = meFC->numIntPoints_;
            
            sierra::nalu::GenericFieldType *exposedAreaVec_
                = &(meta_.declare_field<sierra::nalu::GenericFieldType>(static_cast<stk::topology::rank_t>( \
                                                                                 meta_.side_rank()), "exposed_area_vector"));
            stk::mesh::put_field_on_mesh(*exposedAreaVec_, *part, nDim*numScsIp , nullptr);
        }
    }


  sierra::nalu::AlgorithmType algType = sierra::nalu::INTERIOR;  
  for (const auto& targetName : interiorPartNameVec_) {
      auto* part = meta_.get_part(targetName);
      realm_->interiorPartVec_.push_back(part);
      const stk::topology the_topo = part->topology();
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
      auto* targetPart = meta_.get_part(targetName);
      const std::vector<stk::mesh::Part*> & sub_parts = targetPart->subsets();
      for( std::vector<stk::mesh::Part*>::const_iterator i = sub_parts.begin();
           i != sub_parts.end(); ++i )
      {
          stk::mesh::Part * const part = *i ;

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
    
    

}

void MeshGeometry::initialize()
{

    computeGeometryAlgDriver_->execute();

    const int ndim = meta_.spatial_dimension();
    ScalarFieldType* dual_nodal_volume = meta_.get_field<ScalarFieldType>(
        stk::topology::NODE_RANK, "dual_nodal_volume");
    ScalarFieldType* ref_dual_nodal_volume = meta_.get_field<ScalarFieldType>(
        stk::topology::NODE_RANK, "ref_dual_nodal_volume");

    stk::mesh::Selector sel = stk::mesh::selectField(*dual_nodal_volume);
    const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);

    for (auto b: bkts) {
        for (size_t in=0; in < b->size(); in++) {
            auto node = (*b)[in];
            double* vol = stk::mesh::field_data(*dual_nodal_volume, node);
            double* refVol = stk::mesh::field_data(*ref_dual_nodal_volume, node);
            *refVol = *vol;
        }
    }
  
}

void MeshGeometry::execute()
{

    computeGeometryAlgDriver_->execute();
    
}

}

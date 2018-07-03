
#include "MeshMotion.h"
#include "MeshRotation.h"
#include "OpenfastFSI.h"

#include "stk_mesh/base/Field.hpp"

namespace tioga_nalu {

MeshMotion::MeshMotion(
    stk::mesh::MetaData& meta,
    stk::mesh::BulkData& bulk,
    const YAML::Node& node
) : meta_(meta),
    bulk_(bulk)
{
    load(node);
}

void MeshMotion::load(const YAML::Node& node)
{
    const auto& minfo = node["motion_group"];

    const int num_groups = minfo.size();
    meshMotionVec_.resize(num_groups);

    std::cout << "Number of groups = " << num_groups << std::endl ;
    for (int i=0; i < num_groups; i++) {
        const auto& motion_def = minfo[i];
        std::string type = "rotation";
        if (motion_def["type"]) type = motion_def["type"].as<std::string>();

        if (type == "rotation") {
            std::cout << "Rotation type found " << std::endl ;
            meshMotionVec_[i].reset(new MeshRotation(meta_, bulk_, motion_def));
        } else if (type == "openfastfsi") {
            std::cout << "OpenFAST FSI found " << std::endl ;
            meshMotionVec_[i].reset(new OpenfastFSI(meta_, bulk_, motion_def));
        } else {
            throw std::runtime_error("MeshMotion: Invalid mesh motion type: " + type);
        }
    }

    if (node["start_time"])
        startTime_ = node["start_time"].as<double>();
    if (node["num_time_steps"]) {
        numSteps_ = node["num_time_steps"].as<int>();
        deltaT_ = node["delta_t"].as<double>();
    }
}

void MeshMotion::setup()
{
    VectorFieldType& coordinates = meta_.declare_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    VectorFieldType& current_coordinates = meta_.declare_field<VectorFieldType>(
        stk::topology::NODE_RANK, "current_coordinates");
    VectorFieldType& mesh_displacement = meta_.declare_field<VectorFieldType>(
        stk::topology::NODE_RANK, "mesh_displacement");

    stk::mesh::put_field(coordinates, meta_.universal_part());
    stk::mesh::put_field(current_coordinates, meta_.universal_part());
    stk::mesh::put_field(mesh_displacement, meta_.universal_part());

    const int nDim = meta_.spatial_dimension();
    std::vector<std::string> partNameVec = {"blade1", "blade2", "blade3", "hub", "nacelle", "tower"};
    VectorFieldType *fsiForce = &(meta_.declare_field<VectorFieldType>(stk::topology::NODE_RANK, "fsi_force"));
    for (std::vector<std::string>::iterator it = partNameVec.begin() ; it != partNameVec.end(); ++it) {
        auto * part = meta_.get_part(*it);
        stk::mesh::put_field(*fsiForce, *part, nDim);
    }
    
    for (auto& mm: meshMotionVec_)
        mm->setup();
}

void MeshMotion::initialize()
{
    init_coordinates();
    //TODO: Create some if condition when OpenfastFSI is chosen as a part of MeshMotion
    create_sample_force_field();
    for (auto& mm: meshMotionVec_)
        mm->initialize(startTime_);
}

void MeshMotion::execute(const int istep)
{
    const double curr_time = startTime_ + (istep + 1) * deltaT_;
    currentTime_ = curr_time;

    for (auto& mm: meshMotionVec_)
        mm->execute(curr_time);
}

void MeshMotion::init_coordinates()
{
    const int ndim = meta_.spatial_dimension();
    VectorFieldType* modelCoords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    VectorFieldType* currCoords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "current_coordinates");

    stk::mesh::Selector sel = meta_.universal_part();
    const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);

    for (auto b: bkts) {
        for (size_t in=0; in < b->size(); in++) {
            auto node = (*b)[in];
            double* oldxyz = stk::mesh::field_data(*modelCoords, node);
            double* xyz = stk::mesh::field_data(*currCoords, node);

            for (int d=0; d < ndim; d++)
                xyz[d] = oldxyz[d];
        }
    }
}


void MeshMotion::create_sample_force_field() {
    
    VectorFieldType *fsiForce = meta_.get_field<VectorFieldType>(stk::topology::NODE_RANK, "fsi_force");

    //For a 100m blade with uniform square cross section of 1m x 1m and a resolution of 21 (around cross section) x 101 (span), set a force corresponding to a uniform load of 5e-4 N/m along the blade.
    std::vector<std::string> bladePartNameVec = {"blade1", "blade2", "blade3"};
    for (std::vector<std::string>::iterator it = bladePartNameVec.begin() ; it != bladePartNameVec.end(); ++it) {
        auto * part = meta_.get_part(*it);
        stk::mesh::Selector sel(*part);
        std::cout << "Setting force at nodes " << std::endl ;        
        const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);
        for (auto b: bkts) {
            for (size_t in=0; in < b->size(); in++) {
                auto node = (*b)[in];
                double *fsiForceNode = stk::mesh::field_data(*fsiForce, node);
                std::cout << "Setting force at node " << in << std::endl ;
                fsiForceNode[0] = 2.3573785950023575e-05;
            }
        }
    }
    
}
    
} // tioga_nalu

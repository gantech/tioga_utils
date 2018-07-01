#include "fsiTurbine.h"

#include "stk_mesh/base/Field.hpp"

namespace tioga_nalu {

fsiTurbine::fsiTurbine(
    int iTurb,
    const YAML::Node & node,
    stk::mesh::MetaData & meta,
    stk::mesh::BulkData & bulk
) :
meta_(meta),
bulk_(bulk),
turbineProc_(-1),
turbineInProc_(false),
twrLoadMap_(NULL),
twrDispMap_(NULL),
fsiForce_(NULL)
{

    if(node["tower_part"]) {
        //TODO: Get rid of restriction to allow only one part name for hub
        twrPartName_ = node["tower_part"].as<std::string>();
    } else {
        std::cout << "Tower part name not specified for turbine " << iTurb << std::endl;
    }
    if(node["nacelle_part"]) {
        //TODO: Get rid of restriction to allow only one part name for hub
        nacellePartName_ = node["nacelle_part"].as<std::string>();
    } else {
        std::cout << "Nacelle part name not specified for turbine " << iTurb << std::endl;
    }        
    if(node["hub_part"]) {
        //TODO: Get rid of restriction to allow only one part name for hub
        hubPartName_ = node["hub_part"].as<std::string>();
    } else {
        std::cout << "Hub part names not specified for turbine " << iTurb << std::endl;
    }
    if(node["blade_parts"]) {
        //TODO: Get rid of restriction to allow only one part name per blade
        bladePartNames_ = node["blade_parts"].as<std::vector<std::string>>();
        nBlades_ = bladePartNames_.size();
        bladePartVec_.resize(nBlades_);
        bldLoadMap_.resize(nBlades_);        
        bldDispMap_.resize(nBlades_);
        bldDispMapInterp_.resize(nBlades_);
    } else {
        std::cout << "Blade part names not specified for turbine " << iTurb << std::endl;
    }
    
}

fsiTurbine::~fsiTurbine() {

    //Nothing to do here so far
    
}

void fsiTurbine::setup() {

    //TODO: Check if any part of the turbine surface is on this processor and set turbineInProc_ to True/False

    //TODO:: Figure out a way to check the consistency between the number of blades specified in the Nalu input file and the number of blades in the OpenFAST model.
    
    fsiForce_ = meta_.get_field<VectorFieldType>(stk::topology::NODE_RANK, "fsi_force");
    if (fsiForce_ == NULL)
        throw std::runtime_error("fsiTurbine:: Cannot find NODE_RANK field named 'fsi_force'.");
    

    twrPart_ = meta_.get_part(twrPartName_);
    if (twrPart_ == NULL)
        throw std::runtime_error("fsiTurbine:: No part found for mesh part corresponding to " + twrPartName_);
    
    hubPart_ = meta_.get_part(hubPartName_);
    if (hubPart_ == NULL)
        throw std::runtime_error("fsiTurbine:: No part found for mesh part corresponding to " + hubPartName_);
    
    nacellePart_ = meta_.get_part(nacellePartName_);
    if (nacellePart_ == NULL)
        throw std::runtime_error("fsiTurbine:: No part found for mesh part corresponding to " + nacellePartName_);
    
    for (int i=0; i < nBlades_; i++) {
        bladePartVec_[i] = meta_.get_part(nacellePartName_);
        if (bladePartVec_[i] == NULL)
            throw std::runtime_error("fsiTurbine:: No part found for mesh part corresponding to " + bladePartNames_[i]);
        
    }
        
}

void fsiTurbine::initialize() {

    //Allocate memory for loads and deflections data
}

//! Convert pressure and viscous/turbulent stress on the turbine surface CFD mesh into a "fsiForce" field on the turbine surface CFD mesh
void fsiTurbine::computeFSIforce() {

    
}

//! Map loads from the "fsiForce" field on the turbine surface CFD mesh into point load array that gets transferred to openfast
void fsiTurbine::mapLoads() {

    //To implement this function - assume that 'bldLoadMap_' field contains the node id along the blade or the tower that will accumulate the load corresponding to the node on the CFD surface mesh

    computeFSIforce();
    
}

//! Map the deflections from the openfast nodes to the turbine surface CFD mesh. Will call 'computeDisplacement' for each node on the turbine surface CFD mesh.
void fsiTurbine::mapDisplacements() {
    
   //To implement this function - assume that 'bldDispMap_' field contains the lower node id of the corresponding element along the blade or the tower along with the 'bldDispMapInterp_' field that contains the non-dimensional location of the CFD surface mesh node on that element.

    // For e.g., for blade 'k' if the lower node id from 'bldDisMap_' is 'j' and the non-dimenional location from 'bldDispMapInterp_' is 'm', then the translational displacement for the CFD surface mesh is
    // (1-m) * bld_def[k][j*6+0] + m * bld_def[k][(j+1)*6+0]
    // (1-m) * bld_def[k][j*6+1] + m * bld_def[k][(j+1)*6+1]
    // (1-m) * bld_def[k][j*6+2] + m * bld_def[k][(j+1)*6+2]

    //TODO: When the turbine is rotating, displacement of the surface of the blades and hub (not the nacelle and tower) should be split into a rigid body motion due to the rotation of the turbine, yawing of the turbine and a deflection of the structure itself. The yaw rate and rotation rate will vary over time. OpenFAST always stores the displacements of all nodes with respect to the reference configuration. When using a sliding mesh interface (yaw = 0), the mesh blocks inside the rotating part of the sliding interface should be moved with rigid body motion corresponding to the rotation rate first and then a second mesh deformation procedure should be performed to apply the remaining structural deflection. Figure out how to do this.
}

//! Convert one array of 6 deflections (transX, transY, transZ, wmX, wmY, wmZ) into one vector of translational displacement at a given node on the turbine surface CFD mesh.
void fsiTurbine::computeDisplacement() {

    
}

//! Map each node on the turbine surface CFD mesh to 
void fsiTurbine::computeMapping() {

    
}

double fsiTurbine::projectPt2Line(std::vector<double> & pt, std::vector<double> & lStart, std::vector<double> & lEnd) {
    
}
    
    
}

    

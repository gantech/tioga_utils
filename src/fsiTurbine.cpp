#include "fsiTurbine.h"

namespace tioga_nalu {

fsiTurbine::fsiTurbine(int iTurb, const YAML::Node & node) :
twrLoadMap_(NULL),
twrDispMap_(NULL),
fsiForce_(NULL)
{

    if(node["tower_part"]) {
        twrPartName_ = node["tower_part"].as<std::string>();
    } else {
        std::cout << "Tower part name not specified for turbine " << iTurb << std::endl;
    }
    if(node["nacelle_part"]) {
        nacellePartName_ = node["nacelle_part"].as<std::string>();
    } else {
        std::cout << "Nacelle part name not specified for turbine " << iTurb << std::endl;
    }        
    if(node["hub_part"]) {
        hubPartName_ = node["hub_part"].as<std::string>();
    } else {
        std::cout << "Hub part names not specified for turbine " << iTurb << std::endl;
    }
    if(node["blade_parts"]) {
        bladePartNames_ = node["blade_parts"].as<std::vector<std::string>>();
    } else {
        std::cout << "Blade part names not specified for turbine " << iTurb << std::endl;
    }
    
}

fsiTurbine::~fsiTurbine() {

    //Nothing to do here so far
    
}

//! Convert pressure and viscous/turbulent stress on the turbine surface CFD mesh into a "fsiForce" field on the turbine surface CFD mesh
void fsiTurbine::computeFSIforce() {

    
}

//! Transfer loads from the "fsiForce" field on the turbine surface CFD mesh into point load array that gets transferred to openfast
void fsiTurbine::transferLoads() {

    //To implement this function - assume that 'bldLoadMap_' field contains the node id along the blade or the tower that will accumulate the load corresponding to the node on the CFD surface mesh
    
}

//! Transfer the deflections from the openfast nodes to the turbine surface CFD mesh. Will call 'computeDisplacement' for each node on the turbine surface CFD mesh.
void fsiTurbine::transferDisplacements() {

   //To implement this function - assume that 'bldDispMap_' field contains the lower node id of the corresponding element along the blade or the tower along with the 'bldDispMapInterp_' field that contains the non-dimensional location of the CFD surface mesh node on that element.

    // For e.g., for blade 'k' if the lower node id from 'bldDisMap_' is 'j' and the non-dimenional location from 'bldDispMapInterp_' is 'm', then the translational displacement for the CFD surface mesh is
    // (1-m) * bld_def[k][j*6+0] + m * bld_def[k][j*6+0]
    // (1-m) * bld_def[k][j*6+1] + m * bld_def[k][j*6+1]
    // (1-m) * bld_def[k][j*6+2] + m * bld_def[k][j*6+2]
}

//! Convert one array of 6 deflections (transX, transY, transZ, wmX, wmY, wmZ) into one vector of translational displacement at a given node on the turbine surface CFD mesh.
void fsiTurbine::computeDisplacement() {

    
}

//! Map each node on the turbine surface CFD mesh to 
void fsiTurbine::computeMapping() {

    
}

}

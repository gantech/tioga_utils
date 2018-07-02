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

    const int ndim = meta_.spatial_dimension();
    VectorFieldType* modelCoords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");

    // Do the tower first
    stk::mesh::Selector sel(*twrPart_);
    const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);

    for (auto b: bkts) {
        for (size_t in=0; in < b->size(); in++) {
            auto node = (*b)[in];
            double* xyz = stk::mesh::field_data(*modelCoords, node);
            int* loadMapNode = stk::mesh::field_data(*twrLoadMap_, node);
            int* dispMapNode = stk::mesh::field_data(*twrDispMap_, node);
            double* dispMapInterpNode = stk::mesh::field_data(*twrDispMapInterp_, node);
            std::vector<double> ptCoords(0.0, ndim);
            for(int i=0; i < ndim; i++)
                ptCoords[i] = xyz[i];
            double nDimCoord = -1.0;
            int nPtsTwr = params_.nBRfsiPtsTwr;
            for (int i=0; i < nPtsTwr-1; i++) {
                std::vector<double> lStart = {brFSIdata_.twr_ref_pos[i*6], brFSIdata_.twr_ref_pos[i*6+1], brFSIdata_.twr_ref_pos[i*6+2]};
                std::vector<double> lEnd = {brFSIdata_.twr_ref_pos[(i+1)*6], brFSIdata_.twr_ref_pos[(i+1)*6+1], brFSIdata_.twr_ref_pos[(i+1)*6+2]};
                nDimCoord = projectPt2Line(ptCoords, lStart, lEnd);

                if ((nDimCoord > 0) && (nDimCoord < 1.0)) {
                    *dispMapInterpNode = nDimCoord;
                    *dispMapNode = i;
                    *loadMapNode = i + std::round(nDimCoord);
                }
            }

            //If no element in the OpenFAST mesh contains this node do some sanity check on the perpendicular distance between the surface mesh node and the line joining the ends of the tower
            if (nDimCoord < 1.0)  {
                std::vector<double> lStart = {brFSIdata_.twr_ref_pos[0], brFSIdata_.twr_ref_pos[1], brFSIdata_.twr_ref_pos[2]};
                std::vector<double> lEnd = {brFSIdata_.twr_ref_pos[(nPtsTwr-1)*6], brFSIdata_.twr_ref_pos[(nPtsTwr-1)*6+1], brFSIdata_.twr_ref_pos[(nPtsTwr-1)*6+2]};
                double perpDist = perpProjectDist_Pt2Line(ptCoords, lStart, lEnd);
                if (perpDist > 0.2) {// Something's wrong if a node on the surface mesh of the tower is more than 20% of the tower length away from the tower axis. 
                    throw std::runtime_error("Can't find a projection for point (" + std::to_string(ptCoords[0]) + "," + std::to_string(ptCoords[1]) + "," + std::to_string(ptCoords[2]) + ") on the tower on turbine " + std::to_string(params_.TurbID) + ". The tower extends from " + std::to_string(lStart[0]) + "," + std::to_string(lStart[1]) + "," + std::to_string(lStart[2]) + ") to " + std::to_string(lEnd[0]) + "," + std::to_string(lEnd[1]) + "," + std::to_string(lEnd[2]) + "). Are you sure the initial position and orientation of the mesh is consistent with the input file parameters and the OpenFAST model.");
                } else { //Assign this node to the last point and element of the OpenFAST mesh
                    *dispMapInterpNode = 1.0;
                    *dispMapNode = nPtsTwr-2;
                    *loadMapNode = nPtsTwr-1;
                }
            }
        }
    }

    // Now the blades
    int nBlades = params_.numBlades;
    for (int iBlade=0; iBlade < nBlades; iBlade++) {
        int iStart = 0;
        int nPtsBlade = params_.nBRfsiPtsBlade[iBlade];        
        stk::mesh::Selector sel(*bladePartVec_[iBlade]);
        const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);
        
        for (auto b: bkts) {
            for (size_t in=0; in < b->size(); in++) {
                auto node = (*b)[in];
                double* xyz = stk::mesh::field_data(*modelCoords, node);
                int* loadMapNode = stk::mesh::field_data(*bldLoadMap_, node);
                int* dispMapNode = stk::mesh::field_data(*bldDispMap_, node);
                double* dispMapInterpNode = stk::mesh::field_data(*bldDispMapInterp_, node);
                std::vector<double> ptCoords(0.0, ndim);
                for(int i=0; i < ndim; i++)
                    ptCoords[i] = xyz[i];
                double nDimCoord = -1.0;
                for (int i=0; i < nPtsBlade-1; i++) {
                    std::vector<double> lStart = {brFSIdata_.bld_ref_pos[(iStart+i)*6], brFSIdata_.bld_ref_pos[(iStart+i)*6+1], brFSIdata_.bld_ref_pos[(iStart+i)*6+2]};
                    std::vector<double> lEnd = {brFSIdata_.bld_ref_pos[(iStart+i+1)*6], brFSIdata_.bld_ref_pos[(iStart+i+1)*6+1], brFSIdata_.bld_ref_pos[(iStart+i+1)*6+2]};
                    nDimCoord = projectPt2Line(ptCoords, lStart, lEnd);
                    
                    if ((nDimCoord > 0) && (nDimCoord < 1.0)) {
                        *dispMapInterpNode = nDimCoord;
                        *dispMapNode = i;
                        *loadMapNode = i + std::round(nDimCoord);
                    }
                }
                
                //If no element in the OpenFAST mesh contains this node do some sanity check on the perpendicular distance between the surface mesh node and the line joining the ends of the blade
                if (nDimCoord < 1.0)  {
                    std::vector<double> lStart = {brFSIdata_.bld_ref_pos[iStart*6], brFSIdata_.bld_ref_pos[iStart*6+1], brFSIdata_.bld_ref_pos[iStart*6+2]};
                    std::vector<double> lEnd = {brFSIdata_.bld_ref_pos[(iStart+nPtsBlade-1)*6], brFSIdata_.bld_ref_pos[(iStart+nPtsBlade-1)*6+1], brFSIdata_.twr_ref_pos[(iStart+nPtsBlade-1)*6+2]};
                    double perpDist = perpProjectDist_Pt2Line(ptCoords, lStart, lEnd);
                    if (perpDist > 0.2) {// Something's wrong if a node on the surface mesh of the blade is more than 20% of the blade length away from the blade axis. 
                        throw std::runtime_error("Can't find a projection for point (" + std::to_string(ptCoords[0]) + "," + std::to_string(ptCoords[1]) + "," + std::to_string(ptCoords[2]) + ") on blade " + std::to_string(iBlade) + " on turbine " + std::to_string(params_.TurbID) + ". The blade extends from " + std::to_string(lStart[0]) + "," + std::to_string(lStart[1]) + "," + std::to_string(lStart[2]) + ") to " + std::to_string(lEnd[0]) + "," + std::to_string(lEnd[1]) + "," + std::to_string(lEnd[2]) + "). Are you sure the initial position and orientation of the mesh is consistent with the input file parameters and the OpenFAST model.");
                    } else { //Assign this node to the last point and element of the OpenFAST mesh                        
                        *dispMapInterpNode = 1.0;
                        *dispMapNode = nPtsBlade-2;
                        *loadMapNode = nPtsBlade-1;
                    }
                        
                }
            }
        }
        iStart += nPtsBlade;
    }
    

    
}

/** Project a point 'pt' onto a line from 'lStart' to 'lEnd' and return the non-dimensional location of the projected point along the line in [0-1] coordinates
    \f[ 
    nonDimCoord = \frac{ (\vec{pt} - \vec{lStart}) \cdot ( \vec{lEnd} - \vec{lStart} ) }{ (\vec{lEnd} - \vec{lStart}) \cdot (\vec{lEnd} - \vec{lStart}) }
    \f]
*/
double fsiTurbine::projectPt2Line(std::vector<double> & pt, std::vector<double> & lStart, std::vector<double> & lEnd) {
    
    double nonDimCoord = 0.0;
    
    double num = 0.0;
    double denom = 0.0;
    
    for (int i=0; i < 3; i++) {
        num += (pt[i] - lStart[i]) * (lEnd[i] - lStart[i]) ;
        denom += (lEnd[i] - lStart[i]) * (lEnd[i] - lStart[i]) ;
    }
    
    nonDimCoord = num/denom;
    return nonDimCoord;
}

/** Project a point 'pt' onto a line from 'lStart' to 'lEnd' and return the non-dimensional distance of 'pt' from the line w.r.t the distance from 'lStart' to 'lEnd'
    \f[ 
    \vec{perp} &= (\vec{pt} - \vec{lStart}) - \frac{ (\vec{pt} - \vec{lStart}) \cdot ( \vec{lEnd} - \vec{lStart} ) }{ (\vec{lEnd} - \vec{lStart}) \cdot (\vec{lEnd} - \vec{lStart}) } ( \vec{lEnd} - \vec{lStart} ) \ \
    nonDimPerpDist = \frac{\lvert \vec{perp} \rvert}{ \lvert  (\vec{lEnd} - \vec{lStart}) \rvert }
    \f]
*/
double fsiTurbine::perpProjectDist_Pt2Line(std::vector<double> & pt, std::vector<double> & lStart, std::vector<double> & lEnd) {
    
    double nonDimCoord = 0.0;
    double num = 0.0;
    double denom = 0.0;
    for (int i=0; i < 3; i++) {
        num += (pt[i] - lStart[i]) * (lEnd[i] - lStart[i]) ;
        denom += (lEnd[i] - lStart[i]) * (lEnd[i] - lStart[i]) ;
    }
    nonDimCoord = num/denom;

    double nonDimPerpDist = 0.0;
    for(int i=0; i < 3; i++) {
        double tmp = (pt[i] - lStart[i]) - nonDimCoord * (lEnd[i] - lStart[i]) ;
        nonDimPerpDist += tmp * tmp ;
    }
    nonDimPerpDist = sqrt(nonDimPerpDist/denom) ;
    
    return nonDimPerpDist;
}
    
    
}

    

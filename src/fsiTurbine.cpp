#include "fsiTurbine.h"

#include "stk_util/parallel/ParallelReduce.hpp"
#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_mesh/base/Field.hpp"
#include <cmath>

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
        bladePartVec_[i] = meta_.get_part(bladePartNames_[i]);
        if (bladePartVec_[i] == NULL)
            throw std::runtime_error("fsiTurbine:: No part found for mesh part corresponding to " + bladePartNames_[i]);
        
    }

    twrLoadMap_ = meta_.get_field<ScalarIntFieldType>(stk::topology::NODE_RANK, "twr_load_map");
    if (twrLoadMap_ == NULL)
        twrLoadMap_ =  &(meta_.declare_field<ScalarIntFieldType>(stk::topology::NODE_RANK, "twr_load_map"));
    stk::mesh::put_field(*twrLoadMap_, *twrPart_);

    twrDispMap_ = meta_.get_field<ScalarIntFieldType>(stk::topology::NODE_RANK, "twr_disp_map");
    if (twrDispMap_ == NULL)
        twrDispMap_ =  &(meta_.declare_field<ScalarIntFieldType>(stk::topology::NODE_RANK, "twr_disp_map"));
    stk::mesh::put_field(*twrDispMap_, *twrPart_);

    twrDispMapInterp_ = meta_.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "twr_disp_map_interp");
    if (twrDispMapInterp_ == NULL)
        twrDispMapInterp_ =  &(meta_.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "twr_disp_map_interp"));
    stk::mesh::put_field(*twrDispMapInterp_, *twrPart_);
    
    bldLoadMap_ = meta_.get_field<ScalarIntFieldType>(stk::topology::NODE_RANK, "bld_load_map");
    if (bldLoadMap_ == NULL)
        bldLoadMap_ =  &(meta_.declare_field<ScalarIntFieldType>(stk::topology::NODE_RANK, "bld_load_map"));
    for (int i=0; i < nBlades_; i++)
        stk::mesh::put_field(*bldLoadMap_, *bladePartVec_[i]);

    bldDispMap_ = meta_.get_field<ScalarIntFieldType>(stk::topology::NODE_RANK, "bld_disp_map");
    if (bldDispMap_ == NULL)
        bldDispMap_ =  &(meta_.declare_field<ScalarIntFieldType>(stk::topology::NODE_RANK, "bld_disp_map"));
    for (int i=0; i < nBlades_; i++)
        stk::mesh::put_field(*bldDispMap_, *bladePartVec_[i]);

    bldDispMapInterp_ = meta_.get_field<ScalarFieldType>(stk::topology::NODE_RANK, "bld_disp_map_interp");
    if (bldDispMapInterp_ == NULL)
        bldDispMapInterp_ =  &(meta_.declare_field<ScalarFieldType>(stk::topology::NODE_RANK, "bld_disp_map_interp"));
    for (int i=0; i < nBlades_; i++)
        stk::mesh::put_field(*bldDispMapInterp_, *bladePartVec_[i]);
    
}

void fsiTurbine::initialize() {

    //Allocate memory for loads and deflections data

    int nTwrPts = params_.nBRfsiPtsTwr;
    int nBlades = params_.numBlades;
    int nTotBldPts = 0;
    for (int i=0; i < nBlades; i++)
        nTotBldPts += params_.nBRfsiPtsBlade[i];
    brFSIdata_.twr_ref_pos.resize(6*nTwrPts);
    brFSIdata_.twr_def.resize(6*nTwrPts);
    brFSIdata_.twr_vel.resize(6*nTwrPts);
    brFSIdata_.twr_ld.resize(6*nTwrPts);
    brFSIdata_.bld_ref_pos.resize(6*nTotBldPts);
    brFSIdata_.bld_def.resize(6*nTotBldPts);
    brFSIdata_.bld_vel.resize(6*nTotBldPts);
    brFSIdata_.bld_ld.resize(6*nTotBldPts);
    brFSIdata_.hub_ref_pos.resize(6);
    brFSIdata_.hub_def.resize(6);
    brFSIdata_.hub_vel.resize(6);
    brFSIdata_.nac_ref_pos.resize(6);
    brFSIdata_.nac_def.resize(6);
    brFSIdata_.nac_vel.resize(6);
    
}

//! Convert pressure and viscous/turbulent stress on the turbine surface CFD mesh into a "fsiForce" field on the turbine surface CFD mesh
void fsiTurbine::computeFSIforce() {

    
}

//! Map loads from the "fsiForce" field on the turbine surface CFD mesh into point load array that gets transferred to openfast
void fsiTurbine::mapLoads() {

    //To implement this function - assume that 'bldLoadMap_' field contains the node id along the blade or the tower that will accumulate the load corresponding to the node on the CFD surface mesh

    computeFSIforce();


    //First zero out forces on the OpenFAST mesh
    for (size_t i=0; i < params_.nBRfsiPtsTwr; i++) {
        for (size_t j=0; j < 6; j++)
            brFSIdata_.twr_ld[i*6+j] = 0.0;
    }
        
    int nBlades = params_.numBlades;
    int iRunTot = 0;
    for (size_t iBlade=0; iBlade < nBlades; iBlade++) {
        int nPtsBlade = params_.nBRfsiPtsBlade[iBlade];
        for (size_t i=0; i < nPtsBlade; i++) {
            for (size_t j=0; j < 6; j++)
                brFSIdata_.bld_ld[iRunTot*6+j] = 0.0;
            iRunTot++;
        }
            
    }

    // Now map loads
    const int ndim = meta_.spatial_dimension();
    VectorFieldType* modelCoords = meta_.get_field<VectorFieldType>(stk::topology::NODE_RANK, "coordinates");

    std::vector<double> tmpNodePos(3,0.0); // Vector to temporarily store a position vector
    
    // Do the tower first
    stk::mesh::Selector sel(meta_.locally_owned_part() & (*twrPart_) );
    const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);
    for (auto b: bkts) {
        for (size_t in=0; in < b->size(); in++) {            
            auto node = (*b)[in];
            double* fsiForceNode = stk::mesh::field_data(*fsiForce_, node);
            double* xyz = stk::mesh::field_data(*modelCoords, node);
            int* mapNode = stk::mesh::field_data(*twrDispMap_, node);
            double *mapInterpNode = stk::mesh::field_data(*twrDispMapInterp_, node);

            //Find the interpolated reference position first
            linInterpVec(&brFSIdata_.twr_ref_pos[(*mapNode)*6], &brFSIdata_.twr_ref_pos[(*mapNode + 1)*6], *mapInterpNode, tmpNodePos.data());

            std::vector<double> tmpForceMoment(6,0.0); // Temporarily store total force and moment as (fX, fY, fZ, mX, mY, mZ)
            //Now compute the force and moment on the interpolated reference position
            computeEffForceMoment(fsiForceNode, xyz, tmpForceMoment.data(), tmpNodePos.data());

            //Split the force and moment into the two surrounding nodes in a variationally consistent manner using the interpolation factor
            splitForceMoment(tmpForceMoment.data(), *mapInterpNode, &(brFSIdata_.twr_ld[(*mapNode)*6]), &(brFSIdata_.twr_ld[(*mapNode+1)*6]));
            
        }
    }

    // Now the blades
    int iStart = 0;
    for (int iBlade=0; iBlade < nBlades; iBlade++) {
        stk::mesh::Selector sel(meta_.locally_owned_part() & (*bladePartVec_[iBlade]) );
        const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);
        for (auto b: bkts) {
            for (size_t in=0; in < b->size(); in++) {            
                auto node = (*b)[in];
                double* fsiForceNode = stk::mesh::field_data(*fsiForce_, node);
                double* xyz = stk::mesh::field_data(*modelCoords, node);
                int* mapNode = stk::mesh::field_data(*bldDispMap_, node);
                double* mapInterpNode = stk::mesh::field_data(*bldDispMapInterp_, node);

                //Find the interpolated reference position first
                linInterpVec(&brFSIdata_.bld_ref_pos[(*mapNode + iStart)*6], &brFSIdata_.bld_ref_pos[(*mapNode + iStart + 1)*6], *mapInterpNode, tmpNodePos.data());

                std::vector<double> tmpForceMoment(6,0.0); // Temporarily store total force and moment as (fX, fY, fZ, mX, mY, mZ)
                //Now compute the force and moment on the interpolated reference position
                computeEffForceMoment(fsiForceNode, xyz, tmpForceMoment.data(), tmpNodePos.data());

                //Now split the force and moment on the interpolated reference position into the 'left' and 'right' nodes
                splitForceMoment(tmpForceMoment.data(), *mapInterpNode, &(brFSIdata_.bld_ld[(*mapNode + iStart)*6]), &(brFSIdata_.bld_ld[(*mapNode + iStart + 1)*6]) );
                
            }
        }


        int nPtsBlade = params_.nBRfsiPtsBlade[iBlade];
        
        //Compute the total force and moment at the hub from this blade
        std::vector<double> hubForceMoment(6,0.0);
        computeHubForceMomentForPart(hubForceMoment, brFSIdata_.hub_ref_pos, bladePartVec_[iBlade]);

        //Now compute total force and moment at the hub from the loads mapped to the 
        std::vector<double> l_hubForceMomentMapped(6,0.0);
        for (size_t i=0 ; i < nPtsBlade; i++) {
            computeEffForceMoment( &(brFSIdata_.bld_ld[(i+iStart)*6]), &(brFSIdata_.bld_ref_pos[(i+iStart)*6]), l_hubForceMomentMapped.data(), brFSIdata_.hub_ref_pos.data() );
            for(size_t j=0; j < ndim; j++) // Add the moment manually
                l_hubForceMomentMapped[3+j] += brFSIdata_.bld_ld[(i+iStart)*6+3+j];
        }
        std::vector<double> hubForceMomentMapped(6,0.0);
        stk::all_reduce_sum(bulk_.parallel(), l_hubForceMomentMapped.data(), hubForceMomentMapped.data(), 6);
        
        
        if (bulk_.parallel_rank() == turbineProc_) {
            
            std::cout << "Total force moment on the hub due to blade " << iBlade << std::endl;
            std::cout << "Force = (";
            for(size_t j=0; j < ndim; j++)
                std::cout << hubForceMoment[j] << ", ";
            std::cout << ") Moment = (" ;
            for(size_t j=0; j < ndim; j++)
                std::cout << hubForceMoment[3+j] << ", ";
            std::cout << ")" << std::endl;
            std::cout << "Total force moment on the hub from mapped load due to blade " << iBlade << std::endl;
            std::cout << "Force = (";
            for(size_t j=0; j < ndim; j++)
                std::cout << hubForceMomentMapped[j] << ", ";
            std::cout << ") Moment = (" ;
            for(size_t j=0; j < ndim; j++)
                std::cout << hubForceMomentMapped[3+j] << ", ";
            std::cout << ")" << std::endl;
        }
        
        iStart += nPtsBlade;
    }
    
}

//! Split a force and moment into the surrounding 'left' and 'right' nodes in a variationally consistent manner using 
void fsiTurbine::splitForceMoment(double *totForceMoment, double interpFac, double *leftForceMoment, double *rightForceMoment) {

    for(size_t i=0; i<6; i++) {
        leftForceMoment[i] += (1.0 - interpFac) * totForceMoment[i];
        rightForceMoment[i] += interpFac * totForceMoment[i];
    }

    
}

void fsiTurbine::computeHubForceMomentForPart(std::vector<double> & hubForceMoment, std::vector<double> & hubPos, stk::mesh::Part * part) {

    const int ndim = meta_.spatial_dimension();
    VectorFieldType* modelCoords = meta_.get_field<VectorFieldType>(stk::topology::NODE_RANK, "coordinates");

    std::vector<double> l_hubForceMoment(6,0.0);
    
    stk::mesh::Selector sel(meta_.locally_owned_part() & (*part) );
    const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);
    for (auto b: bkts) {
        for (size_t in=0; in < b->size(); in++) {            
            auto node = (*b)[in];
            double* xyz = stk::mesh::field_data(*modelCoords, node);
            double* fsiForceNode = stk::mesh::field_data(*fsiForce_, node);
            computeEffForceMoment(fsiForceNode, xyz, l_hubForceMoment.data(), hubPos.data());
        }
    }

    stk::all_reduce_sum(bulk_.parallel(), l_hubForceMoment.data(), hubForceMoment.data(), 6);
    
}
    
//! Compute the effective force and moment at the OpenFAST mesh node for a given force at the CFD surface mesh node
void fsiTurbine::computeEffForceMoment(double *forceCFD, double *xyzCFD, double *forceMomOF, double *xyzOF) {

    const int ndim=3; //I don't see this ever being used in other situations
    for(size_t j=0; j < ndim; j++) 
        forceMomOF[j] += forceCFD[j];
    forceMomOF[3] += (xyzCFD[1]-xyzOF[1])*forceCFD[2] - (xyzCFD[2]-xyzOF[2])*forceCFD[1] ;
    forceMomOF[4] += (xyzCFD[2]-xyzOF[2])*forceCFD[0] - (xyzCFD[0]-xyzOF[0])*forceCFD[2] ;
    forceMomOF[5] += (xyzCFD[0]-xyzOF[0])*forceCFD[1] - (xyzCFD[1]-xyzOF[1])*forceCFD[0] ;
    
}

//! Set sample displacement on the OpenFAST mesh before mapping to the turbine blade surface mesh
void fsiTurbine::setSampleDisplacement(double curTime) {


    
    //Turbine rotates at 12.1 rpm
    double omega=(12.1/60.0)*2.0*M_PI; //12.1 rpm
    double theta=omega*curTime;

    double sinOmegaT = std::sin(omega*curTime);
    
    //Rotate the hub first
    double hubRot = 4.0*tan(0.25*theta);
    std::vector<double> wmHubRot = {hubRot, 0.0, 0.0};
//    std::vector<double> wmHubRot = {0.0, 0.0, 0.0};    
    for (size_t i=0; i<3; i++)
        brFSIdata_.hub_def[3+i] = -wmHubRot[i];

    //For each node on the openfast blade1 mesh - compute distance from the blade root node. Apply a rotation varying as the square of the distance between 0 - 45 degrees about the [0 1 0] axis. Apply a translation displacement that produces a tip displacement of 5m
    int iStart = 0;
    int nBlades = params_.numBlades;;
    for (size_t iBlade=0; iBlade < nBlades; iBlade++) {
        std::vector<double> wmRotBlade_ref = {4.0*tan(0.25 * iBlade * 120.0 * M_PI / 180.0), 0.0, 0.0};
        std::vector<double> wmRotBlade(3,0.0);
        composeWM(wmHubRot.data(), wmRotBlade_ref.data(), wmRotBlade.data());
        
        int nPtsBlade = params_.nBRfsiPtsBlade[iBlade];        
        for (size_t i=0; i < nPtsBlade; i++) {
                
            double rDistSq = calcDistanceSquared(&(brFSIdata_.bld_ref_pos[(iStart+i)*6]), &(brFSIdata_.bld_ref_pos[(iStart)*6]) )/10000.0;
            double sinRdistSq = std::sin(rDistSq);
            //Set rotational displacement
            std::vector<double> wmRot1 = {1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0)};
            std::vector<double> wmRot(3,0.0);
            applyWMrotation(wmRotBlade.data(), wmRot1.data(), wmRot.data());
            double rot = 4.0*tan(0.25 * (45.0 * M_PI / 180.0) * sinRdistSq  * sinOmegaT ); // 4.0 * tan(phi/4.0) parameter for Wiener-Milenkovic
            for(size_t j= 0; j < 3; j++) {
                wmRot[j] *= rot;
            }
            
            std::vector<double> finalRot(3,0.0);
            composeWM(wmRot.data(), wmRotBlade.data(), finalRot.data()); //Compose with hub orientation to account for rotation of turbine

            std::vector<double> origZaxis = {0.0, 0.0, 1.0};
            std::vector<double> rotZaxis(3,0.0);
            applyWMrotation(finalRot.data(), origZaxis.data(), rotZaxis.data());

            //Finally transpose the whole thing
            for(size_t j=0; j < 3; j++)
                brFSIdata_.bld_def[(iStart+i)*6+3+j] = -finalRot[j];

            
            //Set translational displacement
            double xDisp = sinRdistSq * 15.0 * sinOmegaT;

            std::vector<double> r(3,0.0);
            for(size_t j=0; j < 3; j++)
                r[j] = brFSIdata_.bld_ref_pos[(iStart+i)*6+j] - brFSIdata_.hub_ref_pos[j];

            std::vector<double> rRot(3,0.0);

            std::vector<double> transDisp = {xDisp, xDisp, xDisp};
            std::vector<double> transDispRot(3,0.0);
            
            applyWMrotation(wmRotBlade.data(), transDisp.data(), transDispRot.data());
            
            applyWMrotation(wmHubRot.data(), r.data(), rRot.data());
            brFSIdata_.bld_def[(iStart+i)*6+0] = transDispRot[0] + rRot[0] - r[0];
            brFSIdata_.bld_def[(iStart+i)*6+1] = transDispRot[1] + rRot[1] - r[1];
            brFSIdata_.bld_def[(iStart+i)*6+2] = transDispRot[2] + rRot[2] - r[2];
            
        }
        iStart += nPtsBlade;
    }
}


//! Set reference displacement on the turbine blade surface mesh, for comparison with Sample displacement set in setSampleDisplacement
void fsiTurbine::setRefDisplacement() {

  // some hard coded parameters
  std::vector<double> bl_ref_pos = {0.0, 0.0, 131.5}; // root of blade 1
  double  tipdisp = 5.0;

  // extract the vector field type set by this function
  const int ndim = meta_.spatial_dimension();
  VectorFieldType* modelCoords = meta_.get_field<VectorFieldType>(
  stk::topology::NODE_RANK, "coordinates");
  VectorFieldType* refDisp = meta_.get_field<VectorFieldType>(
  stk::topology::NODE_RANK, "mesh_displacement_ref");

  int iBlade = 0; // for testing purposes, we are doing this only for the first blade
  stk::mesh::Selector sel(*bladePartVec_[iBlade]); // extract blade
  const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel); // extract buckets for the blade

  for (auto b: bkts) { // loop over number of buckets
    for (size_t in=0; in < b->size(); in++) { // loop over all nodes in the bucket
      auto node = (*b)[in];
      double* xyz = stk::mesh::field_data(*modelCoords, node);
      double* vecRefNode = stk::mesh::field_data(*refDisp, node); // extract

      // Compute position of current node relative to blade root
      double rDistSq = calcDistanceSquared( xyz, bl_ref_pos.data() )/10000.0;
      double sinRdistSq = std::sin(rDistSq);

      //Set translational displacement
      vecRefNode[0] = sinRdistSq * tipdisp;
      vecRefNode[1] = 0.0;
      vecRefNode[2] = 0.0;

      double rot = 4.0*tan(0.25 * (45.0 * M_PI / 180.0) * sinRdistSq); // 4.0 * tan(phi/4.0) parameter for Wiener-Milenkovic
      std::vector<double> wmRot = {0.0, rot, 0.0}; //Wiener-Milenkovic parameter

      std::vector<double> r = {xyz[0], xyz[1], 0.0}; //Wiener-Milenkovic parameter
      std::vector<double> r_rot(3,0.0);
      applyWMrotation(wmRot.data(), r.data(), r_rot.data());

      for(size_t i=0; i < ndim; i++ )
        vecRefNode[i] = vecRefNode[i] + r_rot[i] - r[i];
    }
  }
}

//! Calculate the distance between 3-dimensional vectors 'a' and 'b'
double fsiTurbine::calcDistanceSquared(double * a, double * b) {

    double dist = 0.0;
    for(size_t i=0; i < 3; i++)
        dist += (a[i]-b[i])*(a[i]-b[i]);
    return dist;
    
}

//! Map the deflections from the openfast nodes to the turbine surface CFD mesh. Will call 'computeDisplacement' for each node on the turbine surface CFD mesh.
void fsiTurbine::mapDisplacements() {
    
   //To implement this function - assume that 'bldDispMap_' field contains the lower node id of the corresponding element along the blade or the tower along with the 'bldDispMapInterp_' field that contains the non-dimensional location of the CFD surface mesh node on that element.

    // For e.g., for blade 'k' if the lower node id from 'bldDisMap_' is 'j' and the non-dimenional location from 'bldDispMapInterp_' is 'm', then the translational displacement for the CFD surface mesh is
    // (1-m) * bld_def[k][j*6+0] + m * bld_def[k][(j+1)*6+0]
    // (1-m) * bld_def[k][j*6+1] + m * bld_def[k][(j+1)*6+1]
    // (1-m) * bld_def[k][j*6+2] + m * bld_def[k][(j+1)*6+2]

    //TODO: When the turbine is rotating, displacement of the surface of the blades and hub (not the nacelle and tower) should be split into a rigid body motion due to the rotation of the turbine, yawing of the turbine and a deflection of the structure itself. The yaw rate and rotation rate will vary over time. OpenFAST always stores the displacements of all nodes with respect to the reference configuration. When using a sliding mesh interface (yaw = 0), the mesh blocks inside the rotating part of the sliding interface should be moved with rigid body motion corresponding to the rotation rate first and then a second mesh deformation procedure should be performed to apply the remaining structural deflection. Figure out how to do this.

    const int ndim = meta_.spatial_dimension();
    VectorFieldType* modelCoords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    VectorFieldType* displacement = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "mesh_displacement");
    VectorFieldType* refDisplacement = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "mesh_displacement_ref");
    VectorFieldType* meshVelocity = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "mesh_velocity");

   
    //Do the tower first
    // stk::mesh::Selector sel(*twrPart_);
    // const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);

    // for (auto b: bkts) {
    //     for (size_t in=0; in < b->size(); in++) {
    //         auto node = (*b)[in];
    //         double* oldxyz = stk::mesh::field_data(*modelCoords, node);
    //         double *dx = stk::mesh::field_data(*displacement, node);
    //         int* dispMapNode = stk::mesh::field_data(*twrDispMap_, node);
    //         double* dispMapInterpNode = stk::mesh::field_data(*twrDispMapInterp_, node);
    //         double *mVel = stk::mesh::field_data(*meshVelocity, node);
            
    //         //Find the interpolated reference position first
    //         linInterpTotDisplacement(&brFSIdata_.twr_ref_pos[(*dispMapNode)*6], &brFSIdata_.twr_ref_pos[(*dispMapNode + 1)*6], *dispMapInterpNode, tmpNodePos.data());

    //         //Now linearly interpolate the deflections to the intermediate location
    //         linInterpTotDisplacement(&brFSIdata_.twr_def[(*dispMapNode)*6], &brFSIdata_.twr_def[(*dispMapNode + 1)*6], *dispMapInterpNode, totDispNode.data());

    //         //Now transfer the interpolated displacement to the CFD mesh node
    //         computeDisplacement(totDispNode.data(), tmpNodePos.data(), dx, oldxyz);

    //         //Now linearly interpolate the velocity to the intermediate location
    //         linInterpTotVelocity(&brFSIdata_.twr_vel[(*dispMapNode)*6], &brFSIdata_.twr_vel[(*dispMapNode + 1)*6], *dispMapInterpNode, totVelNode.data());

    //         //Now transfer the interpolated translational and rotational velocity to an equivalent translational velocity on the CFD mesh node
    //         computeMeshVelocity(totVelNode.data(), totDispNode.data(), tmpNodePos.data(), mVel, oldxyz);
            
    //     }
    // }


    // Now the blades
    int nBlades = params_.numBlades;
    int iStart = 0;
    for (int iBlade=0; iBlade < nBlades; iBlade++) {
        std::cout << "Mapping blade " << iBlade << std::endl ;
        int nPtsBlade = params_.nBRfsiPtsBlade[iBlade];        
        stk::mesh::Selector sel(*bladePartVec_[iBlade]);
        const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);

        for (auto b: bkts) {
            for (size_t in=0; in < b->size(); in++) {
                auto node = (*b)[in];
                double* oldxyz = stk::mesh::field_data(*modelCoords, node);
                double *dx = stk::mesh::field_data(*displacement, node);
                int* dispMapNode = stk::mesh::field_data(*bldDispMap_, node);
                int* loadMapNode = stk::mesh::field_data(*bldLoadMap_, node);                
                double* dispMapInterpNode = stk::mesh::field_data(*bldDispMapInterp_, node);
                double *mVel = stk::mesh::field_data(*meshVelocity, node);

                std::vector<double> totDispNode(6,0.0); // Total displacement at any node in (transX, transY, transZ, rotX, rotY, rotZ)
                std::vector<double> totVelNode(6,0.0); // Total velocity at any node in (transX, transY, transZ, rotX, rotY, rotZ)
                std::vector<double> tmpNodePos(6,0.0); // Vector to temporarily store a position and orientation vector

                //Find the interpolated reference position first
                linInterpTotDisplacement(&brFSIdata_.bld_ref_pos[(*dispMapNode + iStart)*6], &brFSIdata_.bld_ref_pos[(*dispMapNode + iStart + 1)*6], *dispMapInterpNode, tmpNodePos.data());
                
                //Now linearly interpolate the deflections to the intermediate location
                linInterpTotDisplacement(&brFSIdata_.bld_def[(*dispMapNode + iStart)*6], &brFSIdata_.bld_def[(*dispMapNode + iStart + 1)*6], *dispMapInterpNode, totDispNode.data());
                
                //Now transfer the interpolated displacement to the CFD mesh node
                computeDisplacement(totDispNode.data(), tmpNodePos.data(), dx, oldxyz);

                //Now linearly interpolate the velocity to the intermediate location
                linInterpTotVelocity(&brFSIdata_.bld_vel[(*dispMapNode + iStart)*6], &brFSIdata_.bld_vel[(*dispMapNode + iStart + 1)*6], *dispMapInterpNode, totVelNode.data());

                //Now transfer the interpolated translational and rotational velocity to an equivalent translational velocity on the CFD mesh node
                computeMeshVelocity(totVelNode.data(), totDispNode.data(), tmpNodePos.data(), mVel, oldxyz);
                
            }
        }

        // double errorNorm = compute_error_norm(displacement, refDisplacement, bladePartVec_[iBlade]);

        // if (!bulk_.parallel_rank()) {
        //     std::cout << "Error in displacement for blade " << iBlade << " = " << errorNorm << std::endl ;
        // }
        iStart += nPtsBlade;
    }


    //Now the hub
    stk::mesh::Selector hubsel(*hubPart_);
    const auto& hubbkts = bulk_.get_buckets(stk::topology::NODE_RANK, hubsel);
    for (auto b: hubbkts) {
        for (size_t in=0; in < b->size(); in++) {
            auto node = (*b)[in];
            double* oldxyz = stk::mesh::field_data(*modelCoords, node);
            double *dx = stk::mesh::field_data(*displacement, node);
            double *mVel = stk::mesh::field_data(*meshVelocity, node);
            
            //Now transfer the displacement to the CFD mesh node
            computeDisplacement(brFSIdata_.hub_def.data(), brFSIdata_.hub_ref_pos.data(), dx, oldxyz);

            //Now transfer the translational and rotational velocity to an equivalent translational velocity on the CFD mesh node
            computeMeshVelocity(brFSIdata_.hub_vel.data(), brFSIdata_.hub_def.data(), brFSIdata_.hub_ref_pos.data(), mVel, oldxyz);
            
        }
    }

    //Now the nacelle
    stk::mesh::Selector nacsel(*nacellePart_);
    const auto& nacbkts = bulk_.get_buckets(stk::topology::NODE_RANK, nacsel);
    for (auto b: nacbkts) {
        for (size_t in=0; in < b->size(); in++) {
            auto node = (*b)[in];
            double* oldxyz = stk::mesh::field_data(*modelCoords, node);
            double *dx = stk::mesh::field_data(*displacement, node);
            double *mVel = stk::mesh::field_data(*meshVelocity, node);
            
            //Now transfer the displacement to the CFD mesh node
            computeDisplacement(brFSIdata_.nac_def.data(), brFSIdata_.nac_ref_pos.data(), dx, oldxyz);

            //Now transfer the translational and rotational velocity to an equivalent translational velocity on the CFD mesh node
            computeMeshVelocity(brFSIdata_.nac_vel.data(), brFSIdata_.nac_def.data(), brFSIdata_.nac_ref_pos.data(), mVel, oldxyz);
            
        }
    }
    
}

//! Linearly interpolate dispInterp = dispStart + interpFac * (dispEnd - dispStart). Special considerations for Wiener-Milenkovic parameters
void fsiTurbine::linInterpTotDisplacement(double *dispStart, double *dispEnd, double interpFac, double * dispInterp) {

    // Handle the translational displacement first
    linInterpVec(dispStart, dispEnd, interpFac, dispInterp);
    // Now deal with the rotational displacement
    linInterpRotation( &dispStart[3], &dispEnd[3], interpFac, &dispInterp[3]);
    
}

//! Linearly interpolate velInterp = velStart + interpFac * (velEnd - velStart). 
void fsiTurbine::linInterpTotVelocity(double *velStart, double *velEnd, double interpFac, double * velInterp) {

    // Handle the translational velocity first
    linInterpVec(velStart, velEnd, interpFac, velInterp);
    // Now deal with the rotational velocity
    linInterpVec(&velStart[3], &velEnd[3], interpFac, &velInterp[3]);
    
}
    
//! Linearly interpolate between 3-dimensional vectors 'a' and 'b' with interpolating factor 'interpFac'
void fsiTurbine::linInterpVec(double * a, double * b, double interpFac, double * aInterpb) {

    for (size_t i=0; i < 3; i++)
        aInterpb[i] = a[i] + interpFac * (b[i] - a[i]);
    
}
    
/* Linearly interpolate the Wiener-Milenkovic parameters between 'qStart' and 'qEnd' into 'qInterp' with an interpolating factor 'interpFac'
    see O.A.Bauchau, 2011, Flexible Multibody Dynamics p. 649, section 17.2, Algorithm 1'
*/
void fsiTurbine::linInterpRotation(double * qStart, double * qEnd, double interpFac, double * qInterp) {

    std::vector<double> intermedQ(3,0.0);
    composeWM(qStart, qEnd, intermedQ.data(), -1.0); //Remove rigid body rotation of qStart
    for(size_t i=0; i < 3; i++)
        intermedQ[i] = interpFac * intermedQ[i]; // Do the interpolation
    composeWM(qStart, intermedQ.data(), qInterp); // Add rigid body rotation of qStart back
    
}

//! Compose Wiener-Milenkovic parameters 'p' and 'q' into 'pPlusq'. If a transpose of 'p' is required, set tranposeP to '-1', else leave blank or set to '+1'
void fsiTurbine::composeWM(double * p, double * q, double * pPlusq, double transposeP, double transposeQ) {
        
    double p0 = 2.0 - 0.125*dot(p,p);
    double q0 = 2.0 - 0.125*dot(q,q);
    std::vector<double> pCrossq(3,0.0);
    cross(p, q, pCrossq.data());

    double delta1 = (4.0-p0)*(4.0-q0);
    double delta2 = p0*q0 - transposeP*dot(p,q);
    double premultFac = 0.0;
    if (delta2 < 0)
        premultFac = -4.0/(delta1 - delta2);
    else
        premultFac = 4.0/(delta1 + delta2);

    for (size_t i=0; i < 3; i++)
        pPlusq[i] = premultFac * (transposeQ * p0 * q[i] + transposeP * q0 * p[i] + transposeP * transposeQ * pCrossq[i] );
    
}
    
double fsiTurbine::dot(double * a, double * b) {

    return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
    
}

void fsiTurbine::cross(double * a, double * b, double * aCrossb) {

    aCrossb[0] = a[1]*b[2] - a[2]*b[1];
    aCrossb[1] = a[2]*b[0] - a[0]*b[2];
    aCrossb[2] = a[0]*b[1] - a[1]*b[0];
    
}
    
//! Convert one array of 6 deflections (transX, transY, transZ, wmX, wmY, wmZ) into one vector of translational displacement at a given node on the turbine surface CFD mesh.
void fsiTurbine::computeDisplacement(double *totDispNode, double * totPosOF,  double *transDispNode, double * xyzCFD) {

    //Get the relative distance between totPosOF and xyzCFD in the inertial frame
    std::vector<double> p(3,0.0); 
    for (size_t i=0; i < 3; i++)
        p[i] = xyzCFD[i] - totPosOF[i];
    //Convert 'p' vector to the local frame of reference
    std::vector<double> pLoc(3,0.0);
    applyWMrotation(&(totPosOF[3]), p.data(), pLoc.data());

    std::vector<double> pRot(3,0.0);
    applyWMrotation(&(totDispNode[3]), pLoc.data(), pRot.data(),-1); // Apply the rotation corresponding to the final orientation to bring back to inertial frame

    for (size_t i=0; i < 3; i++)
        transDispNode[i] = totDispNode[i] + pRot[i] - p[i];
    
}

//! Convert one array of 6 velocities (transX, transY, transZ, wmX, wmY, wmZ) into one vector of translational velocity at a given node on the turbine surface CFD mesh.
void fsiTurbine::computeMeshVelocity(double *totVelNode, double * totDispNode, double * totPosOF,  double *transVelNode, double * xyzCFD) {

    //Get the relative distance between totPosOF and xyzCFD in the inertial frame
    std::vector<double> p(3,0.0); 
    for (size_t i=0; i < 3; i++)
        p[i] = xyzCFD[i] - totPosOF[i];
    //Convert 'p' vector to the local frame of reference
    std::vector<double> pLoc(3,0.0);
    applyWMrotation(&(totPosOF[3]), p.data(), pLoc.data(),-1.0);
    
    std::vector<double> pRot(3,0.0);
    applyWMrotation(&(totDispNode[3]), pLoc.data(), pRot.data()); // Apply the rotation corresponding to the final orientation to bring back to inertial frame

    std::vector<double> omegaCrosspRot(3,0.0);
    cross(&(totVelNode[3]), pRot.data(), omegaCrosspRot.data());
    
    for (size_t i=0; i < 3; i++)
        transVelNode[i] = totVelNode[i] + omegaCrosspRot[i];
    
}
    
double fsiTurbine::compute_error_norm(VectorFieldType * vec, VectorFieldType * vec_ref, stk::mesh::Part * part) {

    const int ndim = meta_.spatial_dimension();
    std::vector<double> errorNorm(3,0);
    int nNodes = 0;
    
    stk::mesh::Selector sel(*part);
    const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);

    for (auto b: bkts) {
        for (size_t in=0; in < b->size(); in++) {
            auto node = (*b)[in];
            double* vecNode = stk::mesh::field_data(*vec, node);
            double* vecRefNode = stk::mesh::field_data(*vec_ref, node);
            for(size_t i=0; i < ndim; i++)
                errorNorm[i] += (vecNode[i] - vecRefNode[i])*(vecNode[i] - vecRefNode[i]);

            nNodes++;
        }
    }

    std::vector<double> g_errorNorm(3,0.0);
    stk::all_reduce_sum(bulk_.parallel(), errorNorm.data(), g_errorNorm.data(), 3);
    
    int g_nNodes = 0;
    stk::all_reduce_sum(bulk_.parallel(), &nNodes, &g_nNodes, 1);

    for (size_t i=0; i < ndim; i++)
        g_errorNorm[i] = sqrt(g_errorNorm[i]/g_nNodes);

    return g_errorNorm[0];
    
}

//! Apply a Wiener-Milenkovic rotation 'wm' to a vector 'r' into 'rRot'. To optionally transpose the rotation, set 'tranpose=-1.0'.
void fsiTurbine::applyWMrotation(double * wm, double * r, double *rRot, double transpose) {

    double wm0 = 2.0-0.125*dot(wm, wm);
    double nu = 2.0/(4.0-wm0);
    double cosPhiO2 = 0.5*wm0*nu;
    std::vector<double> wmCrossR(3,0.0);
    cross(wm, r, wmCrossR.data());
    std::vector<double> wmCrosswmCrossR(3,0.0);
    cross(wm, wmCrossR.data(), wmCrosswmCrossR.data());
    
    for(size_t i=0; i < 3; i++)
        rRot[i] = r[i] + transpose * nu * cosPhiO2 * wmCrossR[i] + 0.5 * nu * nu * wmCrosswmCrossR[i];
    
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
            std::vector<double> ptCoords(ndim, 0.0);
            for(int i=0; i < ndim; i++)
                ptCoords[i] = xyz[i];
            bool foundProj = false;
            double nDimCoord = -1.0;
            int nPtsTwr = params_.nBRfsiPtsTwr;
            if (nPtsTwr > 0) {
                for (int i=0; i < nPtsTwr-1; i++) {
                    std::vector<double> lStart = {brFSIdata_.twr_ref_pos[i*6], brFSIdata_.twr_ref_pos[i*6+1], brFSIdata_.twr_ref_pos[i*6+2]};
                    std::vector<double> lEnd = {brFSIdata_.twr_ref_pos[(i+1)*6], brFSIdata_.twr_ref_pos[(i+1)*6+1], brFSIdata_.twr_ref_pos[(i+1)*6+2]};
                    nDimCoord = projectPt2Line(ptCoords, lStart, lEnd);
                    
                    if ((nDimCoord >= 0) && (nDimCoord <= 1.0)) {
                        *dispMapInterpNode = nDimCoord;
                        *dispMapNode = i;
//                        *loadMapNode = i + std::round(nDimCoord);
                        foundProj = true;
                        break;
                    } 
                }
                
                //If no element in the OpenFAST mesh contains this node do some sanity check on the perpendicular distance between the surface mesh node and the line joining the ends of the tower
                if (!foundProj) {
                    std::vector<double> lStart = {brFSIdata_.twr_ref_pos[0], brFSIdata_.twr_ref_pos[1], brFSIdata_.twr_ref_pos[2]};
                    std::vector<double> lEnd = {brFSIdata_.twr_ref_pos[(nPtsTwr-1)*6], brFSIdata_.twr_ref_pos[(nPtsTwr-1)*6+1], brFSIdata_.twr_ref_pos[(nPtsTwr-1)*6+2]};
                    double perpDist = perpProjectDist_Pt2Line(ptCoords, lStart, lEnd);
                    if (perpDist > 0.2) {// Something's wrong if a node on the surface mesh of the tower is more than 20% of the tower length away from the tower axis. 
                        throw std::runtime_error("Can't find a projection for point (" + std::to_string(ptCoords[0]) + "," + std::to_string(ptCoords[1]) + "," + std::to_string(ptCoords[2]) + ") on the tower on turbine " + std::to_string(params_.TurbID) + ". The tower extends from " + std::to_string(lStart[0]) + "," + std::to_string(lStart[1]) + "," + std::to_string(lStart[2]) + ") to " + std::to_string(lEnd[0]) + "," + std::to_string(lEnd[1]) + "," + std::to_string(lEnd[2]) + "). Are you sure the initial position and orientation of the mesh is consistent with the input file parameters and the OpenFAST model.");
                    }                
                    if (nDimCoord < 0.0)  {
                        //Assign this node to the first point and element of the OpenFAST mesh
                        *dispMapInterpNode = 0.0;
                        *dispMapNode = 0;
//                        *loadMapNode = 0;
                    } else if (nDimCoord > 1.0) { //Assign this node to the last point and element of the OpenFAST mesh
                        *dispMapInterpNode = 1.0;
                        *dispMapNode = nPtsTwr-2;
//                        *loadMapNode = nPtsTwr-1;
                    }
                }
            }
        }
    }
    
    // Now the blades
    int nBlades = params_.numBlades;
    int iStart = 0;
    for (int iBlade=0; iBlade < nBlades; iBlade++) {        
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
                std::vector<double> ptCoords(ndim, 0.0);
                for(int i=0; i < ndim; i++)
                    ptCoords[i] = xyz[i];
                bool foundProj = false;
                double nDimCoord = -1.0;
                for (int i=0; i < nPtsBlade-1; i++) {
                    std::vector<double> lStart = {brFSIdata_.bld_ref_pos[(iStart+i)*6], brFSIdata_.bld_ref_pos[(iStart+i)*6+1], brFSIdata_.bld_ref_pos[(iStart+i)*6+2]};
                    std::vector<double> lEnd = {brFSIdata_.bld_ref_pos[(iStart+i+1)*6], brFSIdata_.bld_ref_pos[(iStart+i+1)*6+1], brFSIdata_.bld_ref_pos[(iStart+i+1)*6+2]};
                    nDimCoord = projectPt2Line(ptCoords, lStart, lEnd);
                    
                    if ((nDimCoord >= 0) && (nDimCoord <= 1.0)) {
                        foundProj = true;
                        *dispMapInterpNode = nDimCoord;
                        *dispMapNode = i;
//                        *loadMapNode = i + std::round(nDimCoord);
                        break;
                    }
                }
                
                //If no element in the OpenFAST mesh contains this node do some sanity check on the perpendicular distance between the surface mesh node and the line joining the ends of the blade
                if (!foundProj) {

                    std::vector<double> lStart = {brFSIdata_.bld_ref_pos[iStart*6], brFSIdata_.bld_ref_pos[iStart*6+1], brFSIdata_.bld_ref_pos[iStart*6+2]};
                    std::vector<double> lEnd = {brFSIdata_.bld_ref_pos[(iStart+nPtsBlade-1)*6], brFSIdata_.bld_ref_pos[(iStart+nPtsBlade-1)*6+1], brFSIdata_.bld_ref_pos[(iStart+nPtsBlade-1)*6+2]};

                    std::cout << "Can't find a projection for point (" + std::to_string(ptCoords[0]) + "," + std::to_string(ptCoords[1]) + "," + std::to_string(ptCoords[2]) + ") on blade " + std::to_string(iBlade) + " on turbine " + std::to_string(params_.TurbID) + ". The blade extends from " + std::to_string(lStart[0]) + "," + std::to_string(lStart[1]) + "," + std::to_string(lStart[2]) + ") to " + std::to_string(lEnd[0]) + "," + std::to_string(lEnd[1]) + "," + std::to_string(lEnd[2]) + ")." << std::endl ;                    
                    double perpDist = perpProjectDist_Pt2Line(ptCoords, lStart, lEnd);
                    if (perpDist > 0.2) {// Something's wrong if a node on the surface mesh of the blade is more than 20% of the blade length away from the blade axis. 
                        throw std::runtime_error("Can't find a projection for point (" + std::to_string(ptCoords[0]) + "," + std::to_string(ptCoords[1]) + "," + std::to_string(ptCoords[2]) + ") on blade " + std::to_string(iBlade) + " on turbine " + std::to_string(params_.TurbID) + ". The blade extends from " + std::to_string(lStart[0]) + "," + std::to_string(lStart[1]) + "," + std::to_string(lStart[2]) + ") to " + std::to_string(lEnd[0]) + "," + std::to_string(lEnd[1]) + "," + std::to_string(lEnd[2]) + "). Are you sure the initial position and orientation of the mesh is consistent with the input file parameters and the OpenFAST model.");
                    }
                    
                    if (nDimCoord < 0.0)  {
                        //Assign this node to the first point and element of the OpenFAST mesh
                            *dispMapInterpNode = 0.0;
                            *dispMapNode = 0;
//                            *loadMapNode = 0;
                    } else if (nDimCoord > 1.0) { //Assign this node to the last point and element of the OpenFAST mesh
                        *dispMapInterpNode = 1.0;
                        *dispMapNode = nPtsBlade-2;
//                        *loadMapNode = nPtsBlade-1;
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

    

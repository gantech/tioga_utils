#ifndef FSITURBINEN_H
#define FSITURBINEN_H

#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/CoordinateSystems.hpp"
#include "stk_mesh/base/Field.hpp"

#include <vector>
#include <string>

#include "yaml-cpp/yaml.h"

namespace tioga_nalu {

typedef stk::mesh::Field<double, stk::mesh::Cartesian> VectorFieldType;
typedef stk::mesh::Field<double> ScalarFieldType;
typedef stk::mesh::Field<int> ScalarIntFieldType;

class fsiTurbine
{

public:
    
    fsiTurbine(int iTurb, const YAML::Node &);

    virtual ~fsiTurbine();

    //! Convert pressure and viscous/turbulent stress on the turbine surface CFD mesh into a "fsiForce" field on the turbine surface CFD mesh
    virtual void computeFSIforce() ;

    //! Transfer loads from the "fsiForce" field on the turbine surface CFD mesh into point load array that gets transferred to openfast
    virtual void transferLoads();
    
    //! Transfer the deflections from the openfast nodes to the turbine surface CFD mesh. Will call 'computeDisplacement' for each node on the turbine surface CFD mesh.
    virtual void transferDisplacements();
    
    //! Convert one array of 6 deflections (transX, transY, transZ, wmX, wmY, wmZ) into one vector of translational displacement at a given node on the turbine surface CFD mesh.
    virtual void computeDisplacement();

    //! Map each node on the turbine surface CFD mesh to 
    virtual void computeMapping();
    
    //!Array of 6 deflections (transX, transY, transZ, wmX, wmY, wmZ) for each input node from OpenFAST on the tower
    std::vector<double> twr_def_;
    //!Array of 6 deflections (transX, transY, transZ, wmX, wmY, wmZ) for each input node from OpenFAST on each blade        
    std::vector<std::vector<double>> bld_def_;
    //! Deflections of the hub point
    std::vector<double> hub_def_;
    //! Deflections of the nacelle point
    std::vector<double> nac_def_;

private:

    fsiTurbine() = delete;
    fsiTurbine(const fsiTurbine&) = delete;
    
    int turbineProc_; // The MPI rank containing the OpenFAST instance of the turbine
    bool turbineInProc_; // A boolean flag to determine if the processor contains any part of the turbine

    ScalarIntFieldType * twrLoadMap_; // Maps every node on the tower surface to the closest node of the openfast tower mesh element 
    ScalarIntFieldType * twrDispMap_; // Maps every node on the tower surface to the lower node of the openfast mesh element containing the projection of the tower surface node on to the openfast mesh tower element
    int nBlades_; // Number of blades in the turbine    
    std::vector<ScalarIntFieldType *> bldLoadMap_; // Maps every node on the blade surface to the closest node of the openfast blade mesh element 
    std::vector<ScalarIntFieldType *> bldDispMap_; // Maps every node on the blade surface to the lower node of the openfast mesh element containing the projection of the blade surface node on to the openfast mesh blade element
    std::vector<ScalarFieldType *> bldDispMapInterp_; // The location of the CFD surface mesh node projected along the OpenFAST mesh element in non-dimensional [0,1] co-ordinates.

    //! Field containing the FSI force at all nodes on the turbine surface
    VectorFieldType * fsiForce_;

    //! Part name of the tower
    std::string twrPartName_;
    //! Pointer to tower part
    stk::mesh::Part * twrPart_;
    //! Part name of the hub
    std::string hubPartName_;
    //! Pointer to hub part
    stk::mesh::Part * hubPart_;
    //! Part name of the nacelle
    std::string nacellePartName_;
    //! Pointer to the nacelle part
    stk::mesh::Part * nacellePart_;
    //! Part names of the blades
    std::vector<std::string> bladePartNames_;
    //! Pointers to the blade parts
    stk::mesh::PartVector bladePartVec_;
    
    
};

}

#endif /* FSITURBINE_H */


#include "OpenfastFSI.h"
#include "fsiTurbine.h"

#include <cassert>
#include <cmath>

namespace tioga_nalu {

OpenfastFSI::OpenfastFSI(
    stk::mesh::MetaData& meta,
    stk::mesh::BulkData& bulk,
    const YAML::Node& node
) : MotionBase(meta, bulk)
{    
    load(node);
}

OpenfastFSI::~OpenfastFSI() {
    
    FAST.end() ;
    
}
    
void OpenfastFSI::read_turbine_data(int iTurb, fast::fastInputs & fi, YAML::Node turbNode) {

    //Read turbine data for a given turbine using the YAML node
    if (turbNode["turb_id"]) {
        fi.globTurbineData[iTurb].TurbID = turbNode["turb_id"].as<int>();
    } else {
        turbNode["turb_id"] = iTurb;
    }
    if (turbNode["sim_type"]) {
        if (turbNode["sim_type"].as<std::string>() == "ext-loads") {
            fi.globTurbineData[iTurb].sType = fast::EXTLOADS;
            fsiTurbineData_[iTurb] = new fsiTurbine(iTurb, turbNode, meta_, bulk_);
        } else {
            fi.globTurbineData[iTurb].sType = fast::EXTINFLOW;
        }
    } else {
        fi.globTurbineData[iTurb].sType = fast::EXTINFLOW;
    }
    if (turbNode["FAST_input_filename"]) {
        fi.globTurbineData[iTurb].FASTInputFileName = turbNode["FAST_input_filename"].as<std::string>() ;  
    } else {
        fi.globTurbineData[iTurb].FASTInputFileName = "";
    }
    if (turbNode["restart_filename"]) {
        fi.globTurbineData[iTurb].FASTRestartFileName = turbNode["restart_filename"].as<std::string>() ;
    } else {
        fi.globTurbineData[iTurb].FASTRestartFileName = "";
    }
    if ( (fi.globTurbineData[iTurb].FASTRestartFileName == "") && (fi.globTurbineData[iTurb].FASTInputFileName == "") )
        throw std::runtime_error("Both FAST_input_filename and restart_filename are empty or not specified for Turbine " + std::to_string(iTurb));
    if (turbNode["turbine_base_pos"].IsSequence() ) {
        fi.globTurbineData[iTurb].TurbineBasePos = turbNode["turbine_base_pos"].as<std::vector<float> >() ;
    } else {
        fi.globTurbineData[iTurb].TurbineBasePos = std::vector<float>(3,0.0);
    }
    if (turbNode["turbine_hub_pos"].IsSequence() ) {
        fi.globTurbineData[iTurb].TurbineHubPos = turbNode["turbine_hub_pos"].as<std::vector<double> >() ;
    } else {
        fi.globTurbineData[iTurb].TurbineHubPos =  std::vector<double>(3,0.0);
    }
    if (turbNode["num_force_pts_blade"]) {
        fi.globTurbineData[iTurb].numForcePtsBlade = turbNode["num_force_pts_blade"].as<int>();
    } else {
        fi.globTurbineData[iTurb].numForcePtsBlade = 0;
    }
    if (turbNode["num_force_pts_tower"]) {
        fi.globTurbineData[iTurb].numForcePtsTwr = turbNode["num_force_pts_tower"].as<int>();
    } else {
        fi.globTurbineData[iTurb].numForcePtsTwr = 0;
    }
    if (turbNode["nacelle_cd"]) {
        fi.globTurbineData[iTurb].nacelle_cd = turbNode["nacelle_cd"].as<float>();
    } else {
        fi.globTurbineData[iTurb].nacelle_cd = 0.0;
    }
    if (turbNode["nacelle_area"]) {
        fi.globTurbineData[iTurb].nacelle_area = turbNode["nacelle_area"].as<float>();
    } else {
        fi.globTurbineData[iTurb].nacelle_area = 0.0;
    }
    if (turbNode["air_density"]) {
        fi.globTurbineData[iTurb].air_density = turbNode["air_density"].as<float>();
    } else {
        fi.globTurbineData[iTurb].air_density = 0.0;
    }
}

void OpenfastFSI::load(const YAML::Node& node)
{

    fi.comm = MPI_COMM_WORLD;
    
    fi.nTurbinesGlob = node["n_turbines_glob"].as<int>();
    
    if (fi.nTurbinesGlob > 0) {

        fsiTurbineData_.resize(fi.nTurbinesGlob);
        
        if(node["dry_run"]) {
            fi.dryRun = node["dry_run"].as<bool>();
        } 
        
        if(node["debug"]) {
            fi.debug = node["debug"].as<bool>();
        }
        
        if(node["sim_start"]) {
            if (node["sim_start"].as<std::string>() == "init") {
                fi.simStart = fast::init;
            } else if(node["sim_start"].as<std::string>() == "trueRestart") {
                fi.simStart = fast::trueRestart;
            } else if(node["sim_start"].as<std::string>() == "restartDriverInitFAST") {
                fi.simStart = fast::restartDriverInitFAST;
            } else {
                throw std::runtime_error("sim_start is not well defined in the input file");
            }
        }
        
        if(node["t_start"]) { 
            fi.tStart = node["t_start"].as<double>();
        } else {
            throw std::runtime_error("t_start is missing in the input file");
        }
        if(node["n_checkpoint"]) {
            fi.nEveryCheckPoint = node["n_checkpoint"].as<int>();
        } else {
            throw std::runtime_error("n_checkpoint is missing in the input file");
        }
        if (node["dt_FAST"]) {
            fi.dtFAST = node["dt_FAST"].as<double>();
        } else {
            throw std::runtime_error("dt_FAST is missing in the input file");
        }
        if (node["n_substeps"]) {
            fi.nSubsteps = node["n_substeps"].as<int>();
        } else {
            std::cout << "n_substeps is missing in the input file. Proceeding with assumption n_substeps=1." << std::endl ;
            fi.nSubsteps = 1;
        }
        if (node["t_max"]) {
            fi.tMax = node["t_max"].as<double>(); // t_max is the total duration to which you want to run FAST. This should be the same or greater than the max time given in the FAST fst file. Choose this carefully as FAST writes the output file only at this point if you choose the binary file output.
        } else {
            fi.tMax = 0;
        }
        
        fi.globTurbineData.resize(fi.nTurbinesGlob);
        for (int iTurb=0; iTurb < fi.nTurbinesGlob; iTurb++) {
            if (node["Turbine" + std::to_string(iTurb)]) {
                read_turbine_data(iTurb, fi, node["Turbine" + std::to_string(iTurb)] );
            } else {
                throw std::runtime_error("Node for Turbine" + std::to_string(iTurb) + " not present in input file or I cannot read it");
            }
        }
        
    } else {
        throw std::runtime_error("Number of turbines <= 0 ");
    }
    
    
    const auto& fparts = node["mesh_parts"];
    if (fparts.Type() == YAML::NodeType::Scalar) {
        partNames_.push_back(fparts.as<std::string>());
    } else {
        partNames_ = fparts.as<std::vector<std::string>>();
    }
    
    FAST.setInputs(fi);

}

void OpenfastFSI::setup() {

    int nTurbinesGlob = FAST.get_nTurbinesGlob();
    for (int i=0; i < nTurbinesGlob; i++) {
        if(fsiTurbineData_[i] != NULL) // This may not be a turbine intended for blade-resolved simulation
            fsiTurbineData_[i]->setup();
    }
    
}

void OpenfastFSI::initialize(double initial_time)
{
    
    FAST.allocateTurbinesToProcsSimple();
    std::cout << "Proc " << bulk_.parallel_rank() << " Called allocate" << std::endl ;
    FAST.init();
    //TODO: Check here on the processor containing the turbine that the number of blades on the turbine is the same as the number of blade parts specified in the input file.

    //TODO: In the documentation, mention that the CFD mesh must always be created for the turbine in the reference position defined in OpenFAST, i.e. with blade1 pointing up and the other blades following it in order as the turbine rotates clockwise facing downwind. If the mesh is not created this way, the mapping won't work. Any non-zero initial azimuth and/or initial yaw must be only specified in the OpenFAST input file and the mesh will automatically be deformed after calling solution0. Requiring the initial CFD mesh to be in the reference configuration may not always work if the mesh domain and initial yaw setting does not align with the reference configuration. May be this is isn't an issue because unlike AeroDyn, ExtLoads does not create the initial mesh independent of the ElastoDyn/BeamDyn. May be ExtLoads already has the correct yaw and azimuth setting from OpenFAST after the init call. In this case, the CFD mesh must start in the correct azimuth and yaw configuration. In which case, the initial yaw and azimuth must be obtained from OpenFAST and the mesh around the turbine must be deformed through rigid body motion first before starting any mapping.

    int nTurbinesGlob = FAST.get_nTurbinesGlob();
    for (int i=0; i < nTurbinesGlob; i++) {
        if(fsiTurbineData_[i] != NULL) {// This may not be a turbine intended for blade-resolved simulation
            int turbProc = FAST.getProc(i);
            fsiTurbineData_[i]->setProc(turbProc);
            if (bulk_.parallel_rank() == turbProc) {
                FAST.get_turbineParams(i, fsiTurbineData_[i]->params_);
            }
            bcast_turbine_params(i);
            fsiTurbineData_[i]->initialize();
         }
    }
    
    compute_mapping();
    std::cout << "Proc " << bulk_.parallel_rank() << " Called computeMapping" << std::endl ;
    send_loads();
    std::cout << "Proc " << bulk_.parallel_rank() << " Called send loads" << std::endl ;    
    FAST.solution0();
    deform_mesh(initial_time);
    
}

void OpenfastFSI::bcast_turbine_params(int iTurb) {

    std::vector<int> tIntParams(7,0); //Assumes a max number of blades of 3
    std::vector<double> tDoubleParams(6,0.0);
    int turbProc = fsiTurbineData_[iTurb]->getProc();
    if (bulk_.parallel_rank() == turbProc) {
        tIntParams[0] = fsiTurbineData_[iTurb]->params_.TurbID;
        tIntParams[1] = fsiTurbineData_[iTurb]->params_.numBlades;
        tIntParams[2] = fsiTurbineData_[iTurb]->params_.nTotBRfsiPtsBlade;
        tIntParams[3] = fsiTurbineData_[iTurb]->params_.nBRfsiPtsTwr;
        for (size_t i=0; i < fsiTurbineData_[iTurb]->params_.numBlades; i++) {
            tIntParams[4+i] = fsiTurbineData_[iTurb]->params_.nBRfsiPtsBlade[i];
        }

        for(size_t i=0; i < 3; i++) {
            tDoubleParams[i] = fsiTurbineData_[iTurb]->params_.TurbineBasePos[i];
            tDoubleParams[3+i] = fsiTurbineData_[iTurb]->params_.TurbineHubPos[i];
        }

    }
    int iError = MPI_Bcast(tIntParams.data(), 7, MPI_INT, turbProc, bulk_.parallel());
    iError = MPI_Bcast(tDoubleParams.data(), 6, MPI_DOUBLE, turbProc, bulk_.parallel());


    if (bulk_.parallel_rank() != turbProc) {
        fsiTurbineData_[iTurb]->params_.TurbID = tIntParams[0];
        fsiTurbineData_[iTurb]->params_.numBlades = tIntParams[1] ;
        fsiTurbineData_[iTurb]->params_.nTotBRfsiPtsBlade = tIntParams[2] ;
        fsiTurbineData_[iTurb]->params_.nBRfsiPtsTwr = tIntParams[3] ;
        fsiTurbineData_[iTurb]->params_.nBRfsiPtsBlade.resize(fsiTurbineData_[iTurb]->params_.numBlades);
        for (size_t i=0; i < fsiTurbineData_[iTurb]->params_.numBlades; i++) {
            fsiTurbineData_[iTurb]->params_.nBRfsiPtsBlade[i] = tIntParams[4+i] ;
        }

        fsiTurbineData_[iTurb]->params_.TurbineBasePos.resize(3);
        fsiTurbineData_[iTurb]->params_.TurbineHubPos.resize(3);
        for(size_t i=0; i < 3; i++) {
            fsiTurbineData_[iTurb]->params_.TurbineBasePos[i] = tDoubleParams[i];
            fsiTurbineData_[iTurb]->params_.TurbineHubPos[i] = tDoubleParams[3+i];
        }
        
    }

    
        

    
    
}

void OpenfastFSI::compute_mapping() {

    int nTurbinesGlob = FAST.get_nTurbinesGlob();
    for (int i=0; i < nTurbinesGlob; i++) {
        if(fsiTurbineData_[i] != NULL) {// This may not be a turbine intended for blade-resolved simulation
            int turbProc = fsiTurbineData_[i]->getProc();
            if (bulk_.parallel_rank() == turbProc) {
                FAST.getTowerRefPositions(fsiTurbineData_[i]->brFSIdata_.twr_ref_pos, i);
                FAST.getBladeRefPositions(fsiTurbineData_[i]->brFSIdata_.bld_ref_pos, i);
                FAST.getHubRefPosition(fsiTurbineData_[i]->brFSIdata_.hub_ref_pos, i);
                FAST.getNacelleRefPosition(fsiTurbineData_[i]->brFSIdata_.nac_ref_pos, i);
            }
            int iError = MPI_Bcast(fsiTurbineData_[i]->brFSIdata_.twr_ref_pos.data(), (fsiTurbineData_[i]->params_.nBRfsiPtsTwr)*6, MPI_DOUBLE, turbProc, bulk_.parallel());
            int nTotBldNodes = fsiTurbineData_[i]->params_.nTotBRfsiPtsBlade;
            iError = MPI_Bcast(fsiTurbineData_[i]->brFSIdata_.bld_ref_pos.data(), nTotBldNodes*6, MPI_DOUBLE, turbProc, bulk_.parallel());
            iError = MPI_Bcast(fsiTurbineData_[i]->brFSIdata_.hub_ref_pos.data(), 6, MPI_DOUBLE, turbProc, bulk_.parallel());
            iError = MPI_Bcast(fsiTurbineData_[i]->brFSIdata_.nac_ref_pos.data(), 6, MPI_DOUBLE, turbProc, bulk_.parallel());

            fsiTurbineData_[i]->computeMapping();
        }
    }
    
}

void OpenfastFSI::execute(double current_time)
{

    deform_mesh(current_time);
    //In Nalu - the CFD time step should ideally be performed here. 
    send_loads();
    FAST.update_states_driver_time_step();
    FAST.advance_to_next_driver_time_step();
    
}
  

void OpenfastFSI::send_loads() {

    int nTurbinesGlob = FAST.get_nTurbinesGlob();
    for (int i=0; i < nTurbinesGlob; i++) {
        if(fsiTurbineData_[i] != NULL) {// This may not be a turbine intended for blade-resolved simulation
            int turbProc = fsiTurbineData_[i]->getProc();
            fsiTurbineData_[i]->mapLoads();

            int nTotBldNodes = fsiTurbineData_[i]->params_.nTotBRfsiPtsBlade;            
            if (bulk_.parallel_rank() == turbProc) {
                int iError = MPI_Reduce(MPI_IN_PLACE, fsiTurbineData_[i]->brFSIdata_.twr_ld.data(), (fsiTurbineData_[i]->params_.nBRfsiPtsTwr)*6, MPI_DOUBLE, MPI_SUM, turbProc, bulk_.parallel());
                iError = MPI_Reduce(MPI_IN_PLACE, fsiTurbineData_[i]->brFSIdata_.bld_ld.data(), nTotBldNodes*6, MPI_DOUBLE, MPI_SUM, turbProc, bulk_.parallel());
                FAST.setTowerForces(fsiTurbineData_[i]->brFSIdata_.twr_ld, i, fast::np1);
                FAST.setBladeForces(fsiTurbineData_[i]->brFSIdata_.bld_ld, i, fast::np1);
            } else {
                int iError = MPI_Reduce(fsiTurbineData_[i]->brFSIdata_.twr_ld.data(), NULL, (fsiTurbineData_[i]->params_.nBRfsiPtsTwr)*6, MPI_DOUBLE, MPI_SUM, turbProc, bulk_.parallel());
                iError = MPI_Reduce(fsiTurbineData_[i]->brFSIdata_.bld_ld.data(), NULL, (nTotBldNodes)*6, MPI_DOUBLE, MPI_SUM, turbProc, bulk_.parallel());
            }
        }
    }
    
}

void OpenfastFSI::get_displacements() {

    int nTurbinesGlob = FAST.get_nTurbinesGlob();
    for (int i=0; i < nTurbinesGlob; i++) {
        if(fsiTurbineData_[i] != NULL) {// This may not be a turbine intended for blade-resolved simulation
            int turbProc = fsiTurbineData_[i]->getProc();
            if (bulk_.parallel_rank() == turbProc) {            
                FAST.getTowerDisplacements(fsiTurbineData_[i]->brFSIdata_.twr_def, fsiTurbineData_[i]->brFSIdata_.twr_vel, i, fast::np1);
                FAST.getBladeDisplacements(fsiTurbineData_[i]->brFSIdata_.bld_def, fsiTurbineData_[i]->brFSIdata_.bld_vel, i, fast::np1);
                FAST.getHubDisplacement(fsiTurbineData_[i]->brFSIdata_.hub_def, fsiTurbineData_[i]->brFSIdata_.hub_vel, i, fast::np1);
                FAST.getNacelleDisplacement(fsiTurbineData_[i]->brFSIdata_.nac_def, fsiTurbineData_[i]->brFSIdata_.nac_vel, i, fast::np1);
            }
            int iError = MPI_Bcast(fsiTurbineData_[i]->brFSIdata_.twr_ref_pos.data(), (fsiTurbineData_[i]->params_.nBRfsiPtsTwr)*6, MPI_DOUBLE, turbProc, bulk_.parallel());
            int nTotBldNodes = fsiTurbineData_[i]->params_.nTotBRfsiPtsBlade;
            iError = MPI_Bcast(fsiTurbineData_[i]->brFSIdata_.bld_ref_pos.data(), nTotBldNodes*6, MPI_DOUBLE, turbProc, bulk_.parallel());
            iError = MPI_Bcast(fsiTurbineData_[i]->brFSIdata_.hub_ref_pos.data(), 6, MPI_DOUBLE, turbProc, bulk_.parallel());
            iError = MPI_Bcast(fsiTurbineData_[i]->brFSIdata_.nac_ref_pos.data(), 6, MPI_DOUBLE, turbProc, bulk_.parallel());
            
            
        }
    }
    
}

void OpenfastFSI::deform_mesh(double current_time)
{

    std::cout << "Getting displacements from OpenFAST" << std::endl ;
    get_displacements(); // Get displacements from the OpenFAST - C++ API

    int nTurbinesGlob = FAST.get_nTurbinesGlob();
    for (int i=0; i < nTurbinesGlob; i++) {
        if(fsiTurbineData_[i] != NULL) {// This may not be a turbine intended for blade-resolved simulation {
            fsiTurbineData_[i]->setSampleDisplacement();
            fsiTurbineData_[i]->setRefDisplacement();
            fsiTurbineData_[i]->mapDisplacements();
        }
    }
    
    const int ndim = meta_.spatial_dimension();
    VectorFieldType* modelCoords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "coordinates");
    VectorFieldType* currCoords = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "current_coordinates");
    VectorFieldType* displacement = meta_.get_field<VectorFieldType>(
        stk::topology::NODE_RANK, "mesh_displacement");

    stk::mesh::Selector sel = stk::mesh::selectUnion(partVec_);
    const auto& bkts = bulk_.get_buckets(stk::topology::NODE_RANK, sel);

    // double mag = 0.0;
    // for (int d=0; d < ndim; d++)
    //     mag += axis_[d] * axis_[d];
    // mag = std::sqrt(mag);
    // const double angle = omega_ * current_time;
    // const double cosang = std::cos(0.5*angle);
    // const double sinang = std::sin(0.5*angle);
    // const double q0 = cosang;
    // const double q1 = sinang * axis_[0] / mag;
    // const double q2 = sinang * axis_[1] / mag;
    // const double q3 = sinang * axis_[2] / mag;

    // for (auto b: bkts) {
    //     for (size_t in=0; in < b->size(); in++) {
    //         auto node = (*b)[in];
    //         double* oldxyz = stk::mesh::field_data(*modelCoords, node);
    //         double* xyz = stk::mesh::field_data(*currCoords, node);
    //         double* dx = stk::mesh::field_data(*displacement, node);

    //         const double cx = oldxyz[0] - origin_[0];
    //         const double cy = oldxyz[1] - origin_[1];
    //         const double cz = oldxyz[2] - origin_[2];

    //         xyz[0] = (q0*q0 + q1*q1 - q2*q2 - q3*q3) * cx +
    //             2.0 * (q1*q2 - q0*q3) * cy +
    //             2.0 * (q0*q2 + q1*q3) * cz + origin_[0];

    //         xyz[1] = 2.0 * (q1*q2 + q0*q3) * cx +
    //             (q0*q0 - q1*q1 + q2*q2 - q3*q3) * cy +
    //             2.0 * (q2*q3 - q0*q1) * cz + origin_[1];

    //         xyz[2] = 2.0 * (q1*q3 - q0*q2) * cx +
    //             2.0 * (q0*q1 + q2*q3) * cy +
    //             (q0*q0 - q1*q1 - q2*q2 + q3*q3) * cz + origin_[2];

    //         for (int d=0; d < ndim; d++)
    //             dx[d] = xyz[d] - oldxyz[d];
    //     }
    // }
}

} // tioga_nalu

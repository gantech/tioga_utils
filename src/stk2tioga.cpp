

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/BroadcastArg.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/perf_util.hpp>

#include <stk_mesh/base/FindRestriction.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/Stencils.hpp>
#include <stk_mesh/base/TopologyDimensions.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <Ionit_Initializer.h>


#include <cmath>

#include "TiogaSTKIface.h"
#include "MeshMotion.h"
#include "MockRealm.h"
#include "MeshGeometry.h"
#include "Timer.h"
#include "tioga.h"

typedef stk::mesh::Field<double, stk::mesh::Cartesian> VectorFieldType;
typedef stk::mesh::Field<double> ScalarFieldType;
typedef stk::mesh::Field<int> ScalarIntFieldType;


const double pi = std::acos(-1.0);

void print_memory_diag(const stk::mesh::BulkData& bulk)
{
    const double factor = 1024.0;
    size_t curr_max, curr_min, curr_avg;
    stk::get_current_memory_usage_across_processors(
        bulk.parallel(), curr_max, curr_min, curr_avg);

    if (bulk.parallel_rank() == 0)
        std::cout << "Memory usage (KB): Avg. = "
                  << 1.0 * curr_avg / factor << "; Min. = "
                  << 1.0 * curr_min / factor << "; Max. = "
                  << 1.0 * curr_max / factor << std::endl;
}

void print_hwm_memory_diag(const stk::mesh::BulkData& bulk)
{
    std::ios::fmtflags cflags(std::cout.flags());
    const double factor = 1024.0 * 1024;
    size_t curr_max, curr_min, curr_avg;
    stk::get_memory_high_water_mark_across_processors(
        bulk.parallel(), curr_max, curr_min, curr_avg);

    if (bulk.parallel_rank() == 0)
      std::cout << "Memory high-water mark: Avg. = " << std::setw(6)
                << std::fixed << std::setprecision(4) << curr_avg / factor
                << " MB; Min. = " << std::setw(6) << std::fixed
                << std::setprecision(4) << curr_min / factor
                << " MB; Max. = " << std::setw(6) << std::fixed
                << std::setprecision(4) << curr_max / factor << " MB"
                << std::endl;
    std::cout.flags(cflags);
}

void tag_procs(stk::mesh::MetaData& meta, stk::mesh::BulkData& bulk)
{
  int iproc = bulk.parallel_rank();
  ScalarFieldType *ipnode = meta.get_field<ScalarFieldType>
    (stk::topology::NODE_RANK, "pid_node");
  ScalarFieldType *ipelem = meta.get_field<ScalarFieldType>
    (stk::topology::ELEM_RANK, "pid_elem");

  stk::mesh::Selector msel = meta.locally_owned_part();
  const stk::mesh::BucketVector& nbkts = bulk.get_buckets(
    stk::topology::NODE_RANK, msel);

  for (auto b: nbkts) {
    double* ip = stk::mesh::field_data(*ipnode, *b);
    for(size_t in=0; in < b->size(); in++) {
      ip[in] = iproc;
    }
  }

  const stk::mesh::BucketVector& ebkts = bulk.get_buckets(
    stk::topology::ELEM_RANK, msel);

  for (auto b: ebkts) {
    double* ip = stk::mesh::field_data(*ipelem, *b);
    for(size_t in=0; in < b->size(); in++) {
      ip[in] = iproc;
    }
  }

  std::vector<const stk::mesh::FieldBase*> fieldVec{ipnode, ipelem};
  stk::mesh::communicate_field_data(bulk, fieldVec);
}

void move_mesh(stk::mesh::MetaData& meta, stk::mesh::BulkData& bulk)
{
  double omega = 4.0 * 6.81041368647038;
  VectorFieldType* coords = meta.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, "coordinates");

  stk::mesh::PartVector pvec(6);
  for (int i=0; i < 6; i++) {
    std::string pname = "Unspecified-" + std::to_string(i+2) + "-HEX";
    pvec[i] = meta.get_part(pname);
  }

  stk::mesh::Selector mselect = stk::mesh::selectUnion(pvec);
  stk::mesh::BucketVector bkts = bulk.get_buckets(
    stk::topology::NODE_RANK, mselect);
  for (auto b: bkts) {
    for (size_t in=0; in < b->size(); in++) {
      double *pts = stk::mesh::field_data(*coords, (*b)[in]);
      double xold = pts[0];
      double zold = pts[2];
      pts[0] = xold * std::cos(omega) + zold * std::sin(omega);
      pts[2] = - xold * std::sin(omega) + zold * std::cos(omega);
    }
  }
}

size_t init_write_mesh(
    const YAML::Node& inpfile,
    stk::mesh::MetaData& meta,
    stk::mesh::BulkData& bulk,
    stk::io::StkMeshIoBroker& stkio,
    const double time=0.0)
{
  bool do_write = true;
  if (inpfile["write_outputs"])
    do_write = inpfile["write_outputs"].as<bool>();
  if (!do_write) return -1;

  bool has_motion = false;
  if (inpfile["motion_info"])
    has_motion = true;

  tag_procs(meta, bulk);

  std::string out_mesh = inpfile["output_mesh"].as<std::string>();
  size_t fh = stkio.create_output_mesh(out_mesh, stk::io::WRITE_RESTART);
  
  if (bulk.parallel_rank() == 0)
      std::cout << "Writing output file: " << out_mesh << std::endl;

  
  if (has_motion) {
      VectorFieldType* mesh_disp = meta.get_field<VectorFieldType>(
          stk::topology::NODE_RANK, "mesh_displacement");
      stkio.add_field(fh, *mesh_disp);
      VectorFieldType* mesh_disp_ref = meta.get_field<VectorFieldType>(
          stk::topology::NODE_RANK, "mesh_displacement_ref");
      stkio.add_field(fh, *mesh_disp_ref);

      VectorFieldType* mesh_vel = meta.get_field<VectorFieldType>(
          stk::topology::NODE_RANK, "mesh_velocity");
      stkio.add_field(fh, *mesh_vel);
      VectorFieldType* mesh_vel_ref = meta.get_field<VectorFieldType>(
          stk::topology::NODE_RANK, "mesh_velocity_ref");
      stkio.add_field(fh, *mesh_vel_ref);
      
  }

  ScalarIntFieldType* twrLoadMap = meta.get_field<ScalarIntFieldType>(
      stk::topology::NODE_RANK, "twr_load_map");
  ScalarIntFieldType* bldLoadMap = meta.get_field<ScalarIntFieldType>(
      stk::topology::NODE_RANK, "bld_load_map");
  ScalarIntFieldType* twrDispMap = meta.get_field<ScalarIntFieldType>(
      stk::topology::NODE_RANK, "twr_disp_map");
  ScalarIntFieldType* bldDispMap = meta.get_field<ScalarIntFieldType>(
      stk::topology::NODE_RANK, "bld_disp_map");
  ScalarFieldType* twrDispMapInterp = meta.get_field<ScalarFieldType>(
      stk::topology::NODE_RANK, "twr_disp_map_interp");
  ScalarFieldType* bldDispMapInterp = meta.get_field<ScalarFieldType>(
      stk::topology::NODE_RANK, "bld_disp_map_interp");
  stkio.add_field(fh, *twrLoadMap);
  stkio.add_field(fh, *bldLoadMap);
  stkio.add_field(fh, *twrDispMap);
  stkio.add_field(fh, *bldDispMap);
  stkio.add_field(fh, *twrDispMapInterp);
  stkio.add_field(fh, *bldDispMapInterp);

  ScalarFieldType* dual_nodal_volume = meta.get_field<ScalarFieldType>(
      stk::topology::NODE_RANK, "dual_nodal_volume");
  stkio.add_field(fh, *dual_nodal_volume);
  ScalarFieldType* ref_dual_nodal_volume = meta.get_field<ScalarFieldType>(
      stk::topology::NODE_RANK, "ref_dual_nodal_volume");
  stkio.add_field(fh, *ref_dual_nodal_volume);
  VectorFieldType* currCoords = meta.get_field<VectorFieldType>(
      stk::topology::NODE_RANK, "current_coordinates");
  stkio.add_field(fh, *currCoords);
  
  return fh;
}

void write_mesh(
    size_t fh,
    stk::io::StkMeshIoBroker& stkio,
    double time)
{
    stkio.begin_output_step(fh, time);
    stkio.write_defined_output_fields(fh);
    stkio.end_output_step(fh);
    
}


int main(int argc, char** argv)
{
  stk::ParallelMachine comm = stk::parallel_machine_init(&argc, &argv);

  {
      auto timerTotal = tioga_nalu::get_timer("stk2tioga::zz_total_time");
      int iproc = stk::parallel_machine_rank(comm);
      int nproc = stk::parallel_machine_size(comm);
      stk::mesh::MetaData meta;
      stk::mesh::BulkData bulk(meta, comm, stk::mesh::BulkData::NO_AUTO_AURA);

      std::string yaml_filename;
      if (argc == 2) {
          yaml_filename = argv[1];
      } else {
          throw std::runtime_error("Need input file");
      }
      YAML::Node inpfile = YAML::LoadFile(yaml_filename);

      stk::io::StkMeshIoBroker stkio(comm);

      if ((nproc > 1) && inpfile["decomposition_method"]) {
          auto decomp_method = inpfile["decomposition_method"].as<std::string>();
          stkio.property_add(Ioss::Property("DECOMPOSITION_METHOD", decomp_method));
      }

      if (iproc == 0)
          std::cout << "Preparing mesh meta data... " << std::endl;
      std::string inp_mesh = inpfile["input_mesh"].as<std::string>();
      {
          auto timeMon1 = tioga_nalu::get_timer("stk2tioga::init_meta_data");
          stkio.add_mesh_database(inp_mesh, stk::io::READ_MESH);
          stkio.set_bulk_data(bulk);
          stkio.create_input_mesh();
          stkio.add_all_mesh_fields_as_input_fields();
      }

      sierra::nalu::Simulation* sim_ = create_simulation();
      sierra::nalu::Realm& realm_ = create_mock_realm(meta, bulk, sim_);
      if (has_motion) realm_.solutionOptions_->meshMotion_ = true;
      tioga_nalu::MeshGeometry *meshGeometry_ = new tioga_nalu::MeshGeometry(&realm_);
      meshGeometry_->setup();
      
      bool has_motion = false;
      std::unique_ptr<tioga_nalu::MeshMotion> mesh_motion;
      if (inpfile["motion_info"]) {
          has_motion = true;
      }
      
      std::string coords_name = "coordinates";
      if (has_motion) {
          mesh_motion.reset(
              new tioga_nalu::MeshMotion(meta, bulk, inpfile["motion_info"]));
          coords_name = "current_coordinates";
      }
      
      const YAML::Node& oset_info = inpfile["overset_info"];

      if (iproc == 0)
          std::cout << "Calling TIOGA setup... " << std::endl;
      if (has_motion) mesh_motion->setup();

      ScalarFieldType& ipnode = meta.declare_field<ScalarFieldType>
          (stk::topology::NODE_RANK, "pid_node");
      ScalarFieldType& ipelem = meta.declare_field<ScalarFieldType>
          (stk::topology::ELEM_RANK, "pid_elem");
      stk::mesh::put_field(ipnode, meta.universal_part());
      stk::mesh::put_field(ipelem, meta.universal_part());
      
      if (iproc == 0)
          std::cout << "Loading mesh... " << std::endl;
      {
          auto timeMon2 = tioga_nalu::get_timer("stk2tioga::populate_bulk_data");
          stkio.populate_bulk_data();
      }

      if (iproc == 0)
          std::cout << "Initializing TIOGA... " << std::endl;
      if (has_motion) mesh_motion->initialize();

      meshGeometry_->initialize();
      
      size_t ofileID = init_write_mesh(inpfile, meta, bulk, stkio, 0.0);
      write_mesh(ofileID, stkio, 0.0);
      
      if (has_motion) {
          int nsteps = mesh_motion->num_steps();
          if (iproc == 0)
              std::cout << "Execution mesh motion for num steps = "
                        << nsteps << std::endl << std::endl;

          for (int nt =0; nt < nsteps; nt++) {
              mesh_motion->execute(nt);
              meshGeometry_->execute();
              if (iproc == 0)
                  std::cout << "--------------------------------------------------\n"
                            << "Time step/time = " << (nt + 1)
                            << "; " << mesh_motion->current_time()
                            << " s" << std::endl;
              print_memory_diag(bulk);

              {
                  auto timeMon3 = tioga_nalu::get_timer("stk2tioga::write_mesh");
                  double curr_time = has_motion ? mesh_motion->current_time() : 0.0;
                  write_mesh(ofileID, stkio, curr_time);
              }
              
          }

          
      }

      stk::parallel_machine_barrier(bulk.parallel());

      bool dump_partitions = false;
      if (inpfile["dump_tioga_partitions"])
          dump_partitions = inpfile["dump_tioga_partitions"].as<bool>();
      if (dump_partitions) {
          auto timeMon4 = tioga_nalu::get_timer("stk2tioga::dump_tioga_partitions");
          if (iproc == 0)
              std::cout << "Dumping tioga partitions... " << std::endl;
      }
      stk::parallel_machine_barrier(bulk.parallel());
      print_hwm_memory_diag(bulk);
  }
  Teuchos::TimeMonitor::summarize(
      std::cout, false, true, false, Teuchos::Union);
  stk::parallel_machine_finalize();
  return 0;
}

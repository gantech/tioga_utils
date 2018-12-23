#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include "Simulation.h"
#include "Realm.h"
#include "SolutionOptions.h"


const std::string naluDefaultInputs =
  "Simulations:                                                            \n"
  "  - name: sim1                                                          \n"
  "    time_integrator: ti_1                                               \n"
  "    optimizer: opt1                                                     \n"
  "                                                                        \n"
  "linear_solvers:                                                         \n"
  "                                                                        \n"
  "  - name: solve_scalar                                                  \n"
  "    type: tpetra                                                        \n"
  "    method: gmres                                                       \n"
  "    preconditioner: sgs                                                 \n"
  "    tolerance: 1e-5                                                     \n"
  "    max_iterations: 50                                                  \n"
  "    kspace: 50                                                          \n"
  "    output_level: 0                                                     \n"
  "                                                                        \n"
  "  - name: solve_cont                                                    \n"
  "    type: tpetra                                                        \n"
  "    method: gmres                                                       \n"
  "    preconditioner: muelu                                               \n"
  "    tolerance: 1e-5                                                     \n"
  "    max_iterations: 50                                                  \n"
  "    kspace: 50                                                          \n"
  "    output_level: 0                                                     \n"
  "    recompute_preconditioner: no                                        \n"
  "    muelu_xml_file_name: milestone.xml                                  \n"
  "                                                                        \n"
  "Time_Integrators:                                                       \n"
  "  - StandardTimeIntegrator:                                             \n"
  "      name: ti_1                                                        \n"
  "      start_time: 0                                                     \n"
  "      time_step: 0.1                                                    \n"
  "      termination_time: 1.5                                             \n"
  "      time_stepping_type: adaptive                                      \n"
  "      time_step_count: 0                                                \n"
  "      second_order_accuracy: yes                                        \n"
  "                                                                        \n"
  "      realms: []                                                        \n"
  "                                                                        \n"
  ;

const std::string realmDefaultSettings =
  "- name: unitTestRealm                                                  \n"
  "  use_edges: no                                                        \n"
  "                                                                       \n"
  "  equation_systems:                                                    \n"
  "    name: theEqSys                                                     \n"
  "    max_iterations: 2                                                  \n"
  "                                                                       \n"
  "    solver_system_specification:                                       \n"
  "      temperature: solve_scalar                                        \n"
  "                                                                       \n"
  "    systems:                                                           \n"
  "      - HeatConduction:                                                \n"
  "          name: myHC                                                   \n"
  "          max_iterations: 1                                            \n"
  "          convergence_tolerance: 1e-5                                  \n"
  "                                                                       \n"
  "  time_step_control:                                                   \n"
  "    target_courant: 2.0                                                \n"
  "    time_step_change_factor: 1.2                                       \n"
  "                                                                       \n"
  "  solution_options:                                                    \n"
  "    name: unitTestRealmOptions                                         \n"
  "    turbulence_model: laminar                                          \n"
  "    interp_rhou_together_for_mdot: yes                                 \n"
  "    use_consolidated_solver_algorithm: yes                              \n"
  "    reduced_sens_cvfem_poisson: no                                     \n"
  "                                                                       \n"
  "    options:                                                           \n"
  "      - laminar_prandtl:                                               \n"
  "          enthalpy: 0.7                                                \n"
  "      - turbulent_prandtl:                                             \n"
  "          enthalpy: 1.0                                                \n"
  "      - shifted_gradient_operator:                                     \n"
  "          velocity: no                                                 \n"
  "          pressure: no                                                 \n"
  "          mixture_fraction: no                                         \n"
  ;


sierra::nalu::Simulation* create_simulation()
{
    sierra::nalu::Simulation* sim = nullptr;
    YAML::Node sim_node = YAML::Load(naluDefaultInputs);
    sim = new sierra::nalu::Simulation(sim_node);

    return sim;
}


sierra::nalu::Realm& create_mock_realm(stk::mesh::MetaData& meta, stk::mesh::BulkData& bulk, sierra::nalu::Simulation* sim)
{
    sierra::nalu::Realm* realm = nullptr;
    YAML::Node realm_node = YAML::Load(realmDefaultSettings);
    realm = new sierra::nalu::Realm(*sim->realms_, realm_node);

    // Populate solution options
    realm->solutionOptions_->load(realm_node);

    // Set-up mesh metadata and bulkdata
    realm->metaData_ = &meta;
    realm->bulkData_ = &bulk;
    
    return *realm;
}

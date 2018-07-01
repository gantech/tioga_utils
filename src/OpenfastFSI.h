#ifndef OPENFASTFSI_H
#define OPENFASTFSI_H

#include "MotionBase.h"
#include "fsiTurbine.h"
#include "OpenFAST.H"
#include "yaml-cpp/yaml.h"

namespace tioga_nalu {

class OpenfastFSI : public MotionBase
{
public:
    OpenfastFSI(
        stk::mesh::MetaData&,
        stk::mesh::BulkData&,
        const YAML::Node&);

    ~OpenfastFSI();

    virtual void setup();
        
    virtual void initialize(double);

    virtual void execute(double);

private:
    OpenfastFSI() = delete;
    OpenfastFSI(const OpenfastFSI&) = delete;

    void load(const YAML::Node&);

    void deform_mesh(double);

    void compute_mapping();

    void send_loads();

    void get_displacements();

    std::vector<double> origin_{0.0, 0.0, 0.0};

    //Data for coupling to Openfast
    
    fast::OpenFAST FAST;
    
    fast::fastInputs fi ;

    std::vector<fsiTurbine*> fsiTurbineData_;
    
    void read_turbine_data(int iTurb, fast::fastInputs & fi, YAML::Node turbNode);

    void read_inputs(fast::fastInputs & fi, YAML::Node & ofNode);
    
};


} // tioga_nalu

#endif /* OPENFASTFSI_H */

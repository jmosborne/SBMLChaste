#ifndef VANLEEUWEN2007NONDIMODESYSTEMANDSRNMODEL_HPP_
#define VANLEEUWEN2007NONDIMODESYSTEMANDSRNMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <cmath>
#include <iostream>
#include "AbstractOdeSystem.hpp"

class VanLeeuwen2007NonDimOdeSystem : public AbstractOdeSystem
{
private:

    /* Initialise compartments and values. */
    double cytosolmembraneandnucleus;

    /* Initialise model parameters. */
    double K_T;
    double K_C;
    double K_D;
    double p_u;
    double wnt_level;
    double xi_D;
    double xi_Dx;
    double xi_X;
    double xi_C;
    double d_Dx;
    double d_X;
    double d_D;
    double p_c;
    double gamma1;
    double gamma2;
    double s_D;
    double s_X;
    double d_u;
    double s_c;
    double d_c;
    double s_CA;
    double d_CA;
    double s_CT;
    double d_CT;
    double s_A;
    double d_A;
    double s_T;
    double d_T;
    double s_Y;
    double d_Y;
    double ComplexTransitThreshold;


    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

public:

    /* Default constructor. */
    VanLeeuwen2007NonDimOdeSystem(std::vector<double> stateVariables=std::vector<double>());

    /* Destructor. */
    ~VanLeeuwen2007NonDimOdeSystem();


    void Init();

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY);

};

namespace
{
namespace serialization
{
/* Serialize information required to construct a VanLeeuwen2007NonDimOdeSystem. */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const VanLeeuwen2007NonDimOdeSystem * t, const unsigned int file_version)
{
    const std::vector<double> state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}
/* De-serialize constructor parameters and intiialise a VanLeeuwen2007NonDimOdeSystem. */
template<class Archive>
inline void load_construct_data(
    Archive & ar, VanLeeuwen2007NonDimOdeSystem * t, const unsigned int file_version)
{
    std::vector<double> state_variables;
    ar & state_variables;
    
    // Invoke inplace constructor to initialise instance
    ::new(t)VanLeeuwen2007NonDimOdeSystem(state_variables);
}
}
} // namespace ...

/* Define SRN model using Wrappers. */
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<VanLeeuwen2007NonDimOdeSystem,14> VanLeeuwen2007NonDimSrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(VanLeeuwen2007NonDimOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, VanLeeuwen2007NonDimOdeSystem, 14)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(VanLeeuwen2007NonDimSrnModel)

#endif /* VANLEEUWEN2007NONDIMODESYSTEMANDSRNMODEL_HPP_ */
#ifndef TAN2014BETACATENINMODELTWOSBMLTEMPLATEODESYSTEMANDSRNMODEL_HPP_
#define TAN2014BETACATENINMODELTWOSBMLTEMPLATEODESYSTEMANDSRNMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include <cmath>
#include <iostream>
#include "AbstractOdeSystem.hpp"
class Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem : public AbstractOdeSystem
{
 private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeSystem>(*this);
    }

public:

    /* Default constructor. */
    Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem(std::vector<double> stateVariables=std::vector<double>());

    /* Destructor. */
    ~Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem();

    void EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY); 
};

namespace
{
namespace serialization
{
/* Serialize information required to construct a Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem. */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem * t, const unsigned int file_version)
{
    const std::vector<double> state_variables = t->rGetConstStateVariables();
    ar & state_variables;
}
/* De-serialize constructor parameters and intiialise a Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem. */
template<class Archive>
inline void load_construct_data(
    Archive & ar, Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem * t, const unsigned int file_version)
{
    std::vector<double> state_variables;
    ar & state_variables;
    
    // Invoke inplace constructor to initialise instance
    ::new(t)Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem(state_variables);
}
}
} // namespace ...

/* Define Srn model using Wrappers. */
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem,7> Tan2014BetaCateninModelTwoSbmlTemplateSrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem, 7)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Tan2014BetaCateninModelTwoSbmlTemplateSrnModel)

#endif /* TAN2014BETACATENINMODELTWOSBMLTEMPLATEODESYSTEMANDSRNMODEL_HPP_ */
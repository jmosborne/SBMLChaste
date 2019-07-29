#include "Tan2014BetaCateninModelTwoSbmlTemplateSrnModel.hpp"
#include "CellwiseOdeSystemInformation.hpp"
/* SBML ODE System */
Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem::Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem (std::vector<double> stateVariables)
    : AbstractOdeSystem(7)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem>);

    SetDefaultInitialCondition(0, 46.6);
    SetDefaultInitialCondition(1, 581.1);
    SetDefaultInitialCondition(2, 418.9);
    SetDefaultInitialCondition(3, 32.6);
    SetDefaultInitialCondition(4, 516.8);
    SetDefaultInitialCondition(5, 483.2);
    SetDefaultInitialCondition(6, 1.0);

    //// NEW ser intial valiue for parameters, soon overriden by correct values
    this->mParameters.push_back(0.0); // Do in Constructor
    this->mParameters.push_back(0.0); // Do in Constructor

    if (stateVariables != std::vector<double>())    {
        SetStateVariables(stateVariables);
    }
}

Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem::~Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem()
{
}

void Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    /* Define state variables */
    double bcat_cm = rY[0];
    double ligand_cm = rY[1];
    double complex_cm = rY[2];
    double bcat_nu = rY[3];
    double ligand_nu = rY[4];
    double complex_nu = rY[5];
    double bcat_pool = rY[6];

    //// Define Parrameters NOT GENERATED AT THE MOMENT
    double wnt_level = this->mParameters[0]; 
    double degredation_multiplier = this->mParameters[1]; 

    /* Define compartments and their sizes */
/*
    const double compartment = 1.0;
*/
    const double CytosolMembrane = 1.16;
    const double nucleus = 0.65;

    /* Define constant parameters and their values */
 /*   
    const double Bsyn = 1.306;
    const double K_n_active_k = 17.16;
    const double kC_k1 = 1e-05;
    const double kC_k2 = 0.000647;
    const double kN_k1 = 0.0001;
    const double kN_k2 = 0.00349;
    const double kdiffusion_k = 39.13;
    const double K_c_active_k = 4.5;
*/
    // RESCALED PARAMETERS TO MAKE THINGS HAPPEN FASTER 
    const double kdegradation = 0.0163*20.0;
    const double Bsyn = 1.306*20.0;
    const double K_n_active_k = 17.16*20.0;
    const double kC_k1 = (1e-05)*20.0;
    const double kC_k2 = 0.000647*20.0;
    const double kN_k1 = 0.0001*20.0;
    const double kN_k2 = 0.00349*20.0;
    const double kdiffusion_k = 39.13*20.0;
    const double K_c_active_k = 4.5*20.0;

    /* Define the rules required for reactions. */

    /* Define the reactions in this model. */
    double Bsynthesis = Bsyn * bcat_pool;

    //double kDegradation = kdegradation * bcat_cm;
    double kDegradation = degredation_multiplier * (1.0 - 0.5 * wnt_level) * kdegradation * bcat_cm; // Modified by Ozzy to look like equations in pdf

    double kC = kC_k1 * bcat_cm * ligand_cm - kC_k2 * complex_cm;
 
    double kN = kN_k1 * bcat_nu * ligand_nu - kN_k2 * complex_nu;

    double kdiffusion = kdiffusion_k * (bcat_cm - bcat_nu);

    double K_c_active = K_c_active_k * bcat_cm;

    double K_n_active = K_n_active_k * bcat_nu;

    /* Generated RHS
    rDY[0] = (Bsynthesis - kDegradation - kC - kdiffusion - K_c_active + K_n_active) / CytosolMembrane; // dbcat_cm/dt
    rDY[1] = (kC) / CytosolMembrane; // dligand_cm/dt
    rDY[2] = (kC) / CytosolMembrane; // dcomplex_cm/dt
    rDY[3] = (kN + kdiffusion + K_c_active - K_n_active) / nucleus; // dbcat_nu/dt
    rDY[4] = (kN) / nucleus; // dligand_nu/dt
    rDY[5] = (kN) / nucleus; // dcomplex_nu/dt
    rDY[6] = (Bsynthesis + kDegradation) / compartment; // dbcat_pool/dt
    */

    // Modified by Ozzy to look like equations in pdf
    rDY[0] = Bsynthesis - kDegradation - kC + (- kdiffusion - K_c_active + K_n_active) / CytosolMembrane; // dbcat_cm/dt
    rDY[1] = -(kC);// / CytosolMembrane; // dligand_cm/dt
    rDY[2] = (kC);// / CytosolMembrane; // dcomplex_cm/dt
    rDY[3] = -kN + (kdiffusion + K_c_active - K_n_active) / nucleus; // dbcat_nu/dt
    rDY[4] = -(kN);// / nucleus; // dligand_nu/dt
    rDY[5] = (kN);// / nucleus; // dcomplex_nu/dt
    rDY[6] = 0.0;//(Bsynthesis + kDegradation) / compartment; // dbcat_pool/dt

}

template<>
void CellwiseOdeSystemInformation<Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("bcat_cm");
    this->mVariableUnits.push_back("non-dim");
    //this->mInitialConditions.push_back(46.6);
    this->mInitialConditions.push_back(160.2457);

    this->mVariableNames.push_back("ligand_cm");
    this->mVariableUnits.push_back("non-dim");
    //this->mInitialConditions.push_back(581.1);
    this->mInitialConditions.push_back(287.6244);

    this->mVariableNames.push_back("complex_cm");
    this->mVariableUnits.push_back("non-dim");
    //this->mInitialConditions.push_back(418.9);
    this->mInitialConditions.push_back(712.3756);

    this->mVariableNames.push_back("bcat_nu");
    this->mVariableUnits.push_back("non-dim");
    //this->mInitialConditions.push_back(32.6);
    this->mInitialConditions.push_back(124.2053);

    this->mVariableNames.push_back("ligand_nu");
    this->mVariableUnits.push_back("non-dim");
    //this->mInitialConditions.push_back(516.8);
    this->mInitialConditions.push_back(219.3515);

    this->mVariableNames.push_back("complex_nu");
    this->mVariableUnits.push_back("non-dim");
    //this->mInitialConditions.push_back(483.2);
    this->mInitialConditions.push_back(780.6485);

    this->mVariableNames.push_back("bcat_pool");
    this->mVariableUnits.push_back("non-dim");
    //this->mInitialConditions.push_back(0.0);
    this->mInitialConditions.push_back(1.0);

    //Adding parameters NOT GENERATED YET
    this->mParameterNames.push_back("wnt_level");
    this->mParameterUnits.push_back("non-dim");
    ////this->mParameters.push_back(0.0); // Do in Constructor

    this->mParameterNames.push_back("degredation_multiplier");
    this->mParameterUnits.push_back("non-dim");
    ////this->mParameters.push_back(0.0); // Do in Constructor

    this->mInitialised = true;
}

/* Define Srn model using Wrappers. */
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem,7> Tan2014BetaCateninModelTwoSbmlTemplateSrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, Tan2014BetaCateninModelTwoSbmlTemplateOdeSystem, 7)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Tan2014BetaCateninModelTwoSbmlTemplateSrnModel)


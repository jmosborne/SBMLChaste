#include "VanLeeuwen2007SrnModel.hpp"
#include "CellwiseOdeSystemInformation.hpp"
/* SBML ODE System */
VanLeeuwen2007OdeSystem::VanLeeuwen2007OdeSystem (std::vector<double> stateVariables)
    : AbstractOdeSystem(14)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<VanLeeuwen2007OdeSystem>);

    Init();

    SetDefaultInitialCondition(0, 0.067);
    SetDefaultInitialCondition(1, 0.67);
    SetDefaultInitialCondition(2, 2.54);
    SetDefaultInitialCondition(3, 0.45);
    SetDefaultInitialCondition(4, 0.0);
    SetDefaultInitialCondition(5, 10.0);
    SetDefaultInitialCondition(6, 18.14);
    SetDefaultInitialCondition(7, 25.0);
    SetDefaultInitialCondition(8, 2.54);
    SetDefaultInitialCondition(9, 0.0);
    SetDefaultInitialCondition(10, 0.48);
    SetDefaultInitialCondition(11, 2.54);
    SetDefaultInitialCondition(12, 2.54);
    SetDefaultInitialCondition(13, 1.0);

    this->mParameters.push_back(0.0);//wnt_level
    this->mParameters.push_back(1.0);//gamma1
    this->mParameters.push_back(1.0);//gamma2
    this->mParameters.push_back(1.0);//ComplexTransitThreshold

    if (stateVariables != std::vector<double>())    {
        SetStateVariables(stateVariables);
    }
}

VanLeeuwen2007OdeSystem::~VanLeeuwen2007OdeSystem()
{
}

void VanLeeuwen2007OdeSystem::Init()
 {
    /* Initialise the parameters. */
    cytosolmembraneandnucleus = 1.0;
    K_T = 50.0;
    K_C = 200.0;
    K_D = 5.0;
    p_u = 100.0;
    wnt_level = 0.0;
    gamma1 = 1.0;
    gamma2 = 1.0;
    xi_D = 5.0;
    xi_Dx = 5.0;
    xi_X = 200.0;
    xi_C = 0.0;
    s_D = 100.0;
    d_Dx = 5.0;
    s_X = 10.0;
    d_X = 100.0;
    d_u = 50.0;
    s_c = 25.0;
    d_c = 1.0;
    s_CA = 250.0;
    d_CA = 350.0;
    s_CT = 30.0;
    d_CT = 750.0;
    p_c = 0.0;
    s_A = 20.0;
    d_A = 2.0;
    s_T = 10.0;
    d_T = 0.4;
    s_Y = 10.0;
    d_Y = 1.0;
    d_D = 5.0;
    ComplexTransitThreshold = 1.0;
}

void VanLeeuwen2007OdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    /* Define state variables */
    double X = rY[0]; // X
    double D = rY[1]; // D
    double C_o = rY[2]; // C_o
    double C_u = rY[3]; // C_u
    double C_c = rY[4]; // C_c
    double A = rY[5]; // A
    double C_A = rY[6]; // C_A
    double T = rY[7]; // T
    double C_oT = rY[8]; // C_oT
    double C_cT = rY[9]; // C_cT
    double Y = rY[10]; // Y
    //double C_F = rY[11]; // C_F
    //double C_T = rY[12]; // C_T
    double drag = rY[13]; // drag

    /* Define state parameters */
    double wnt_level = this->mParameters[0]; // wnt_level
    double gamma1 = this->mParameters[1]; // gamma1
    double gamma2 = this->mParameters[2]; // gamma2
    //double ComplexTransitThreshold = this->mParameters[3]; // ComplexTransitThreshold

     /* Define algebraic rules. */
    //C_F = C_o + C_c;
    //C_T = C_oT + C_cT;
    ////////////////// MODIFIED ///////////////
    //drag = fmax((C_A - 100) / 3, 1);
    drag = fmax((C_A - 98.84798) / 3.98617, 1);
    ///////////////////////////////////////////

    /* Define the reactions in this model. */
    // r1
    double r1 = s_D * gamma1 * X;

    // r2
    double r2 = (d_Dx + wnt_level * xi_Dx) * D;

    // r22
    double r22 = s_X * cytosolmembraneandnucleus;

    // r23
    double r23 = (d_X + wnt_level * xi_X) * X;

    // r7
    double r7 = p_u * gamma2 * C_o * D / (C_o + C_c + K_D);

    // r16
    double r16 = p_u * gamma2 * C_c * D / (C_c + C_o + K_D);

    // r8
    double r8 = d_u * C_u * cytosolmembraneandnucleus;

    // r3
    double r3 = s_c * cytosolmembraneandnucleus;

    // r4
    double r4 = d_c * C_o * cytosolmembraneandnucleus;

    // r17
    double r17 = d_c * C_c * cytosolmembraneandnucleus;

    // r9
    double r9 = s_CA * C_o * A * cytosolmembraneandnucleus;

    // r10
    double r10 = d_CA * C_A * cytosolmembraneandnucleus;

    // r11
    double r11 = s_CT * C_o * T * cytosolmembraneandnucleus;

    // r18
    double r18 = s_CT * C_c * T * cytosolmembraneandnucleus;

    // r12
    double r12 = d_CT * C_oT * cytosolmembraneandnucleus;

    // r19
    double r19 = d_CT * C_cT * cytosolmembraneandnucleus;

    // r15
    double r15 = (p_c + wnt_level * xi_C) * C_o / (C_o + K_C);

    // r5
    double r5 = s_A * cytosolmembraneandnucleus;

    // r6
    double r6 = d_A * A * cytosolmembraneandnucleus;

    // r20
    double r20 = s_T * cytosolmembraneandnucleus;

    // r21
    double r21 = d_T * T * cytosolmembraneandnucleus;

    // r13
    double r13 = s_Y * (C_oT + C_cT) / (C_oT + C_cT + K_T);

    // r14
    double r14 = d_Y * Y * cytosolmembraneandnucleus;

    // r24
    double r24 = (d_D + wnt_level * xi_D) * D;


    rDY[0] = ( - r1 + r2 + r22 - r23) / cytosolmembraneandnucleus; // dX/dt
    rDY[1] = (r1 - r2 + r7 - r7 + r16 - r16 - r24) / cytosolmembraneandnucleus; // dD/dt
    rDY[2] = ( - r7 + r3 - r4 - r9 + r10 - r11 + r12 - r15) / cytosolmembraneandnucleus; // dC_o/dt
    rDY[3] = (r7 + r16 - r8) / cytosolmembraneandnucleus; // dC_u/dt
    rDY[4] = ( - r16 - r17 - r18 + r19 + r15) / cytosolmembraneandnucleus; // dC_c/dt
    rDY[5] = ( - r9 + r10 + r5 - r6) / cytosolmembraneandnucleus; // dA/dt
    rDY[6] = (r9 - r10) / cytosolmembraneandnucleus; // dC_A/dt
    rDY[7] = ( - r11 - r18 + r12 + r19 + r20 - r21) / cytosolmembraneandnucleus; // dT/dt
    rDY[8] = (r11 - r12 + r13 - r13) / cytosolmembraneandnucleus; // dC_oT/dt
    rDY[9] = (r18 - r19) / cytosolmembraneandnucleus; // dC_cT/dt
    rDY[10] = (r13 - r14) / cytosolmembraneandnucleus; // dY/dt
    rDY[11] = (rDY[2]+ rDY[4]) / cytosolmembraneandnucleus; // dC_F/dt
    rDY[12] = (rDY[8]+ rDY[9]) / cytosolmembraneandnucleus; // dC_T/dt
    rDY[13] = (drag - rY[13]) / cytosolmembraneandnucleus; // ddrag/dt

}

template<>
void CellwiseOdeSystemInformation<VanLeeuwen2007OdeSystem>::Initialise()
{
    this->mVariableNames.push_back("X");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.067);

    this->mVariableNames.push_back("D");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.67);

    this->mVariableNames.push_back("C_o");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(2.54);

    this->mVariableNames.push_back("C_u");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.45);

    this->mVariableNames.push_back("C_c");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("A");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(10.0);

    this->mVariableNames.push_back("C_A");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(18.14);

    this->mVariableNames.push_back("T");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(25.0);

    this->mVariableNames.push_back("C_oT");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(2.54);

    this->mVariableNames.push_back("C_cT");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("Y");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.48);

    this->mVariableNames.push_back("C_F");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(2.54);

    this->mVariableNames.push_back("C_T");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(2.54);

    this->mVariableNames.push_back("drag");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(1.0);

    this->mParameterNames.push_back("wnt_level");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("gamma1");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("gamma2");
    this->mParameterUnits.push_back("non-dim");

    this->mParameterNames.push_back("ComplexTransitThreshold");
    this->mParameterUnits.push_back("non-dim");


    this->mInitialised = true;
}

/* Define SRN model using Wrappers. */
#include "SbmlSrnWrapperModel.hpp"
#include "SbmlSrnWrapperModel.cpp"

typedef SbmlSrnWrapperModel<VanLeeuwen2007OdeSystem,14> VanLeeuwen2007SrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(VanLeeuwen2007OdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, VanLeeuwen2007OdeSystem, 14)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(VanLeeuwen2007SrnModel)


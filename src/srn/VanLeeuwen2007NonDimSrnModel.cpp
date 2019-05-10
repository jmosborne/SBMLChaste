#include "VanLeeuwen2007NonDimSrnModel.hpp"
#include "CellwiseOdeSystemInformation.hpp"
/* SBML ODE System */
VanLeeuwen2007NonDimOdeSystem::VanLeeuwen2007NonDimOdeSystem (std::vector<double> stateVariables)
    : AbstractOdeSystem(14)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<VanLeeuwen2007NonDimOdeSystem>);

    Init();

    SetDefaultInitialCondition(0, 0.268);
    SetDefaultInitialCondition(1, 2.68);
    SetDefaultInitialCondition(2, 76.2);
    SetDefaultInitialCondition(3, 13.5);
    SetDefaultInitialCondition(4, 0.0);
    SetDefaultInitialCondition(5, 300.0);
    SetDefaultInitialCondition(6, 544.2);
    SetDefaultInitialCondition(7, 750.0);
    SetDefaultInitialCondition(8, 76.2);
    SetDefaultInitialCondition(9, 0.0);
    SetDefaultInitialCondition(10, 36.0);
    SetDefaultInitialCondition(11, 76.2);
    SetDefaultInitialCondition(12, 76.2);
    SetDefaultInitialCondition(13, 1.0);

    this->mParameters.push_back(0.0);//wnt_level
    this->mParameters.push_back(1.0);//gamma1
    this->mParameters.push_back(1.0);//gamma2
    this->mParameters.push_back(1.0);//ComplexTransitThreshold

    if (stateVariables != std::vector<double>())    {
        SetStateVariables(stateVariables);
    }
}

VanLeeuwen2007NonDimOdeSystem::~VanLeeuwen2007NonDimOdeSystem()
{
}

void VanLeeuwen2007NonDimOdeSystem::Init()
 {
    /* Initialise the parameters. */
    cytosolmembraneandnucleus = 1.0;
    K_T = 1500.0;
    K_C = 6000.0;
    K_D = 150.0;
    p_u = 1.0;
    wnt_level = 0.0;
    xi_D = 0.00666666666667;
    xi_Dx = 0.00666666666667;
    xi_X = 0.266666666667;
    xi_C = 0.0;
    d_Dx = 0.00666666666667;
    d_X = 0.133333333333;
    d_D = 0.00666666666667;
    p_c = 0.0;
    gamma1 = 1.0;
    gamma2 = 1.0;
    s_D = 0.133333333333;
    s_X = 0.0533333333333;
    d_u = 0.0666666666667;
    s_c = 1.0;
    d_c = 0.00133333333333;
    s_CA = 0.0111111111111;
    d_CA = 0.466666666667;
    s_CT = 0.00133333333333;
    d_CT = 1.0;
    s_A = 0.8;
    d_A = 0.00266666666667;
    s_T = 0.4;
    d_T = 0.000533333333333;
    s_Y = 1.0;
    d_Y = 0.00133333333333;
    ComplexTransitThreshold = 1.0;
}

void VanLeeuwen2007NonDimOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
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
    double C_F = rY[11]; // C_F
    double C_T = rY[12]; // C_T
    double drag = rY[13]; // drag

    /* Define state parameters */
    double wnt_level = this->mParameters[0]; // wnt_level
    double gamma1 = this->mParameters[1]; // gamma1
    double gamma2 = this->mParameters[2]; // gamma2
    double ComplexTransitThreshold = this->mParameters[3]; // ComplexTransitThreshold

     /* Define algebraic rules. */
    C_F = C_o + C_c;
    C_T = C_oT + C_cT;
    drag = fmax((C_A - 2300) / 36, 1);
    //drag = 1.0;

    /* Define the reactions in this model. */
    // r1
    double r1 = s_D * gamma1 * X;

    // r2
    double r2 = (d_Dx + wnt_level * xi_Dx) * D;

    // r22
    double r22 = s_X;

    // r23
    double r23 = (d_X + wnt_level * xi_X) * X;

    // r7
    double r7 = p_u * gamma2 * C_o * D / (C_o + C_c + K_D);

    // r16
    double r16 = p_u * gamma2 * C_c * D / (C_c + C_o + K_D);

    // r8
    double r8 = d_u * C_u;

    // r3
    double r3 = s_c;

    // r4
    double r4 = d_c * C_o;

    // r17
    double r17 = d_c * C_c;

    // r9
    double r9 = s_CA * C_o * A;

    // r10
    double r10 = d_CA * C_A;

    // r11
    double r11 = s_CT * C_o * T;

    // r18
    double r18 = s_CT * C_c * T;

    // r12
    double r12 = d_CT * C_oT;

    // r19
    double r19 = d_CT * C_cT;

    // r15
    double r15 = (p_c + wnt_level * xi_C) * C_o / (C_o + K_C);

    // r5
    double r5 = s_A;

    // r6
    double r6 = d_A * A;

    // r20
    double r20 = s_T;

    // r21
    double r21 = d_T * T;

    // r13
    double r13 = s_Y * (C_oT + C_cT) / (C_oT + C_cT + K_T);

    // r14
    double r14 = d_Y * Y;

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
void CellwiseOdeSystemInformation<VanLeeuwen2007NonDimOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("X");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.268);

    this->mVariableNames.push_back("D");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(2.68);

    this->mVariableNames.push_back("C_o");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(76.2);

    this->mVariableNames.push_back("C_u");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(13.5);

    this->mVariableNames.push_back("C_c");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("A");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(300.0);

    this->mVariableNames.push_back("C_A");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(544.2);

    this->mVariableNames.push_back("T");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(750.0);

    this->mVariableNames.push_back("C_oT");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(76.2);

    this->mVariableNames.push_back("C_cT");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(0.0);

    this->mVariableNames.push_back("Y");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(36.0);

    this->mVariableNames.push_back("C_F");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(76.2);

    this->mVariableNames.push_back("C_T");
    this->mVariableUnits.push_back("non-dim");
    this->mInitialConditions.push_back(76.2);

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

typedef SbmlSrnWrapperModel<VanLeeuwen2007NonDimOdeSystem,14> VanLeeuwen2007NonDimSrnModel;

// Declare identifiers for the serializer
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(VanLeeuwen2007NonDimOdeSystem)
EXPORT_TEMPLATE_CLASS2(SbmlSrnWrapperModel, VanLeeuwen2007NonDimOdeSystem, 14)

#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(VanLeeuwen2007NonDimSrnModel)


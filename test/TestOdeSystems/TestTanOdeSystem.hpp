/*

Copyright (c) 2005-2016, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTTANODEDESYSTEM_HPP_
#define TESTTANODEDESYSTEM_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <vector>
#include <iostream>

//#include "Tan2014BetaCateninModelTwoSbmlTemplateSrnModel.hpp"
#include "Tan2014SbmlSrnModel.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "ColumnDataWriter.hpp"
#include "Timer.hpp"
#include "Debug.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestTanOdeSystem : public CxxTest::TestSuite
{
public:

    void TestTanEquation()
    {
        Tan2014SbmlOdeSystem tan_system;

        double time = 0.0;
        std::vector<double> initial_conditions;
        initial_conditions.push_back(46.6);
        initial_conditions.push_back(581.1);
        initial_conditions.push_back(418.9);
        initial_conditions.push_back(32.6);
        initial_conditions.push_back(516.8);
        initial_conditions.push_back(483.2);
        initial_conditions.push_back(1.0);

        std::vector<double> derivs(initial_conditions.size());
        tan_system.EvaluateYDerivatives(time, initial_conditions, derivs);

        // Test derivatives are correct (default wnt-level=0)
        TS_ASSERT_DELTA(derivs[0], -170.2326546448276*60.0, 1e-4);
        TS_ASSERT_DELTA(derivs[1], 0.0002357000000*60.0, 1e-5);
        TS_ASSERT_DELTA(derivs[2], -0.0002357000000*60.0, 1e-5);
        TS_ASSERT_DELTA(derivs[3], 304.7769846153847*60.0, 1e-5);
        TS_ASSERT_DELTA(derivs[4], 0.0016000000000*60.0, 1e-5);
        TS_ASSERT_DELTA(derivs[5], -0.0016000000000*60.0, 1e-5);
        TS_ASSERT_DELTA(derivs[6], 0.0*60.0, 1e-5);

        //Change Wnt level and check still OK
        tan_system.SetParameter("wnt_level",1);
        tan_system.EvaluateYDerivatives(time, initial_conditions, derivs);

        TS_ASSERT_DELTA(derivs[0], -169.85286464482763335581694263965*60.0, 1e-3);
        TS_ASSERT_DELTA(derivs[1], 0.0002357000000*60.0, 1e-5);
        TS_ASSERT_DELTA(derivs[2], -0.0002357000000*60.0, 1e-5);
        TS_ASSERT_DELTA(derivs[3], 304.7769846153847*60.0, 1e-5);
        TS_ASSERT_DELTA(derivs[4], 0.0016000000000*60.0, 1e-5);
        TS_ASSERT_DELTA(derivs[5], -0.0016000000000*60.0, 1e-5);
        TS_ASSERT_DELTA(derivs[6], 0.0*60.0, 1e-5);
    }

    void TestTanWithChasteSolver() throw(Exception)
    {
        Tan2014SbmlOdeSystem tan_system;

        // Solve system using RK4 solver

        double dt = 0.0001;

        //RK4 solver solution worked out
        RungeKutta4IvpOdeSolver rk4_solver;

        std::vector<double> state_variables = tan_system.GetInitialConditions();

        Timer::Reset();
        OdeSolution solutions = rk4_solver.Solve(&tan_system, state_variables, 0.0, 1000.0, dt, dt);
        Timer::Print("1. Tan RK4");
        
        unsigned end = solutions.rGetSolutions().size() - 1;

        // The following code provides nice output for gnuplot
        // use the command
        // plot "tan.dat" u 1:2
        // or
        // plot "tan.dat" u 1:3 etc. for the various proteins...

        // OutputFileHandler handler("");
        // out_stream file=handler.OpenOutputFile("tan.dat");
        // for (unsigned i=0; i<=end; i++)
        // {
        //     (*file) << solutions.rGetTimes()[i]<< "\t" << solutions.rGetSolutions()[i][0] << "\t" << solutions.rGetSolutions()[i][1] << "\t" << solutions.rGetSolutions()[i][2] << "\t" << solutions.rGetSolutions()[i][3] << "\t" << solutions.rGetSolutions()[i][4] << "\t" << solutions.rGetSolutions()[i][5] << "\n" << std::flush;
        // }
        // file->close();

        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0],80.123027932584804489124508108944, 2e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1],446.75578463039073540130630135536, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2],553.24421536960869616450509056449, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3],62.102819644032457802040880778804, 2e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4],359.76798287328375636207056231797, 5e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5],640.23201712671584573399741202593, 5e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][6],1.0, 1e-5);
    }

 void TestTanWithCvodeSolver() throw(Exception)
    {
        Tan2014SbmlOdeSystem ode_system;

        double end_time = 1000;
        double h_value = 0.01;

        CvodeAdaptor solver;
        OdeSolution solutions;

        std::vector<double> state_variables = ode_system.GetInitialConditions();
        
        Timer::Reset();
        solutions = solver.Solve(&ode_system, state_variables, 0.0, end_time, h_value, 0.1);
        int last = solutions.GetNumberOfTimeSteps();
        Timer::Print("1. Tan CVODE");
        
        unsigned end = solutions.rGetSolutions().size() - 1;

        // The following code provides nice output for gnuplot
        // use the command
        // plot "tan.dat" u 1:2
        // or
        // plot "tan.dat" u 1:3 etc. for the various proteins...

        // OutputFileHandler handler("");
        // out_stream file=handler.OpenOutputFile("tan.dat");
        // for (unsigned i=0; i<=end; i++)
        // {
        //     (*file) << solutions.rGetTimes()[i]<< "\t" << solutions.rGetSolutions()[i][0] << "\t" << solutions.rGetSolutions()[i][1] << "\t" << solutions.rGetSolutions()[i][2] << "\t" << solutions.rGetSolutions()[i][3] << "\t" << solutions.rGetSolutions()[i][4] << "\t" << solutions.rGetSolutions()[i][5] << "\n" << std::flush;
        // }
        // file->close();

        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0],80.123027932584804489124508108944, 2e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1],446.75578463039073540130630135536, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][2],553.24421536960869616450509056449, 1e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][3],62.102819644032457802040880778804, 2e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][4],359.76798287328375636207056231797, 5e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][5],640.23201712671584573399741202593, 5e-2);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][6],1.0, 1e-5);
    }

    void noTestArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "tn_ode.arch";

        {
            Tan2014SbmlOdeSystem ode_system;

            ode_system.SetDefaultInitialCondition(2, 3.25);

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 6u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.0999, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[1], 0.9890, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[2], 3.2500, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[3], 1.4211, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[4], 0.6728, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[5], 0.4854, 1e-4);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive ODE system
            AbstractOdeSystem* const p_const_ode_system = &ode_system;
            output_arch << p_const_ode_system;
        }

        {
            AbstractOdeSystem* p_ode_system;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_ode_system;

            // Check that archiving worked correctly
            std::vector<double> initial_conditions = p_ode_system->GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 6u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.0999, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[1], 0.9890, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[2], 3.2500, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[3], 1.4211, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[4], 0.6728, 1e-4);
            TS_ASSERT_DELTA(initial_conditions[5], 0.4854, 1e-4);

            // Tidy up
            delete p_ode_system;
        }
    }
};

#endif /*TESTTANODEDESYSTEM_HPP_*/

/*

Copyright (c) 2005-2015, University of Oxford.
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

#include "SimpleCellDataComplexNuclearCellCycleModel.hpp"
#include "Exception.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "QuiescentCellLabel.hpp"
#include "Debug.hpp"

SimpleCellDataComplexNuclearCellCycleModel::SimpleCellDataComplexNuclearCellCycleModel()
    : mMinimumCellSize(0.0) //sqrt(3.0)/8.0) // Default value
{
}

SimpleCellDataComplexNuclearCellCycleModel::SimpleCellDataComplexNuclearCellCycleModel(const SimpleCellDataComplexNuclearCellCycleModel& rModel)
   : AbstractSimplePhaseBasedCellCycleModel(rModel),
     mMinimumCellSize(rModel.mMinimumCellSize)
{
    /*
     * Set each member variable of the new cell-cycle model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new cell-cycle model's member variables will already
     * have been correctly initialized in its constructor or parent classes.
     *
     * Note 2: one or more of the new cell-cycle model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new cell-cycle model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     *
     */
}


AbstractCellCycleModel* SimpleCellDataComplexNuclearCellCycleModel::CreateCellCycleModel()
{
    return new SimpleCellDataComplexNuclearCellCycleModel(*this);
}

void SimpleCellDataComplexNuclearCellCycleModel::SetG1Duration()
{
    assert(mpCell != NULL);

    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
    {
        mG1Duration = p_gen->NormalRandomDeviate(GetTransitCellG1Duration(), 1.0);
    }
    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }
    else
    {
        NEVER_REACHED;
    }

    // Check that the normal random deviate has not returned a small or negative G1 duration
    if (mG1Duration < mMinimumGapDuration)
    {
        mG1Duration = mMinimumGapDuration;
    }
}

double SimpleCellDataComplexNuclearCellCycleModel::GetComplexLevel()
{
    assert(mpCell != NULL);

    for (unsigned i=1; i < mpCell->GetCellData()->GetNumItems(); i++)
    {
        if( mpCell->GetCellData()->GetKeys()[i] == "complex_nu")
        {
            return mpCell->GetCellData()->GetItem("complex_nu");
        }
        else if (mpCell->GetCellData()->GetKeys()[i] == "C_T")
        {
            return mpCell->GetCellData()->GetItem("C_T");
        }
    }
    return DOUBLE_UNSET;
}

void SimpleCellDataComplexNuclearCellCycleModel::UpdateCellCyclePhase()
{
    // Implement basic contact inhibition

    // // Get cell volume
    // double cell_volume = mpCell->GetCellData()->GetItem("volume");

    // // Removes the cell label
    // mpCell->RemoveCellProperty<QuiescentCellLabel>();

    // if (mCurrentCellCyclePhase == G_ONE_PHASE)
    // {
    //     // Update G1 duration based on cell volume
    //     if (cell_volume < mMinimumCellSize)
    //     {
    //         // Increase the G1 duration so holds cell in G1 if too small
    //         mG1Duration += SimulationTime::Instance()->GetTimeStep();

    //         /*
    //          * This method is usually called within a CellBasedSimulation, after the CellPopulation
    //          * has called CellPropertyRegistry::TakeOwnership(). This means that were we to call
    //          * CellPropertyRegistry::Instance() here when adding the QuiescentCellLabel, we would be creating
    //          * a new CellPropertyRegistry. In this case the QuiescentCellLabel's cell count would be incorrect.
    //          * We must therefore access the QuiescentCellLabel via the cell's CellPropertyCollection.
    //          */
    //         boost::shared_ptr<AbstractCellProperty> p_label =
    //             mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<QuiescentCellLabel>();
    //         mpCell->AddCellProperty(p_label);
    //     }
    // }
    // End contact inhibition

    // The cell can divide if the Complex concentration >= Complex_division_threshold
    double Complex_level = GetComplexLevel();
    // Set the cell type to TransitCellProliferativeType if the Complex stimulus exceeds Complex_division_threshold
    if (Complex_level >= mpCell->GetCellData()->GetItem("ComplexTransitThreshold"))
    {
        boost::shared_ptr<AbstractCellProperty> p_transit_type = 
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_transit_type);

        // Reset the Cell Cycle phase incase been diferentiated for longer than the cell cycle
        mCurrentCellCyclePhase = G_TWO_PHASE;

        if (mG1Duration == DBL_MAX) // If coming out of Differentiated type then reset the G1 duration
        {
            SetG1Duration();
            PRINT_VARIABLE(mG1Duration);
        }
    }
    else
    {
        // The cell is set to have DifferentiatedCellProliferativeType and so in G0 phase
        boost::shared_ptr<AbstractCellProperty> p_diff_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_diff_type);
    }
    AbstractSimplePhaseBasedCellCycleModel::UpdateCellCyclePhase();
}

void SimpleCellDataComplexNuclearCellCycleModel::InitialiseDaughterCell()
{
    AbstractSimplePhaseBasedCellCycleModel::InitialiseDaughterCell();
}

bool SimpleCellDataComplexNuclearCellCycleModel::CanCellTerminallyDifferentiate()
{
    return false;
}

double SimpleCellDataComplexNuclearCellCycleModel::GetMinimumCellSize()
{
    return mMinimumCellSize;
}

void SimpleCellDataComplexNuclearCellCycleModel::SetMinimumCellSize(double minimumCellSize)
{
    mMinimumCellSize = minimumCellSize;
}

void SimpleCellDataComplexNuclearCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MinimumCellSize>" << mMinimumCellSize << "</MinimumCellSize>\n";
 
    // Call method on direct parent class
    AbstractSimplePhaseBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SimpleCellDataComplexNuclearCellCycleModel)

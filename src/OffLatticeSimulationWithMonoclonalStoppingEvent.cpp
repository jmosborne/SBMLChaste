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

#include "OffLatticeSimulationWithMonoclonalStoppingEvent.hpp"
#include "CellLabel.hpp"

bool OffLatticeSimulationWithMonoclonalStoppingEvent::StoppingEventHasOccurred()
{
   boost::shared_ptr<AbstractCellProperty> p_labeled(CellPropertyRegistry::Instance()->Get<CellLabel>());

	 unsigned num_labeled = p_labeled->GetCellCount();
	 unsigned num_non_labeled = this->mrCellPopulation.GetNumRealCells() - num_labeled;

   return  ( (num_labeled ==0 || num_non_labeled == 0 ||  this->mrCellPopulation.GetNumRealCells() > mMaxNumCells) 
             && (SimulationTime::Instance()->GetTime() > mStartTime) );
}

OffLatticeSimulationWithMonoclonalStoppingEvent::OffLatticeSimulationWithMonoclonalStoppingEvent(
            AbstractCellPopulation<2>& rCellPopulation,
            bool deleteCellPopulationInDestructor,
            bool initialiseCells)
    : OffLatticeSimulation<2>(rCellPopulation,
    deleteCellPopulationInDestructor,
    initialiseCells)
{
    mMaxNumCells = 1000;
    mStartTime = 100;
}

void OffLatticeSimulationWithMonoclonalStoppingEvent::SetMaxNumCells(unsigned maxNumCells)
{
    mMaxNumCells = maxNumCells;
}

void OffLatticeSimulationWithMonoclonalStoppingEvent::SetStartTime(double startTime)
{
    mStartTime = startTime;
}



// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(OffLatticeSimulationWithMonoclonalStoppingEvent)

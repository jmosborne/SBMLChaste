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

#ifndef OFFLATTICESIMULATIONWITHMONOCLONALSTOPPINGEVENT_HPP_
#define OFFLATTICESIMULATIONWITHMONOCLONALSTOPPINGEVENT_HPP_

#include "OffLatticeSimulation.hpp"

/**
 * Simple subclass of OffLatticeSimulation which just overloads StoppingEventHasOccurred
 * for stopping a simulation when its become polyclonal of a given degree..
 */
class OffLatticeSimulationWithMonoclonalStoppingEvent : public OffLatticeSimulation<2>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    
    /**
     * Save or restore the simulation.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<OffLatticeSimulation<2> >(*this);
        archive & mMaxNumCells;
        archive & mStartTime;
    }

    unsigned mMaxNumCells;

    double mStartTime;
    
	/** Define a stopping event which says stop when there are no mutant or no healthy cells. */
    bool StoppingEventHasOccurred();

public:
    OffLatticeSimulationWithMonoclonalStoppingEvent(AbstractCellPopulation<2>& rCellPopulation,
                                                    bool deleteCellPopulationInDestructor=false,
                                                    bool initialiseCells=true);

    void SetMaxNumCells(unsigned maxNumCells);      

    void SetStartTime(double startTime);
};

// Serialization for Boost >= 1.36
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(OffLatticeSimulationWithMonoclonalStoppingEvent)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a OffLatticeSimulationWithMyStoppingEvent.
 */
template<class Archive>
inline void save_construct_data(
    Archive & ar, const OffLatticeSimulationWithMonoclonalStoppingEvent * t, unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<2>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a OffLatticeSimulationWithMyStoppingEvent.
 */
template<class Archive>
inline void load_construct_data(
	Archive & ar, OffLatticeSimulationWithMonoclonalStoppingEvent * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<2>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)OffLatticeSimulationWithMonoclonalStoppingEvent(*p_cell_population, true, false);
}
}
} // namespace

#endif /*OFFLATTICESIMULATIONWITHMONOCLONALSTOPPINGEVENT_HPP_*/

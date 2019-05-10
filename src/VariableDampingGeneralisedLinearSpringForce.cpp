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

#include "VariableDampingGeneralisedLinearSpringForce.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "Debug.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VariableDampingGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::VariableDampingGeneralisedLinearSpringForce()
   : GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VariableDampingGeneralisedLinearSpringForce<ELEMENT_DIM,SPACE_DIM>::AddForceContribution(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation)
{
    // Throw an exception message if not using a subclass of AbstractCentreBasedCellPopulation
    if (dynamic_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("Subclasses of AbstractTwoBodyInteractionForce are to be used with subclasses of AbstractCentreBasedCellPopulation only");
    }

    ///\todo this could be tidied by using the rGetNodePairs for all populations and moving the below calculation into the MutableMesh.
    if (bool(dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation)))
    {
        MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_static_cast_cell_population = static_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation);

        // Iterate over all springs and add force contributions
        for (typename MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>::SpringIterator spring_iterator = p_static_cast_cell_population->SpringsBegin();
             spring_iterator != p_static_cast_cell_population->SpringsEnd();
             ++spring_iterator)
        {
            unsigned nodeA_global_index = spring_iterator.GetNodeA()->GetIndex();
            unsigned nodeB_global_index = spring_iterator.GetNodeB()->GetIndex();

            // Calculate the force between nodes
            c_vector<double, SPACE_DIM> force = this->CalculateForceBetweenNodes(nodeA_global_index, nodeB_global_index, rCellPopulation);

double drag_cell_a = p_static_cast_cell_population->GetCellUsingLocationIndex(nodeA_global_index)->GetCellData()->GetItem("drag");
double drag_cell_b = p_static_cast_cell_population->GetCellUsingLocationIndex(nodeB_global_index)->GetCellData()->GetItem("drag");

            // Add the force contribution to each node
            c_vector<double, SPACE_DIM> scaled_force = force/drag_cell_a;
            c_vector<double, SPACE_DIM> negative_force = -1.0*force/drag_cell_b;
            spring_iterator.GetNodeB()->AddAppliedForceContribution(negative_force);
            spring_iterator.GetNodeA()->AddAppliedForceContribution(scaled_force);
        }
    }
    else    // This is a NodeBasedCellPopulation
    {
        AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>* p_static_cast_cell_population = static_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation);

        std::vector< std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>* > >& r_node_pairs = p_static_cast_cell_population->rGetNodePairs();

        for (typename std::vector< std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>* > >::iterator iter = r_node_pairs.begin();
            iter != r_node_pairs.end();
            iter++)
        {
            std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>* > pair = *iter;

            unsigned node_a_index = pair.first->GetIndex();
            unsigned node_b_index = pair.second->GetIndex();

            // Calculate the force between nodes
            c_vector<double, SPACE_DIM> force = this->CalculateForceBetweenNodes(node_a_index, node_b_index, rCellPopulation);
            for (unsigned j=0; j<SPACE_DIM; j++)
            {
                assert(!std::isnan(force[j]));
            }

double drag_cell_a = p_static_cast_cell_population->GetCellUsingLocationIndex(node_a_index)->GetCellData()->GetItem("drag");
double drag_cell_b = p_static_cast_cell_population->GetCellUsingLocationIndex(node_b_index)->GetCellData()->GetItem("drag");

            // Add the force contribution to each node
            c_vector<double, SPACE_DIM> scaled_force = force/drag_cell_a;
            c_vector<double, SPACE_DIM> negative_force = -1.0*force/drag_cell_b;

            pair.first->AddAppliedForceContribution(scaled_force);
            pair.second->AddAppliedForceContribution(negative_force);
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VariableDampingGeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Call direct parent class
    GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class VariableDampingGeneralisedLinearSpringForce<1,1>;
template class VariableDampingGeneralisedLinearSpringForce<1,2>;
template class VariableDampingGeneralisedLinearSpringForce<2,2>;
template class VariableDampingGeneralisedLinearSpringForce<1,3>;
template class VariableDampingGeneralisedLinearSpringForce<2,3>;
template class VariableDampingGeneralisedLinearSpringForce<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(VariableDampingGeneralisedLinearSpringForce)

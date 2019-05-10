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

#ifndef SimpleCellDataComplexNuclearCellCycleModel_HPP_
#define SimpleCellDataComplexNuclearCellCycleModel_HPP_

#include "AbstractSimplePhaseBasedCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * Simple Complex-dependent cell-cycle model.
 */
class SimpleCellDataComplexNuclearCellCycleModel : public AbstractSimplePhaseBasedCellCycleModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimplePhaseBasedCellCycleModel>(*this);

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;

        archive & mMinimumCellSize;
    }

protected:

    /**
     * Smallest size (in CD^2) which a cell can divide at
     */
    double mMinimumCellSize;

    /**
     * Stochastically set the G1 duration. The G1 duration is taken
     * from a normal distribution whose mean is the G1 duration given
     * in AbstractCellCycleModel for the cell type and whose standard deviation
     * is 1.
     *
     * Called on cell creation at the start of a simulation, and for both
     * parent and daughter cells at cell division.
     */
    void SetG1Duration();

    /**
     * Protected copy-constructor for use by CreateCellCycleModel.
     * The only way for external code to create a copy of a cell cycle model
     * is by calling that method, to ensure that a model of the correct subclass is created.
     * This copy-constructor helps subclasses to ensure that all member variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a daughter cell upon cell division.
     * Note that the parent cell cycle model will have had ResetForDivision() called just before CreateCellCycleModel() is called,
     * so performing an exact copy of the parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     *
     * @param rModel the cell cycle model to copy.
     */
    SimpleCellDataComplexNuclearCellCycleModel(const SimpleCellDataComplexNuclearCellCycleModel& rModel);

public:

    /**
     * Constructor - just a default, mBirthTime is now set in the AbstractCellCycleModel class.
     * mG1Duration is set very high, it is set for the individual cells when InitialiseDaughterCell is called.
     */
    SimpleCellDataComplexNuclearCellCycleModel();

    /*
     * Get the amount of bound complex in the CM
     */
    double GetComplexLevel();

    /**
     * Overridden UpdateCellCyclePhase() method.
     */
    virtual void UpdateCellCyclePhase();

    /**
     * Overridden InitialiseDaughterCell() method.
     */
    virtual void InitialiseDaughterCell();

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     * @return the new cell-cycle model
     */
    virtual AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Overridden CanCellTerminallyDifferentiate() method.
     * @return whether cell can terminally differentiate
     */
    virtual bool CanCellTerminallyDifferentiate();

    /**
     * @return mMinimumCellSize
     */
    double GetMinimumCellSize();

    /**
     * Set mMinimumCellSize.
     *
     * @param minimumCellSize the value of mMinimumCellSize
     */
    void SetMinimumCellSize(double minimumCellSize);

    /**
     * Outputs cell-cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(SimpleCellDataComplexNuclearCellCycleModel)

#endif /*SimpleCellDataComplexNuclearCellCycleModel_HPP_*/

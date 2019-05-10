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

#ifndef WNTCONCENTRATIONMODIFIER_HPP_
#define WNTCONCENTRATIONMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"


// When replacing WntConcentration uncomment the enum
#include "WntConcentration.hpp"
///**
// * Possible types of WntConcentration, currently:
// *  NONE - for testing and to remove Wnt dependence
// *  LINEAR - for cylindrical crypt model
// *  RADIAL - for crypt projection model
// */
//typedef enum WntConcentrationType_
//{
//    NONE,
//    LINEAR,
//    RADIAL,
//    EXPONENTIAL
//} WntConcentrationType;

/**
 * A modifier class which at each simulation time step calculates the volume of each cell
 * and stores it in in the CellData property as "volume". To be used in conjunction with
 * contact inhibition cell cycle models.
 */
template<unsigned DIM>
class WntConcentrationModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
    /**
      * The length of the crypt.
      */
     double mCryptLength;

     /**
      * Whether this WntConcentrationModifier object has had its crypt length set.
      */
     bool mLengthSet;

     /**
      * The type of WntConcentration current options are
      * NONE - returns zero everywhere
      * LINEAR - decreases from 1 to zero at height specified by mWntConcentrationParameter
      * RADIAL - decreases from 1 to zero at height specified by mWntConcentrationParameter
      */
     WntConcentrationType mWntType;

     /**
      * Whether this WntConcentrationModifier object has had its type set.
      */
     bool mTypeSet;

     /**
      * A value to return for testing purposes.
      */
     double mConstantWntValueForTesting;

     /**
      * Whether to return the testing value
      * (when false WntConcentration works with CellPopulation).
      */
     bool mUseConstantWntValueForTesting;

     /**
      * For LINEAR or RADIAL Wnt type:
      * The proportion of the crypt that has a Wnt gradient.
      * The Wnt concentration goes from one at the base to zero at this height up the crypt.
      *
      * For EXPONENTIAL Wnt type:
      * The parameter lambda in the Wnt concentration
      * Wnt = exp(-height/lambda)
      */
     double mWntConcentrationParameter;

     /**
      * Parameter a, for use in crypt projection simulations, in which the crypt
      * surface is given in cylindrical polar coordinates by z = a*r^b.
      * mCryptProjectionParameterA has no units
      */
     double mCryptProjectionParameterA;

     /**
      * Parameter b, for use in crypt projection simulations, in which the crypt
      * surface is given in cylindrical polar coordinates by z = a*r^b.
      * mCryptProjectionParameterB has no units
      */
    double mCryptProjectionParameterB;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
        archive & mCryptLength;
        archive & mLengthSet;
        archive & mWntType;
        archive & mTypeSet;
        archive & mConstantWntValueForTesting;
        archive & mUseConstantWntValueForTesting;
        archive & mWntConcentrationParameter;
        archive & mCryptProjectionParameterA;
        archive & mCryptProjectionParameterB;
    }

public:

    /**
     * Default constructor.
     */
    WntConcentrationModifier();

    /**
     * Destructor.
     */
    virtual ~WntConcentrationModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Helper method to compute the wnt concentration of each cell in the population and store these in the CellData.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);


    /**
     * @return mCryptLength
     */
    double GetCryptLength();

    /**
     * Set mCryptLength. Must be called before GetWntLevel().
     *
     * @param cryptLength  the new value of mCryptLength
     */
    void SetCryptLength(double cryptLength);

    /**
     * @return the type of Wnt concentration.
     */
    WntConcentrationType GetType();

    /**
     * Set the type of Wnt concentration. Must be called before GetWntLevel().
     *
     * @param type the type of Wnt concentration
     */
    void SetType(WntConcentrationType type);

    /**
     * Force the Wnt concentration to return a given value for all cells.
     * Only for testing.
     *
     * @param value the constant value to set the Wnt concentration to be
     */
    void SetConstantWntValueForTesting(double value);

    /**
     * @return mWntConcentrationParameter
     */
    double GetWntConcentrationParameter();

    /**
     * Set mWntConcentrationParameter.
     *
     * @param wntConcentrationParameter the new value of mWntConcentrationParameter
     */
    void SetWntConcentrationParameter(double wntConcentrationParameter);

    /**
     * @return mCryptProjectionParameterA
     */
    double GetCryptProjectionParameterA();

    /**
     * @return mCryptProjectionParameterB
     */
    double GetCryptProjectionParameterB();

    /**
     * Set mCryptProjectionParameterA.
     *
     * @param cryptProjectionParameterA  the new value of mCryptProjectionParameterA
     */
    void SetCryptProjectionParameterA(double cryptProjectionParameterA);

    /**
     * Set mCryptProjectionParameterB.
     *
     * @param cryptProjectionParameterB  the new value of mCryptProjectionParameterB
     */
    void SetCryptProjectionParameterB(double cryptProjectionParameterB);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(WntConcentrationModifier)

#endif /*WNTCONCENTRATIONMODIFIER_HPP_*/

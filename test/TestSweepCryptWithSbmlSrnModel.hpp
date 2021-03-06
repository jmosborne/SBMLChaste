/*
 * Run simulations for ... paper 
 *
 * Author: Axel Almet
 * Created on: April 22 2016
 * Last Modified: Ozzy 16th Nopv 2018
 */

#ifndef TESTCRYPTSIMULATIONSWITHSBMLGENERATEDSRNMODEL_HPP_
#define TESTCRYPTSIMULATIONSWITHSBMLGENERATEDSRNMODEL_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files
#include "CellBasedSimulationArchiver.hpp" //Needed if we would like to save/load simulations
#include "AbstractCellBasedWithTimingsTestSuite.hpp" //Needed for cell-based tests: times simulations, generates random numbers and has cell properties
#include "PetscSetupAndFinalize.hpp"
#include "CheckpointArchiveTypes.hpp" //Needed if we use GetIdentifier() method (which we do)
#include "SmartPointers.hpp" //Enables macros to save typing
#include "HoneycombMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"

#include "VanLeeuwen2007SrnModel.hpp" //This will be generated by the SBML translator
#include "Tan2014SbmlSrnModel.hpp" //This will be generated by the SBML translator

#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "OffLatticeSimulation.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellPopulationMutantHeightWriter.hpp"
#include "CellPopulationProlifHeightWriter.hpp"
#include "VariableDampingGeneralisedLinearSpringForce.hpp"
#include "GeneralisedLinearSpringForce.hpp"

#include "WntConcentrationModifier.hpp" //Modifier to impose Wnt gradient
#include "SimpleCellDataComplexNuclearCellCycleModel.hpp" //Bound beta-catenin-based cell cycle model (similar to Simple Cell Data Wnt CCM)
#//include "BetaCateninConcentrationModifier.hpp"
#include "VolumeTrackingModifier.hpp"
#include "CellVolumesWriter.hpp"
#include "CellLabelWriter.hpp"

/* 2D Cylindrical crypt classes */
#include "CryptCellsGenerator.hpp" //Helper class for cell generation in crypt
#include "CryptSimulation2d.hpp"
#include "SloughingCellKiller.hpp"

#include "CellDataItemWriter.hpp"
#include "CellLabel.hpp"

#include "OutputFileHandler.hpp"
#include "FileComparison.hpp"

#include "FakePetscSetup.hpp"
#include "Debug.hpp"
#include "Exception.hpp"
#include "CommandLineArguments.hpp"

class TestCryptWithSBMLSrnModel : public AbstractCellBasedWithTimingsTestSuite
{
private:

	/* Function to generate cells for 2D crypt */
	void Generate2dCryptCells(std::vector<CellPtr>& rCells, Cylindrical2dMesh* pMesh, const std::vector<unsigned> locationIndices, double cryptLength, bool randomBirthTimes, bool isTan)
	{
	    rCells.clear();

	    RandomNumberGenerator* p_random_num_gen = RandomNumberGenerator::Instance();

        //Define the mutation states and cell proliferative types
		boost::shared_ptr<AbstractCellProperty> p_wildtype_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();
		boost::shared_ptr<AbstractCellProperty> p_transit_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();

	    // Loop over the mesh and populate rCells
	    for (unsigned i=0; i<locationIndices.size(); i++)
	    {

	        // Create a cell-cycle model and set the spatial dimension
	    	SimpleCellDataComplexNuclearCellCycleModel* p_cell_cycle_model = new SimpleCellDataComplexNuclearCellCycleModel;
	        p_cell_cycle_model->SetDimension(2);

			AbstractOdeSrnModel* p_srn_model;
			if(isTan)
			{
				p_srn_model = new Tan2014SbmlSrnModel;
			}
			else
			{
				p_srn_model = new VanLeeuwen2007SrnModel;
			}

	        double birth_time = 0.0;
	        if (randomBirthTimes)
	        {
	            birth_time = -p_random_num_gen->ranf();
	        }

	        // Create a cell
	        CellPtr p_cell(new Cell(p_wildtype_state, p_cell_cycle_model,p_srn_model));

	        // Set the cell's proliferative type,all are transit
	        p_cell->SetCellProliferativeType(p_transit_type);
	        
	        // Set the cell's birth time
	        birth_time *= 12.0; // hours
	        p_cell->SetBirthTime(birth_time);

			if(isTan)
	        {
				p_cell->GetCellData()->SetItem("gamma",1.0);
	        	p_cell->GetCellData()->SetItem("ComplexTransitThreshold",739.9489837338057);
			}
			else
			{
				p_cell->GetCellData()->SetItem("gamma1",1.0);
				p_cell->GetCellData()->SetItem("gamma2",1.0);
				p_cell->GetCellData()->SetItem("ComplexTransitThreshold",11.666666666665009);
			}

	        rCells.push_back(p_cell);
	    }
	}

	// Mutation 1: Change gamma1 or gamma2
	// gamma1 corresponds to the mutation that reduces alpha_d (loss of APC allele in paper)
	// gamma2 corresponds to r7 and r16 being reduced (called `double mutants' in paper)

	void MutateCells(AbstractCellPopulation<2,2>& rCellPopulation, c_vector<double,2> blobCentre , double blobRadius, double mutationParameter, bool isTan, bool isGamma1)
	{

		boost::shared_ptr<AbstractCellProperty> p_label(rCellPopulation.GetCellPropertyRegistry()->Get<CellLabel>());

	    for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
	         cell_iter != rCellPopulation.End();
	         ++cell_iter)
	    {
	    	c_vector<double,2> cell_centre = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

	    	if (norm_2(blobCentre - cell_centre) < blobRadius)
			{
				if(isTan)
				{
					cell_iter->GetCellData()->SetItem("gamma", mutationParameter);
				}
				else
				{
					if(isGamma1)
					{
						cell_iter->GetCellData()->SetItem("gamma1", mutationParameter);
					}
					else
					{
						cell_iter->GetCellData()->SetItem("gamma2", mutationParameter);
					}
				}
				cell_iter->AddCellProperty(p_label);
			}
	    }
	}


	void MutateAllCells(AbstractCellPopulation<2,2>& rCellPopulation, double mutationParameter, bool isTan, bool isGamma1)
	{

	    for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
	         cell_iter != rCellPopulation.End();
	         ++cell_iter)
	    {
			if(isTan)
			{
				cell_iter->GetCellData()->SetItem("gamma", mutationParameter);
			}
			else
			{
				if(isGamma1)
				{
					cell_iter->GetCellData()->SetItem("gamma1", mutationParameter);
				}
				else
				{
					cell_iter->GetCellData()->SetItem("gamma2", mutationParameter);
				}
			}
		}
	}

public:

	/*
	 * Test with 2D Cylindrical Crypt (Van Leeuwen et al.,	 2009)
	 */
	void TestCylindricalCryptWithMutantBlob() throw (Exception)
	{
		//Command line argument stuff: Get the seed parameters
		TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-run_index"));
        unsigned start_index = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-run_index");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-num_runs"));
        unsigned num_runs = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-num_runs");

		TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-is_tan"));
        bool is_tan = CommandLineArguments::Instance()->GetBoolCorrespondingToOption("-is_tan");

	 	TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-is_gamma_1"));
        bool is_gamma_1 = CommandLineArguments::Instance()->GetBoolCorrespondingToOption("-is_gamma_1");

// unsigned start_index = 0;
// unsigned num_runs = 1;
// bool is_tan = true;
// bool is_gamma_1 = false;

		double time_to_steady_state = 500;//100.0;
		double time_post_mutation = 50;//50.0; 
		double time_step = 0.005; 

		//Set simulation parameters
		unsigned cells_across = 10;
		unsigned cells_up = 20;
		unsigned ghost_rows = 2;

		double crypt_height = 20.0;

		double gamma_step = 0.5;
		double all_gamma_step = 0.5;

        // Loop over the random seed.
		for(unsigned sim_index=start_index; sim_index < start_index + num_runs; sim_index++)
		{
			std::cout << " Run number " << sim_index << "... \n" << std::flush;

			// Reseed the random number generator
			RandomNumberGenerator::Instance()->Reseed(sim_index);

			//Define mesh of cells
			CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, ghost_rows);
			Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

			//Obtain vector of non-ghost nodes to generate cell population
			std::vector<unsigned> real_indices = generator.GetCellLocationIndices();

			//Initialise vector of cells
			std::vector<CellPtr> cells;
			Generate2dCryptCells(cells, p_mesh, real_indices, crypt_height, true, is_tan);

			//Create cell population
			MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, real_indices);

			//Explicitly request output for Paraview
			cell_population.AddCellWriter<CellVolumesWriter>();
			cell_population.AddCellWriter<CellLabelWriter>();
			cell_population.AddPopulationWriter<VoronoiDataWriter>();
			cell_population.AddPopulationWriter<CellPopulationMutantHeightWriter>();
			cell_population.AddPopulationWriter<CellPopulationProlifHeightWriter>();
			// Make cell data writer so can pass in variable name
			
			if(is_tan)
			{
				boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer_1(new CellDataItemWriter<2,2>("bcat_cm"));
				cell_population.AddCellWriter(p_cell_data_item_writer_1);
				boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer_2(new CellDataItemWriter<2,2>("ligand_cm"));
				cell_population.AddCellWriter(p_cell_data_item_writer_2);
				boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer_3(new CellDataItemWriter<2,2>("complex_cm"));
				cell_population.AddCellWriter(p_cell_data_item_writer_3);
				boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer_4(new CellDataItemWriter<2,2>("bcat_nu"));
				cell_population.AddCellWriter(p_cell_data_item_writer_4);
				boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer_5(new CellDataItemWriter<2,2>("ligand_nu"));
				cell_population.AddCellWriter(p_cell_data_item_writer_5);
				boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer_6(new CellDataItemWriter<2,2>("complex_nu"));
				cell_population.AddCellWriter(p_cell_data_item_writer_6);
				boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer_10(new CellDataItemWriter<2,2>("gamma"));
				cell_population.AddCellWriter(p_cell_data_item_writer_10);
			}
			else 
			{
				boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer_1(new CellDataItemWriter<2,2>("C_F"));
				cell_population.AddCellWriter(p_cell_data_item_writer_1);
				boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer_2(new CellDataItemWriter<2,2>("C_u"));
				cell_population.AddCellWriter(p_cell_data_item_writer_2);
				boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer_3(new CellDataItemWriter<2,2>("C_A"));
				cell_population.AddCellWriter(p_cell_data_item_writer_3);
				boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer_4(new CellDataItemWriter<2,2>("C_T"));
				cell_population.AddCellWriter(p_cell_data_item_writer_4);

				boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer_10(new CellDataItemWriter<2,2>("gamma1"));
				cell_population.AddCellWriter(p_cell_data_item_writer_10);
				boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer_11(new CellDataItemWriter<2,2>("gamma2"));
				cell_population.AddCellWriter(p_cell_data_item_writer_11);
			}

			boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer_8(new CellDataItemWriter<2,2>("wnt_level"));
			cell_population.AddCellWriter(p_cell_data_item_writer_8);
			boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer_9(new CellDataItemWriter<2,2>("drag"));
			cell_population.AddCellWriter(p_cell_data_item_writer_9);

			//Define simulation class
			CryptSimulation2d simulator(cell_population);
			simulator.UseJiggledBottomCells();

			//Add force class
			//MAKE_PTR(VariableDampingGeneralisedLinearSpringForce<2>, p_spring_force);
			MAKE_PTR(GeneralisedLinearSpringForce<2>, p_spring_force);
			p_spring_force->SetCutOffLength(1.5); // Default value
			p_spring_force->SetMeinekeSpringStiffness(50.0); // Usually 15)
			simulator.AddForce(p_spring_force);

			//Add cell killer
			MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&cell_population, crypt_height));
			simulator.AddCellKiller(p_killer);

			//Add Wnt concentration modifier
			MAKE_PTR(WntConcentrationModifier<2>, p_wnt_modifier);
			p_wnt_modifier->SetType(LINEAR);
			p_wnt_modifier->SetCryptLength(crypt_height);
			simulator.AddSimulationModifier(p_wnt_modifier);

			//Add modifier to track volume
			MAKE_PTR(VolumeTrackingModifier<2>, p_volume_modifier);
			simulator.AddSimulationModifier(p_volume_modifier);

			//Simulation pre-amble
			std::string main_directory = "SbmlCylindricalCrypt/Tan/";
			if (!is_tan)
			{
				if (is_gamma_1)
				{
					main_directory = "SbmlCylindricalCrypt/VanLeeuwenGamma1/";
				}
				else
				{
					main_directory = "SbmlCylindricalCrypt/VanLeeuwenGamma2/";
				}
			}
			std::string steady_state_output_directory, output_directory;

			std::stringstream out;
			out << sim_index;
			steady_state_output_directory = main_directory +  out.str() + "/SteadyState" ; 

			simulator.SetOutputDirectory(steady_state_output_directory);
			simulator.SetSamplingTimestepMultiple(1.0/time_step); //simulation at every hour
			simulator.SetEndTime(time_to_steady_state);
			simulator.SetDt(time_step); //This is very important, otherwise solutions won't be correct
			simulator.Solve();

			// Save simulation in steady state
			CellBasedSimulationArchiver<2, CryptSimulation2d>::Save(&simulator);

			// loop over mutation_parameter 
			for (double mutation_parameter=0; mutation_parameter <= 1.001; mutation_parameter+=gamma_step)
			{
				std::cout << "\n mutation_parameter " << mutation_parameter << ", " << std::flush;
			
				// loop over mutation blob heights 
				for (double height=4; height <= 12; height+=4)
				{
					std::cout << "Height " << height << "... " << std::flush;
				
					std::stringstream out;
					out << sim_index << "/MutationParameter_"<< mutation_parameter << "/Height_" << height;
					
					output_directory = main_directory +  out.str(); 

					CryptSimulation2d* p_simulator = CellBasedSimulationArchiver<2, CryptSimulation2d>::Load(steady_state_output_directory,time_to_steady_state);

					// now pick a selection of cells to mutate
					c_vector<double,2> blob_centre;
					blob_centre(0) = (double) cells_across / 2.0;
					blob_centre(1) = height;
					double blob_radius = 2.0;
					MutateCells(p_simulator->rGetCellPopulation(),blob_centre,blob_radius,mutation_parameter,is_tan,is_gamma_1);

					p_simulator->SetOutputDirectory(output_directory);
					p_simulator->SetEndTime(time_post_mutation+time_to_steady_state);
					p_simulator->Solve();

					delete p_simulator;
				}
			}

			for (double mutation_parameter=0; mutation_parameter <= 1.001; mutation_parameter+=all_gamma_step)
			{
				std::cout << "\n mutation_parameter " << mutation_parameter << ", " << std::flush;
			
				
				std::stringstream out;
				out << sim_index << "/AllCellsMutated_"<< mutation_parameter;
					
				output_directory = main_directory +  out.str(); 

				CryptSimulation2d* p_simulator = CellBasedSimulationArchiver<2, CryptSimulation2d>::Load(steady_state_output_directory,time_to_steady_state);


				// now mutate all cells
				MutateAllCells(p_simulator->rGetCellPopulation(),mutation_parameter,is_tan,is_gamma_1);

				p_simulator->SetOutputDirectory(output_directory);
				p_simulator->SetEndTime(time_post_mutation+time_to_steady_state);
				p_simulator->Solve();

				delete p_simulator;
			}

			// Extra Gubbins to get to loop: this is usually done by the SetUp and TearDown methods
			SimulationTime::Instance()->Destroy();
			SimulationTime::Instance()->SetStartTime(0.0);
		}
	}
};
#endif /* TESTCRYPTSIMULATIONSWITHSBMLGENERATEDSRNMODEL_HPP_ */

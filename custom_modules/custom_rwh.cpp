/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

// declare cell definitions here 

Cell_Definition CSC;
Cell_Definition DCC; 

// temp for testing 
int diff_count = 0; 
int div_count = 0; 

void differentiation_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	// am I at the start of a cycle?
	if( phenotype.cycle.data.elapsed_time_in_phase < 0.1 )
	{
//		#pragma omp critical 
//		{ div_count++; } 
	
		if( UniformRandom() <= parameters.doubles( "differentiation_probability" ) )
		{
//			#pragma omp critical 
//			{ diff_count++; } 
			pCell->convert_to_cell_definition( DCC ); 
		}
		
	}

	return; 
} 

void dummy_function( Cell* pCell , Phenotype& phenotype , double dt )
{
	// std::cout << phenotype.cycle.data.transition_rate(0,0) << std::endl; 
}

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// Name the default cell type 
	
	cell_defaults.type = 0; 
	cell_defaults.name = "tumor cell"; 
	
	// set default cell cycle model 

	live.phase_link(0,0).fixed_duration = parameters.bools("fixed_cycle_duration"); 
	cell_defaults.functions.cycle_model = live; 
	
	// set default_cell_functions; 
	
	cell_defaults.functions.update_phenotype = NULL; 
	
	// needed for a 2-D simulation: 
	
	/* grab code from heterogeneity */ 
	
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	// make sure the defaults are self-consistent. 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

	// set the rate terms in the default phenotype 

	// first find index for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 

	int G0G1_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::G0G1_phase );
	int S_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::S_phase );
	int live_index = live.find_phase_index( PhysiCell_constants::live ); 

	// initially no necrosis 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 
	
	// initially no apoptosis 
	cell_defaults.phenotype.death.rates[apoptosis_model_index] = 0.0; 

	// set oxygen uptake / secretion parameters for the default cell type 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 10; 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0; 
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 38; 
	
	// add custom data here, if any 
	

	// Now, let's define another cell type. 
	// It's best to just copy the default and modify it. 
	
	// define CSC type 
	
	CSC = cell_defaults; 
	CSC.type = 1;
	// set birth rate 
	CSC.phenotype.cycle.data.transition_rate(live_index,live_index) = 
		parameters.doubles( "CSC_birth_rate" );  
	// set death rate 
	CSC.phenotype.death.rates[apoptosis_model_index] = 
		parameters.doubles( "CSC_apoptosis_rate" ); 
	// set the phenotype function to non-NULL
	CSC.functions.update_phenotype = differentiation_function; 
	
	// define DCC type  
	
	DCC = cell_defaults; 
	DCC.type = 2; 	
	// set birth rate 
	DCC.phenotype.cycle.data.transition_rate(live_index,live_index) = 
		parameters.doubles( "DCC_birth_rate" );  
	// set death rate 
	DCC.phenotype.death.rates[apoptosis_model_index] = 
		parameters.doubles( "DCC_apoptosis_rate" ); 

/*	
	// make this cell type randomly motile, less adhesive, greater survival, 
	// and less proliferative 
	
	motile_cell = cell_defaults; 
	motile_cell.type = 1; 
	motile_cell.name = "motile tumor cell"; 
	
	// make sure the new cell type has its own reference phenotype
	
	motile_cell.parameters.pReference_live_phenotype = &( motile_cell.phenotype ); 
	
	// enable random motility 
	motile_cell.phenotype.motility.is_motile = true; 
	motile_cell.phenotype.motility.persistence_time = parameters.doubles( "motile_cell_persistence_time" ); // 15.0; 
	motile_cell.phenotype.motility.migration_speed = parameters.doubles( "motile_cell_migration_speed" ); // 0.25 micron/minute 
	motile_cell.phenotype.motility.migration_bias = 0.0;// completely random 
	
	// Set cell-cell adhesion to 5% of other cells 
	motile_cell.phenotype.mechanics.cell_cell_adhesion_strength *= parameters.doubles( "motile_cell_relative_adhesion" ); // 0.05; 
	
	// Set apoptosis to zero 
	motile_cell.phenotype.death.rates[apoptosis_model_index] = parameters.doubles( "motile_cell_apoptosis_rate" ); // 0.0; 
	
	// Set proliferation to 10% of other cells. 
	// Alter the transition rate from G0G1 state to S state
	motile_cell.phenotype.cycle.data.transition_rate(G0G1_index,S_index) *= 
		parameters.doubles( "motile_cell_relative_cycle_entry_rate" ); // 0.1; 
*/
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
/* now this is in XML 
	default_microenvironment_options.X_range = {-1000, 1000}; 
	default_microenvironment_options.Y_range = {-1000, 1000}; 
	default_microenvironment_options.simulate_2D = true; 
*/
	
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
/* now this is in XML 	
	// no gradients need for this example 

	default_microenvironment_options.calculate_gradients = false; 
	
	// set Dirichlet conditions 

	default_microenvironment_options.outer_Dirichlet_conditions = true;
	
	// if there are more substrates, resize accordingly 
	std::vector<double> bc_vector( 1 , 38.0 ); // 5% o2
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
	
	// set initial conditions 
	default_microenvironment_options.initial_condition_vector = { 38.0 }; 
*/
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	// create some cells near the origin
	
	Cell* pC;
	
	double Lx = default_microenvironment_options.X_range[1] - default_microenvironment_options.X_range[0]; 
	double Ly = default_microenvironment_options.Y_range[1] - default_microenvironment_options.Y_range[0]; 
	
	for( int n = 0 ; n < parameters.ints( "number_of_stem_cells" ) ; n++ )
	{
		double x = default_microenvironment_options.X_range[0] + 0.1*Lx + 0.8*UniformRandom()*Lx; 
		double y = default_microenvironment_options.Y_range[0] + 0.1*Ly + 0.8*UniformRandom()*Ly; 
		
		pC = create_cell( CSC ); 
		pC->assign_position( x,y, 0.0 );
		
	}
	
	for( int n = 0 ; n < parameters.ints( "number_of_differentiated_cells" ) ; n++ )
	{
		double x = default_microenvironment_options.X_range[0] + 0.1*Lx + 0.8*UniformRandom()*Lx; 
		double y = default_microenvironment_options.Y_range[0] + 0.1*Ly + 0.8*UniformRandom()*Ly; 
		
		pC = create_cell( DCC ); 
		pC->assign_position( x,y, 0.0 );
		
	}
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with flow cytometry coloring 
	
	std::vector<std::string> output = { "black", "black", "black", "black" } ; // simple_cell_coloring(pCell); 
		
	std::cout << "----- time = " << PhysiCell_globals.current_time << std::endl;
	if( pCell->phenotype.death.dead == false && pCell->type == 1 )
	{
		 output[0] = "red"; 
		 output[2] = "red"; 
		 return output; 
	}
	
	if( pCell->phenotype.death.dead == false && pCell->type == 2 )
	{
		 output[0] = "cyan"; 
		 output[2] = "cyan"; 
		 return output; 
	}
	
	if( pCell->phenotype.death.dead == true )
	{
		 output[0] = "black"; 
		 output[2] = "black"; 
		 return output; 
	}
	
	return output; 
}

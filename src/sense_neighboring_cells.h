// -----------------------------------------------------------------------------
//
// Copyright (C) The BioDynaMo Project.
// All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
//
// See the LICENSE file distributed with this work for details.
// See the NOTICE file distributed with this work for additional information
// regarding copyright ownership.
//
// -----------------------------------------------------------------------------
#ifndef SENSE_NEIGHBORING_CELLS_H_
#define SENSE_NEIGHBORING_CELLS_H_

#include "biodynamo.h"


#include "GRNCellObject.h"
#include "GRNCellModule.h"


// set numer of simulation steps
const int simulation_steps = 2500; // Time between two simulation steps equals: 0.01hours (default)

//const double x_range = 150, y_range = 150, z_range = 4500; // set dims of simulation space
//const double simulation_cube_dim = std::max(x_range, y_range);
const double simulation_cube_dim = 100;
const double z_pos_precursor = 0 - (simulation_cube_dim/2); // (x/y range + z position)

// number of precursor cells
const size_t num_precursor_cells = 2;  // number of precursor cells (S1) in the simulation

namespace bdm {

// Define compile time parameter
    BDM_CTPARAM() { BDM_CTPARAM_HEADER();
        // add GRNObject to simulation
        using SimObjectTypes = CTList<GRNCell>;
        // add GRNModule to simulation
        BDM_CTPARAM_FOR(bdm, GRNCell) { using BiologyModules = CTList<GRNModule>; };
    };

inline int Simulate(int argc, const char** argv) {
    // set space parameters of the simulation
    auto set_param = [](auto* param) {
        param->bound_space_ = true;
        param->min_bound_ = -(simulation_cube_dim/2);
        param->max_bound_ = (simulation_cube_dim/2);
        param->run_mechanical_interactions_ = true;
    };

    Simulation<> simulation(argc, argv, set_param);
    auto* rm = simulation.GetResourceManager();  // get pointer to resource manager
    auto* random = simulation.GetRandom();  // get thread of local random number generator.

    // Since sim_objects in this simulation won't modify neighbors, we can
    // safely disable neighbor guards to improve performance.
    simulation.GetExecutionContext()->DisableNeighborGuard();

    double x_coord, y_coord, z_coord;

    // 2D plate for precursor cells (150x150)
    double x_min = 0 - (simulation_cube_dim/2);  // set position of the plate with (0,0) at the center of the simulation space
    double x_max = 0 + (simulation_cube_dim/2);
    double y_min = 0 - (simulation_cube_dim/2);
    double y_max = 0 + (simulation_cube_dim/2);

    // allocate the correct number of cell in our cells structure before
    // cell creation
    rm->template Reserve<GRNCell>(num_precursor_cells);

    // create 2d Layer of cells
    for (size_t i = 0; i < num_precursor_cells; ++i) {
        // create coordinates for cells in 2D plate
        x_coord = random->Uniform(x_min, x_max);
        y_coord = random->Uniform(y_min, y_max);
        z_coord = z_pos_precursor;

        // creating the cell at position x, y, z
        GRNCell cell({x_coord, y_coord, z_coord});
        // set cell parameters
        cell.SetDiameter(default_cell_diameter);
        cell.SetAdherence(0.0001);
        cell.SetMass(0.1);
        cell.AddBiologyModule(GRNModule());
        rm->push_back(cell);// put the created cell in our cells structure
    }

    // 4. Run simulation for N timesteps
    simulation.GetScheduler()->Simulate(simulation_steps);

    std::cout << "Simulation completed successfully!" << std::endl;
    return 0;
}

}  // namespace bdm

#endif  // SENSE_NEIGHBORING_CELLS_H_

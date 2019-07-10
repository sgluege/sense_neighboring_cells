//
// Created by Stefan Glüge on 2019-02-07.
//
// GRN Cell Module, defines the cell behaviour in diffrent states and when/how to switch staes in the GRN

#ifndef LAYERFORMATION_C_GRNCELLMODULE_H
#define LAYERFORMATION_C_GRNCELLMODULE_H

//#include "biodynamo.h"
const double default_cell_diameter = 6; // default diameter of new cells

namespace bdm {

struct GRNModule : public BaseBiologyModule {
    BDM_STATELESS_BM_HEADER(GRNModule, BaseBiologyModule, 1);

public:
    GRNModule() : BaseBiologyModule(gAllEventIds) {}
    /// a single Cell Module that shall reproduce the behaviour of cells in different stages in the GRN
    /// Empty default event constructor, because GrowthModule does not have state.
    template <typename TEvent, typename TBm>
    GRNModule(const TEvent& event, TBm* other, uint64_t new_oid = 0)
            : BaseBiologyModule(event, other, new_oid) {}

    void Run(SimObject* so) override {
        if (auto *cell = dynamic_cast<GRNCell *>(so)) {
            // read current cell types
            const std::string &current_cell_type = cell->GetCellType();

            // run cell behaviour depending on cell type
            if (current_cell_type == "S1") {
                typeS1behaviour(cell);
            } else {
                std::cout << "no behaviour defined for cell of type: " << current_cell_type << std::endl;
            }
        }
    }

    // define behaviour of cell type S1
    template <typename T, typename TSimulation = Simulation>
    void typeS1behaviour(T* cell){

        // find number of neighboring cells in a certain radius
        const double sense_radius = 100.2;
        int same_type = 0;
        int other_type = 0;

//        std::string act_type = cell->GetCellType();
//        std::cout << "act cell type: " << act_type << std::endl;
//        std::string some_type = "S1";
//
//        if (act_type == some_type) {
//            std::cout << act_type << "==" << some_type << std::endl;
//        }

        // lambda updating counters for cell neighbors
        auto countNeighbours = [&](const auto* neighbor) {
            // if neighbor is a GRNCell
            if (neighbor->template IsSoType<GRNCell>()) {
                auto n_soptr = neighbor->template
                        ReinterpretCast<GRNCell>()->GetSoPtr();
                // if GRNCell have the same type
                if (n_soptr->GetCellType() == cell->GetCellType()) {
                    same_type++;
                } else {
                    other_type++;
                }
            }
        }; // end lambda
//
        auto* ctxt = TSimulation::GetActive()->GetExecutionContext();
        ctxt->ForEachNeighborWithinRadius(countNeighbours, *cell, sense_radius);

//        std::cout << "same_type : " << same_type << std::endl;
//        std::cout << "other_type : " << other_type << std::endl;


    // cell differentiation
//        int sameType = 6;
        runCellCycleDiffStepS1(cell, same_type); // grow cell

    }

    /// define apoptosis 1 behaviour
    template <typename T, typename TSimulation = Simulation>
    void typeA1behaviour(T* cell){
        cell->RemoveFromSimulation();
    }

    /// run cell differentiation step (used in type S1)
    /// grow cell and divide when max_diam is reached
    template <typename T, typename TSimulation = Simulation>
    void runCellCycleDiffStepS1(T* cell, int same_count){
        if (cell->GetDiameter() < cell->GetMaxDiam()) { // max_diam not reached
            cell->SetDiameter(cell->GetDiameter() + 0.2); // grow cell diameter by adding some number
        } else { // max diam reached -> divide cell
            if (same_count > 5){  // check if there are too many other cells around
                std::cout << "kill cell due to too many other cells around" << std::endl;
                typeA1behaviour(cell);
            } else {
                auto &&daughter = cell->Divide();
                std::cout << "can divide as sameCount is: " << same_count << std::endl;

                daughter->SetDiameter(default_cell_diameter);  // init daughter with default diameter
            }
        }
    }
};

}

#endif //LAYERFORMATION_C_GRNCELLMODULE_H

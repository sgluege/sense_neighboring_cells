//
// Created by Stefan Glüge on 2019-02-07.
//
// GRN Cell Module, defines the cell behaviour in diffrent states and when/how to switch staes in the GRN

#ifndef LAYERFORMATION_C_GRNCELLMODULE_H
#define LAYERFORMATION_C_GRNCELLMODULE_H

//#include "biodynamo.h"
const double counter_sub_quant_threshold = 50; // threshold of intracellular substance ... when to switch biomodules
const double default_cell_diameter = 6; // default diameter of new cells

namespace bdm {

/// Define general GRNModule
struct GRNModule : public BaseBiologyModule {
    GRNModule() : BaseBiologyModule(gAllEventIds) {}
    /// a single Cell Module that shall reproduce the behaviour of celles in different stages in the GRN

    /// Empty default event constructor, because GrowthModule does not have state.
    template <typename TEvent, typename TBm>
    GRNModule(const TEvent& event, TBm* other, uint64_t new_oid = 0)
            : BaseBiologyModule(event, other, new_oid) {}

    template <typename T, typename TSimulation = Simulation<>>
    void Run(T* cell) {
       // read current cell types
       const std::string& current_cell_type = cell->GetCellType();

        // run cell behaviour depending on cell type
        if (current_cell_type == "S1"){
            typeS1behaviour(cell);
        } else {
            std::cout << "no behaviour defined for cell of type: " << current_cell_type << std::endl;
        }
    }

    // define behaviour of cell type S1
    template <typename T, typename TSimulation = Simulation<>>
    void typeS1behaviour(T* cell){

        // find number of neighboring cells in a certain radius
        const double sense_radius = 1.2;
        int sameType = 0;
        int otherType = 0;

        // lambda updating counters for cell neighbors
        auto countNeighbours = [&](const auto* neighbor) {
            // if neighbor is a GRNCell
            if (neighbor->template IsSoType<GRNCell>()) {
                auto n_soptr = neighbor->template
                        ReinterpretCast<GRNCell>()->GetSoPtr();
                // if GRNCell have not the same type
                if (!(n_soptr->GetCellType() == cell->GetCellType())) {
                    // if GRNCell got the same type
                    if (n_soptr->GetCellType() == cell->GetCellType()) {
                        sameType++;
                    }
                    else {
                        otherType++;
                    }
                }
                else {
                    sameType--;
                }
            }
        }; // end lambda

        auto* ctxt = TSimulation::GetActive()->GetExecutionContext();
//        auto* sim = TSimulation::GetActive();
        ctxt->ForEachNeighborWithinRadius(countNeighbours, *cell, sense_radius);

//        std::cout << "#same neighbors: " << sameType << std::endl;
//        std::cout << "#other neighbors: " << otherType << std::endl;

        // cell differentiation
        runCellCycleDiffStepS1(cell, sameType); // grow cell

    }

    /// define apoptosis 1 behaviour
    template <typename T, typename TSimulation = Simulation<>>
    void typeA1behaviour(T* cell){
        cell->RemoveFromSimulation();
    }

    /// run cell differentiation step (used in type S1)
    /// grow cell and divide when max_diam is reached
    template <typename T, typename TSimulation = Simulation<>>
    void runCellCycleDiffStepS1(T* cell, int sameCount){
        if (cell->GetDiameter() < cell->GetMaxDiam()) { // max_diam not reached
            cell->SetDiameter(cell->GetDiameter() + 0.2); // grow cell diameter by adding some number
        } else { // max diam reached -> divide cell
            if (sameCount > 10){  // check if there are too many other cells around
                std::cout << "kill cell due to too many other cells around" << std::endl;
                typeA1behaviour(cell);
            } else {
                auto &&daughter = cell->Divide();
                daughter->SetDiameter(default_cell_diameter);  // init daughter with default diameter
            }
        }
    }


private:
    BDM_CLASS_DEF_NV(GRNModule, 1);
};

}

#endif //LAYERFORMATION_C_GRNCELLMODULE_H

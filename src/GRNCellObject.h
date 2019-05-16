//
// Created by Stefan Glüge on 2019-02-08.
//
// BDM_SIM_OBJECT for the GRNCell

#ifndef GRNCELLOBJECT_H_
#define GRNCELLOBJECT_H_

//#include "biodynamo.h"

// set parameters for initial cells (type S1)
const double intracellular_substance_quantity = 100; // const double intracellular_substance_quantity = 100; //
const double s1max_diam = 10;
const double s1sense_radius = 1.2;
const std::string s1type = "S1";
//const double counter_sub_quant_threshold = 5; // threshold of intracellular substance ... when to switch biomodules
//const double default_cell_diameter = 6; // default diameter of new cells

namespace bdm {

/// Define our new GRN cell (BDM_SIM_OBJECT) which extends Cell by adding extra data
BDM_SIM_OBJECT(GRNCell, Cell) {  // our object extends the Cell object
// create the header with our new data member
BDM_SIM_OBJECT_HEADER(GRNCell, Cell, 1, cell_type_, cell_max_diam_, cell_color_);

public:
    GRNCellExt() {}
    explicit GRNCellExt(const std::array<double, 3>& position) : Base(position) {
        // init cell with type S1, color, max_diam, and intracellular substance quantity
        cell_type_[kIdx] = s1type;
        cell_max_diam_[kIdx] = s1max_diam;
        cell_color_[kIdx] = 0;
    }

    /// If MyCell divides, daughter 2 copies the data members from the mother
    template <typename TMother>
    GRNCellExt(const CellDivisionEvent& event, TMother* mother) : Base(event, mother) {
        cell_type_[kIdx] = mother->cell_type_[mother->kIdx];
        cell_color_[kIdx] = mother->cell_color_[mother->kIdx];
    }

    /// If a cell divides, daughter keeps the same state from its mother.
    template <typename TDaughter>
    void EventHandler(const CellDivisionEvent& event, TDaughter* daughter) {
        Base::EventHandler(event, daughter);
    }

    // getter and setter for our new data member
    void SetCellType(const std::string& cell_type) { cell_type_[kIdx] = std::string(cell_type); }
    const std::string& GetCellType() const { return cell_type_[kIdx]; }

    void SetCellColor(int cell_color) { cell_color_[kIdx] = cell_color; }
    int GetCellColor() const { return cell_color_[kIdx]; }

    void SetMaxDiam(double max_diam) { cell_max_diam_[kIdx] = max_diam; }
    double GetMaxDiam() const { return cell_max_diam_[kIdx]; }

private:
    // declare new data member and define their type
    // private data can only be accessed by public function and not directly
    vec<std::string> cell_type_;
    vec<double> cell_max_diam_;
    vec<int> cell_color_;
};

}  // namespace bdm


#endif

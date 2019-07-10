//
// Created by Stefan Glüge on 2019-02-08.
//
// BDM_SIM_OBJECT for the GRNCell

#ifndef GRNCELLOBJECT_H_
#define GRNCELLOBJECT_H_

//#include "biodynamo.h"

// set parameters for initial cells (type S1)
//const double intracellular_substance_quantity = 100; // const double intracellular_substance_quantity = 100; //
const double s1max_diam = 10;
const std::string s1_type = "S1";
//const double counter_sub_quant_threshold = 5; // threshold of intracellular substance ... when to switch biomodules
//const double default_cell_diameter = 6; // default diameter of new cells

namespace bdm {

/// Define our new GRN cell (BDM_SIM_OBJECT) which extends Cell by adding extra data
class GRNCell : public Cell { // our object extends the Cell object
// create the header with our new data member
    BDM_SIM_OBJECT_HEADER(GRNCell, Cell, 1, cell_type_, cell_max_diam_, cell_color_);

    public:
        GRNCell() {}
        explicit GRNCell(const Double3& position) : Base(position) {
        // init cell with type S1, color, max_diam, and intracellular substance quantity
        cell_type_ = s1_type;
        cell_max_diam_ = s1max_diam;
        cell_color_ = 0;
    }
    /// If MyCell divides, daughter 2 copies the data members from the mother
    GRNCell(const Event& event, SimObject* other, uint64_t new_oid = 0)
            : Base(event, other, new_oid) {
        if (auto *mother = dynamic_cast<GRNCell *>(other)) {
            cell_color_ = mother->cell_color_;
            cell_type_ = mother->cell_type_;
            cell_max_diam_ = mother->cell_max_diam_;
        }
    }

    /// If a cell divides, daughter keeps the same state from its mother.
    void EventHandler(const Event& event, SimObject* other1,
                      SimObject* other2 = nullptr) override {
        Base::EventHandler(event, other1, other2);
    }

    // getter and setter for our new data member
    void SetCellType(const std::string& cell_type) { cell_type_ = std::string(cell_type); }
    const std::string& GetCellType() const { return cell_type_; }

    void SetCellColor(int cell_color) { cell_color_ = cell_color; }
    int GetCellColor() const { return cell_color_; }

    void SetMaxDiam(double max_diam) { cell_max_diam_ = max_diam; }
    double GetMaxDiam() const { return cell_max_diam_; }

private:
    // declare new data member and define their type
    // private data can only be accessed by public function and not directly
    std::string cell_type_;
    double cell_max_diam_;
    int cell_color_;
};

}  // namespace bdm


#endif

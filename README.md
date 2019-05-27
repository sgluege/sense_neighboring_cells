# sense_neighboring_cells
Implement cells that can sense their local environment, and based on the neighboring cell types commit to apoptosis.


Simulation is initialised with cells of type "S1". They grow and divide as long as less than 5 cells of the same type are in a radius of 100.2.


`typeS1behaviour`() in `GRNCellModule.h` uses the `countNeighbours()` lambda function to count cells of the same (and other) types for each existing cell within a `sense_radius = 100.2`.

`runCellCycleDiffStepS1()` than grows the cell unti a max. diameter. It is then divided or killed depending on the number of neighbors with the same type within the radious.

---------------

Biodynamo version:
```
BioDynaMo v0.1.0-152-ge3d39e63
```
Introduction

HydrologyModel is a C++ tool for automatic depression (sink) detection and hierarchical tree construction from Digital Elevation Models (DEM). It implements a Barnes‑based depression filling algorithm that identifies depressions without modifying the original DEM. The algorithm outputs:

Depression ID raster (each depression gets a unique ID)

D8 flow direction raster

Catchment (contributing area) raster

A hierarchical tree structure (HashWeightTree) storing depression relationships (parent‑child, upstream/downstream)

Optional: Outlet elevation DEM and outlet ID raster

Unlike traditional depression filling algorithms that modify elevations, this tool only marks depressions and builds a topological tree, preserving the original topography for subsequent hydrological analysis.

Key Features
No DEM modification – Original elevations are kept intact

Efficient Mas algorithm using a hybrid priority queue (hash heap + plain queue)

Multi‑level depression hierarchy – handles nested depressions

Rich output formats – GeoTIFF rasters + binary tree file

Query module – browse depression properties, upstream/downstream relations, outlet coordinates (with geographic or projected coordinates)

Export module – generate outlet elevation DEM and outlet ID raster

GDAL integration – reads/writes GeoTIFF with full georeferencing support

Output Files
After running Module 1, the following files are created in the output directory:

File	Type	Description
depression_id.tif	GeoTIFF (UInt32)	Each cell contains the depression ID (0 = no depression)
flow_direction.tif	GeoTIFF (Byte)	D8 flow direction codes (1,2,4,8,16,32,64,128)
catchment_areas.tif	GeoTIFF (UInt32)	Catchment ID for cells that drain into a depression outlet
hash_weight_tree.bin	Binary	Hierarchical tree structure (see class HashWeightTree)
Module 3 produces two additional rasters:

outlet_dem.tif – Elevation of depression outlet cells (original elevation)

outlet_dem_id.tif – Depression ID at each outlet cell

Dependencies
CMake (≥ 3.10)

C++17 compiler (g++ 7+ or clang 6+)

GDAL (libgdal-dev) – GeoTIFF I/O

On Ubuntu/Debian:

bash
sudo apt update
sudo apt install -y build-essential cmake libgdal-dev
Building
bash
cd HydrologyModel
mkdir build && cd build
cmake ..
make -j$(nproc)
The executable hydro_model will be created in the build/ directory.

Usage
Prepare your DEM as a GeoTIFF file (single band, 32‑bit float).

Edit main.cpp to set the correct input_file and base_output_dir paths.

Run the executable:

bash
./hydro_model
Choose a module when prompted:

1 – Run depression analysis (Module 1)

2 – Query an existing tree (Module 2)

3 – Export outlet DEM raster (Module 3)

Example
cpp
// In main.cpp, modify these lines:
std::string input_file = "/path/to/your/dem.tif";
std::string base_output_dir = "/path/to/output";
Then:

bash
./hydro_model

Module Details
Module 1 – Depression Analysis
Reads input DEM

Executes Mas algorithm using a hybrid queue

Outputs all four raster files + binary tree

Prints timing and depression count

Module 2 – Tree Query
Loads previously saved hash_weight_tree.bin

Asks for a depression ID

Displays:

Outlet row/col and map coordinates (longitude/latitude or easting/northing)

Outlet elevation

Parent (downstream) depression ID

Level in the tree

Direct upstream/downstream lists (limited to 10 items)

All upstream/downstream depressions (full list)

Module 3 – Export Outlet DEM
Creates a GeoTIFF where only outlet cells keep their original elevation (others = NoData)

Creates a second GeoTIFF with depression ID at each outlet cell

Useful for visualizing outlet locations or as input for further routing models

Algorithm Overview
The Mas algorithm works as follows:

Initialize a priority queue with all boundary cells.

Pop the cell with lowest elevation from the priority queue.

For each neighbor:

If neighbor is lower → part of the same depression (push to plain queue).

If neighbor is higher or equal → set flow direction (D8) and push to priority queue.

When a new depression is discovered (a cell lower than its processed neighbor), a new depression ID is assigned, and its outlet is recorded.

The algorithm continues until all cells are processed.

Finally, the depression hierarchy is built by linking each depression to its downstream parent.

No elevation values are modified during this process.

Performance
Time complexity: O(N log N) with N = number of cells

Memory: approx. 5×N × 4 bytes + overhead (DEM, depression ID, flow direction, catchment, processed flags)

Tested on DEMs up to 10,000×10,000 cells (100 million cells) with reasonable RAM usage.

License
This project is released under the MIT License. See LICENSE file for details.

Contact & Contributing
Issues and pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

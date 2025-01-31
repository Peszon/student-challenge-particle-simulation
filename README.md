# Solution to student challenge particle simulation
## Background

This repository contains details for the coding challenge announced during
student fair FARM at Chalmers, November 2024. In the challenge, the goal 
was to find the amount of colliding particles in a provided 130 000 particle 
coordinate file.

## Compilation Instructions

To compile the project, ensure you have `g++` with C++20 and OpenMP support 
installed. Then, use the following commands:

```sh
make
```

This will compile the project using optimization flags (-O3) and include 
necessary headers. The compiled executable will be available in the the 
./bin directory.

## Solution details

This implementation provides a collision detection framework that reads a dataset of particle positions and identifies particle pairs within a specified distance. The system supports multiple algorithms, with the Spatial Subdivision Parallel algorithm being the default. The framework executes the algorithm multiple times, averaging execution time, and reports the number of detected collisions. By leveraging OpenMP for parallel execution, this solution enhances computational efficiency for large datasets.

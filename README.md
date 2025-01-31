# Solution to student challenge particle simulation

This repository provides a particle collision detection framework with both parallel and sequential algorithms.

## Background

This repository contains details for a solution to the IPS coding challenge announced during
student fair FARM at Chalmers, November 2024. In the challenge, the goal was to detect the number of collisions among 130000 particles based on their positions in 3D space, using an efficient algorithm.

## Compilation Instructions

To compile the project, ensure you have `g++` with C++20 and OpenMP support 
installed. Then, use the following commands:

```sh
make -j16
```

This will compile the project using optimization flags `-O3`, multithreading flags `-j16` and include 
necessary headers. Change the number after -j to adjust the number of cores that will be used to execute the program. The compiled executable will be available in the the ./bin directory and can be executed with:

```sh
./bin/output
``` 

## Solution details

This implementation provides a collision detection testing framework that reads a dataset of particle positions and let you pick the setting for the algorithm testing. The system supports multiple algorithms, with the parallel Spatial subdivision algorithm being the default if you don't change the settings in the `main.cpp` file. The framework executes the algorithm ten times by default, averaging the execution time, and then reports the number of detected collisions and average runtime.

The parallelization of the algorithm using openMP reduces the relative time the algorithm takes by between 0.4 and 0.25. The execution of the program with the parallel spatial subdivision takes on average 0.17 s on my laptop. 

Both the parallel and sequential Spatial subdivision algorithm builds upon the same three steps. First, divide the volume into smaller cells and then create a object for each particle that reside in each cell. Since each particle has a volume (radius = collision distance/2), each particle can reside in multiple cells and thus there can be multiple objects for each particle. Each object has a "cell id" (hashed value of the cell coordinates), "particle_id"(unique integer of the particle) and "home_cell" (1 if the cell is the home cell of the particle and 0 otherwise).
 
The second step is to sort the objects by cell_id so that all objects from the same cell appear consecutively. This is done by a 8 bit radix sort, since each cell_id contains a 48 bit value this will result in only 6 passes of the data in the sequential algorithm. The bit value of the Radix sort was found empirically. 

The third step is to calculate the amount of collisions in each cell and then sum these to find the total number of collisions. Due to the stability attribute of the radix sort and the saved home_cell information it is further possible to reduce the amount of collision checks that need to be performed, which speeds up the program.

## Solution limitations
Since the algorithm depends on three 16 bit integer value to hold the cell_id, and the cell size depends on the collision distance (to limit the amount of cells each particle can reside in). There is a lower limit on the collision distance to ensure unqiue hash value for the cell_ids, meaning that:

```sh
	(max_x_coordinate - min_x_coordinate) / collision_dist < 2^16-1 = 65535
``` 

This is not any problem with regards to the given data and collision distance, but could become a problem if the trial data or trial collision distance is quite different from the test data or test collision distance.

## CUDA implementation

I have tried to implement a CUDA version of the code, but since I don't have access to a machine with a GPU this was quite an uphill battle. The experimental, not-working, CUDA implementation of the code can be found in the CUDA folder. I found GPU accelerated computation really interesting and I hope to further my knowledge in the area during the summer. 

Thank you for considering my code!

Allt gott

Felix Persson

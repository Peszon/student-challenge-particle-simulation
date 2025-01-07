#pragma once

#include <vector>
#include <cmath>

/* **********************************************
 *
 * Calculates the distance between two points in 3D.
 * 
 * Input: 
 *     - pos1: Postion of particle 1 in 3D.
 *     - pos2: Postion of particle 2 in 3D.
 * 
 * Returns: the distance between the particles as a double.
 *
 * **********************************************/
double calculate_distance(const std::vector<double>& pos1, const std::vector<double>& pos2);

/* **********************************************
 *
 * Calculates the number of particle pairs in the
 * given list that is within the "collision_distance"
 * of eachother.
 * 
 * Input:
 *      - coordinates: 2D list of the 3D positions of the particles.
 *      - collision_distance: Distance that causes a collision.
 *  
 * Returns: the number of collisions as a int.
 *
 * **********************************************/
int calculate_number_colisions(const std::vector<std::vector<double>>& coordinates, double collision_distance);

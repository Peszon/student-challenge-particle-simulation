#pragma once
#include "Collision/CollisionAlgorithm.h"

#include <cmath>
#include <iostream>
#include <cstdint>
#include <vector>


class SpatialSubdivisionAlgorithm : public CollisionAlgorithm
{
    public:
        int run(const Coordinates& Coordinates, double collision_dist);
    private:
        struct object_id
        {
            int particle_id {};
            int Home_cell {};
        };

        double calculate_distance(const Coordinate& p1, const Coordinate& p2);

        void initializeObjectandCellArray(
            const Coordinates& coordinates, 
            double cellLength, 
            double particleRadiusSq, 
            std::vector<int64_t>& cellIdArray, 
            std::vector<object_id>& objectIdArray);

        void decode_hash(const int64_t& hash, int64_t &grid_x, int64_t &grid_y, int64_t &grid_z);
        int64_t hash_coordinates(int64_t x, int64_t y, int64_t z);
};

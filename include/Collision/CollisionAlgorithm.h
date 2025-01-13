#pragma once

#include <vector>

using Coordinate = std::vector<double>; // Represents a 3D coordinate as a vector of doubles
using Coordinates = std::vector<Coordinate>; // Represents a list of 3D coordinates

class CollisionAlgorithm
{
    public:
        virtual ~CollisionAlgorithm() = default;

        virtual int run(const Coordinates& coordinates, double collision_dist) = 0;
};

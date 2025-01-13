#pragma once
#include "Collision/CollisionAlgorithm.h"

#include <cmath>

class BruteForceAlgorithm : public CollisionAlgorithm
{
    public:
        int run(const Coordinates& Coordinates, double collision_dist);
    private:
        double calculate_distance(const Coordinate& p1, const Coordinate& p2);
};

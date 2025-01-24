#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <memory> // For smart pointers
#include <cassert>

#include "Collision/CollisionAlgorithm.h"
#include "Collision/Algorithms/BruteForceAlgorithm.h"
#include "Collision/Algorithms/SpatialSubdivisionAlgorithm.h"
#include "FileReader/FileReader.h"
#include "Timer/Timer.h"

using Coordinate = std::vector<double>; // Represents a 3D coordinate as a vector of doubles
using Coordinates = std::vector<Coordinate>; // Represents a list of 3D coordinates

class CollisionFramework 
{
    public:
        enum class Algorithm
        {
            BruteForce = 0,
            SpatialSubdivision = 2,
            noAlgorithm = 1
        };

        struct TestSettings
        {
            int numberOfRuns {1};
            double collision_distance {0.05};
            std::string data_file_name {"../data/positions.xyz"};
            Algorithm algorithm {Algorithm::BruteForce};
        };
        
        CollisionFramework();

        static CollisionFramework generateCollisionFrameWork(const TestSettings& testsettings);
        double testAlgorithmForTime(); 
        int testAlgorithmForCollisions(); 

    private:
        void generateCollisionAlgorithm();

        TestSettings m_testSettings;
        Coordinates m_coordinatesData;
        std::unique_ptr<CollisionAlgorithm> m_collisionAlgorithm;
};

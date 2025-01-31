#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <memory> // For smart pointers
#include <cassert>

#include "Collision/CollisionAlgorithm.h"
#include "Collision/Algorithms/BruteForceAlgorithm.h"
#include "Collision/Algorithms/SpatialSubdivisionAlgorithm.h"
#include "Collision/Algorithms/SpatialSubdivisionOpenMPAlgorithm.h"
#include "FileReader/FileReader.h"
#include "Timer/Timer.h"

using Coordinate = std::vector<double>; // Represents a 3D coordinate as a vector of doubles
using Coordinates = std::vector<Coordinate>; // Represents a list of 3D coordinates

/**
 * @class CollisionFramework
 * Manages collision detection testing by configuring and running different collision detection algorithms.
 */
class CollisionFramework 
{
    public:
        /**
         * @enum Algorithm
         * Enumeration of the supported collision detection algorithms.
         */
        enum class Algorithm
        {
            BruteForce = 0,                 ///< Brute-force collision detection.
            SpatialSubdivision = 2,         ///< Spatial subdivision-based collision detection.
            SpatialSubdivisionParallel = 3, ///< Parallel spatial subdivision-based collision detection.
            noAlgorithm = 1                 ///< Placeholder for no algorithm selected.
        };

        /**
         * @struct TestSettings
         * Contains configuration settings for testing collision detection algorithms.
         */
        struct TestSettings
        {
            int numberOfRuns {1};                   ///< Number of test runs to average performance.
            double collision_distance {0.05};      ///< Collision detection threshold distance.
            std::string data_file_name {"../data/positions.xyz"}; ///< Path to the input data file.
            Algorithm algorithm {Algorithm::BruteForce}; ///< Selected collision detection algorithm.
        };
        
        // Default constructor that initializes an empty framework with default settings.
        CollisionFramework();

        /**
         * Factory method to generate a configured CollisionFramework instance.
         *
         * @param testsettings The settings for the collision framework.
         * @return A configured CollisionFramework instance.
         * @throws std::runtime_error If the collision algorithm settings are invalid.
         */
        static CollisionFramework generateCollisionFrameWork(const TestSettings& testsettings);

        /**
         * Tests the selected collision detection algorithm for runtime performance.
         *
         * @return The average runtime in seconds over the number of test runs.
         */
        double testAlgorithmForTime();

        /**
         * Tests the selected collision detection algorithm for the number of detected collisions.
         *
         * @return The number of detected collisions.
         */
        int testAlgorithmForCollisions();

    private:
        /**
         * Generates the collision detection algorithm based on the selected test settings.
         * Sets the appropriate algorithm in `m_collisionAlgorithm`.
         */
        void generateCollisionAlgorithm();

        TestSettings m_testSettings; ///< Stores the test settings for the framework.
        Coordinates m_coordinatesData; ///< Holds the 3D particle data loaded from the file.
        std::unique_ptr<CollisionAlgorithm> m_collisionAlgorithm; ///< Pointer to the selected collision detection algorithm.
};

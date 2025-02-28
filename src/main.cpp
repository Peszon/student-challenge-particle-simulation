#include "../include/CollisionFramework/CollisionFramework.h"

int main()
{
    CollisionFramework::TestSettings testSettings = {
        10,                                                       // number of runs the test algorithm for time will average over.
        0.05,                                                     // Distance at witch collision occurs.
        "../student-challenge-particle-simulation/data/positions_large.xyz",                            // Name of data file.
        CollisionFramework::Algorithm::SpatialSubdivisionParallel // Choice of algorithm.
    };

    CollisionFramework framework = CollisionFramework::generateCollisionFrameWork(testSettings);

    // Testing the algorithm.
    double timeTaken = framework.testAlgorithmForTime();
    int collisions = framework.testAlgorithmForCollisions();

    // Printing the results.
    std::cout << "\nTime taken:          " << timeTaken << " seconds\n";
    std::cout << "Collisions detected: " << collisions << '\n';

    return 0;
}

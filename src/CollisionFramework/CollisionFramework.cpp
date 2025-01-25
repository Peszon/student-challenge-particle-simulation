#include "CollisionFramework/CollisionFramework.h"

CollisionFramework::CollisionFramework() 
    : m_testSettings {},
      m_coordinatesData {},
      m_collisionAlgorithm {}
{    
}


CollisionFramework CollisionFramework::generateCollisionFrameWork(const TestSettings& testsettings) 
{
    CollisionFramework collisionFramework {};

    collisionFramework.m_testSettings = testsettings;

    collisionFramework.m_coordinatesData = DataReader::read_positions_data(collisionFramework.m_testSettings.data_file_name);
    collisionFramework.generateCollisionAlgorithm();

    if (!collisionFramework.m_collisionAlgorithm) 
        throw std::runtime_error("Error: Invalid collision framework algorithm settings.");

    std::cout << "Collision framework successfully created." << std::endl;
    return collisionFramework;
}


double CollisionFramework::testAlgorithmForTime() 
{
    Timer t;
    for(int i = 0; i < m_testSettings.numberOfRuns; i++) {
        m_collisionAlgorithm -> run(m_coordinatesData, m_testSettings.collision_distance);
    }
    return t.elapsed() / m_testSettings.numberOfRuns;
}


int CollisionFramework::testAlgorithmForCollisions() 
{
    return (m_collisionAlgorithm -> run(m_coordinatesData, m_testSettings.collision_distance));
}


void CollisionFramework::generateCollisionAlgorithm()
{
    switch (m_testSettings.algorithm) 
    {
        case Algorithm::BruteForce: 
            m_collisionAlgorithm = std::make_unique<BruteForceAlgorithm>();
            break;
        case Algorithm::SpatialSubdivision: 
            m_collisionAlgorithm = std::make_unique<SpatialSubdivisionAlgorithm>();
            break;
        case Algorithm::noAlgorithm:
            m_collisionAlgorithm = nullptr;
            break;
        default:
            m_collisionAlgorithm = nullptr;
    }
}

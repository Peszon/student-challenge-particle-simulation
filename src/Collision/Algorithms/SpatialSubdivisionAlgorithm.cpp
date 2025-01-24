#include "Collision/Algorithms/SpatialSubdivisionAlgorithm.h"

double SpatialSubdivisionAlgorithm::calculate_distance(const Coordinate& pos1, const Coordinate& pos2) {
    return sqrt(pow(pos1[0] - pos2[0], 2) + pow(pos1[1] - pos2[1], 2) + pow(pos1[2] - pos2[2], 2));
}

int SpatialSubdivisionAlgorithm::run(const Coordinates& coordinates, double collision_distance) 
{
    constexpr double cellLength {0.05 * 1.5}; // The times 1.5 makes it so that the maximum phantom cells is 8.
    constexpr double particleRadiusSq { 0.025*0.025 };

    std::vector<int64_t> cellIdArray {};
    std::vector<object_id> objectIdArray {};

    initializeObjectandCellArray(coordinates, cellLength, particleRadiusSq, cellIdArray, objectIdArray);

    for (int i = 0; i < 10; i++) {
        int64_t grid_x, grid_y, grid_z;
        decode_hash(cellIdArray[i], grid_x, grid_y, grid_z);
        std::cout << "Particle id: " << objectIdArray[i].particle_id << ", Home cell: " << objectIdArray[i].Home_cell << "\n";
        std::cout << "Discretized x: " << grid_x << ", y: " << grid_y << ", z: " << grid_z << "\n" << std::endl;
    }

    return 0;
}

void SpatialSubdivisionAlgorithm::initializeObjectandCellArray(const Coordinates& coordinates, double cellLength, double particleRadiusSq, std::vector<int64_t>& cellIdArray, std::vector<object_id>& objectIdArray) {
    int object_id_counter = 1;
    for (const Coordinate &coordinate : coordinates) {
        // Finds the cell grid coordinate of the home cell.
        int64_t grid_x = static_cast<int64_t>(std::floor(coordinate[0] / cellLength));
        int64_t grid_y = static_cast<int64_t>(std::floor(coordinate[1] / cellLength));
        int64_t grid_z = static_cast<int64_t>(std::floor(coordinate[2] / cellLength));

        // Hashes and adds the home cell to cellIdArray list.
        cellIdArray.push_back(hash_coordinates(grid_x, grid_y, grid_z));

        // Adds corresponding information regarding the particle to the ObjectIdArray.
        objectIdArray.push_back(object_id {object_id_counter, 1});

        // Checks if the particle is present in any neighbourign cells.
        for (int64_t dx = -1; dx <= 1; ++dx) {
            for (int64_t dy = -1; dy <= 1; ++dy) {
                for (int64_t dz = -1; dz <= 1; ++dz) {
                    if (dx == 0 && dy == 0 && dz == 0) continue;

                    double minX = (grid_x + dx) * cellLength;
                    double maxX = minX + cellLength;
                    double minY = (grid_y + dy) * cellLength;
                    double maxY = minY + cellLength;
                    double minZ = (grid_z + dz) * cellLength;
                    double maxZ = minZ + cellLength;

                    // Finding the nearest point of the relevant cell.
                    double nearestX = std::max(minX, std::min(coordinate[0], maxX));
                    double nearestY = std::max(minY, std::min(coordinate[1], maxY));
                    double nearestZ = std::max(minZ, std::min(coordinate[2], maxZ));

                    // Calculates distance to the relevant cell
                    double distSq = (coordinate[0] - nearestX) * (coordinate[0] - nearestX)
                                  + (coordinate[1] - nearestY) * (coordinate[1] - nearestY)
                                  + (coordinate[2] - nearestZ) * (coordinate[2] - nearestZ);
  
                    // Adds the cell to the cellIdArray array if the distance between the particle and the relevant cell is less than the particle radius.
                    if (distSq <= particleRadiusSq) {
                        cellIdArray.push_back(hash_coordinates(grid_x + dx, grid_y + dy, grid_z + dz));
                        objectIdArray.push_back(object_id {object_id_counter, 0});
                    }
                }
            }
        }

        object_id_counter++;
    }
}

int64_t SpatialSubdivisionAlgorithm::hash_coordinates(int64_t grid_x, int64_t grid_y, int64_t grid_z) {
    // Encode into a 64-bit hash
    int64_t hash = 0;
    hash |= (grid_x & 0xFFFF);         // Place X in bits 0–15
    hash |= (grid_y & 0xFFFF) << 16;   // Place Y in bits 16–31
    hash |= (grid_z & 0xFFFF) << 32;   // Place Z in bits 32–47
    return hash;
}

void SpatialSubdivisionAlgorithm::decode_hash(const int64_t& hash, int64_t& grid_x, int64_t& grid_y, int64_t& grid_z) {
    // Extract and reinterpret as signed 16-bit integers
    grid_x = static_cast<int16_t>(hash & 0xFFFF);         // Extract X from bits 0–15
    grid_y = static_cast<int16_t>((hash >> 16) & 0xFFFF); // Extract Y from bits 16–31
    grid_z = static_cast<int16_t>((hash >> 32) & 0xFFFF); // Extract Z from bits 32–47
}

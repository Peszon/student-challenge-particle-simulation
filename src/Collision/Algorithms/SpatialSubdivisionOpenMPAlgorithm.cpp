#include "Collision/Algorithms/SpatialSubdivisionOpenMPAlgorithm.h"

int SpatialSubdivisionOpenMPAlgorithm::run(
    const Coordinates& coordinates, 
    double collisionDistance) 
{
    double const cellLength {collisionDistance * 1.5}; // The times 1.5 makes it so that the maximum phantom cells is 8.
    double const collisionDistSq { collisionDistance * collisionDistance };

    std::vector<int64_t> cellIdArray {};
    std::vector<object_id> objectIdArray {};

    initializeObjectandCellArray(coordinates, cellLength, collisionDistSq, cellIdArray, objectIdArray);

    radixSort(cellIdArray, objectIdArray);

    return calculateNumberOfCollisions(coordinates, collisionDistSq, cellIdArray, objectIdArray);
}


int SpatialSubdivisionOpenMPAlgorithm::calculateNumberOfCollisions(
    const Coordinates& coordinates, 
    double collisionDistSq, 
    std::vector<int64_t>& cellIdArray,  
    std::vector<object_id>& objectIdArray) 
{
    const size_t arrayLength = cellIdArray.size();
    int collisionCounter = 0;

    #pragma omp parallel for reduction(+:collisionCounter)
    for (size_t i = 0; i < arrayLength - 1; i++) {
        for (size_t j = i + 1; j < arrayLength; j++) {
            if (cellIdArray[i] != cellIdArray[j]) {
                break;
            } else if (objectIdArray[j].Home_cell == 0) {
                continue;
            } else {
                Coordinate coord1 = coordinates[objectIdArray[i].particle_id];
                Coordinate coord2 = coordinates[objectIdArray[j].particle_id];
                if (calculateDistanceSq(coord1, coord2) < collisionDistSq) {
                    collisionCounter++;
                }
            }
        }
    }

    return collisionCounter;
}

double SpatialSubdivisionOpenMPAlgorithm::calculateDistanceSq(
    const Coordinate& pos1, 
    const Coordinate& pos2) 
{
    return pow(pos1[0] - pos2[0], 2) + pow(pos1[1] - pos2[1], 2) + pow(pos1[2] - pos2[2], 2);
}


void SpatialSubdivisionOpenMPAlgorithm::bitParallelRadix(
    const std::vector<int64_t>& cellIdArray, 
    const std::vector<object_id>& objectIdArray,
    std::vector<int64_t>& cellIdArrayTemp, 
    std::vector<object_id>& objectIdArrayTemp,
    int arrayLength,
    int numbKeys,            // Number of buckets
    int numbThreads,
    int bitshift)            // Number of threads (parts)
{
    // Using the fact that int division drops the decimals.
    int chunk = (arrayLength + numbThreads - 1) / numbThreads; 

    std::vector<std::vector<int>> Count(numbThreads, std::vector<int>(numbKeys, 0));

    // First parallel phase: local bucket counts
    #pragma omp parallel for num_threads(numbThreads)
    for(int part = 0; part < numbThreads; part++) {
        int start = part * chunk;
        int end   = std::min(start + chunk, arrayLength);
        for(int i = start; i < end; i++) {
            int bucket = (cellIdArray[i] >> bitshift) & 0xFF; 
            Count[part][bucket]++;
        }
    }

    // Prefix sums across parts to compute final positions
    int base = 0;
    for(unsigned int bucket = 0; bucket < numbKeys; bucket++) {
        for(unsigned int part = 0; part < numbThreads; part++) {
            Count[part][bucket] += base;
            base = Count[part][bucket];
        }
    }

    // Second parallel phase: place elements into the output array
    #pragma omp parallel for num_threads(numbThreads)
    for (int part = 0; part < numbThreads; part++) {
        int start = part * chunk;
        int end   = std::min(start + chunk, arrayLength);
        for (int i = end - 1; i >= start; i--) {
            unsigned int bucket = static_cast<unsigned int>((cellIdArray[i] >> bitshift) & 0xFF);
            int pos = Count[part][bucket];
            cellIdArrayTemp[pos - 1] = cellIdArray[i];
            objectIdArrayTemp[pos - 1] = objectIdArray[i];
            Count[part][bucket]--;
        }
    }
}

void SpatialSubdivisionOpenMPAlgorithm::radixSort(
    std::vector<int64_t>& cellIdArray,  
    std::vector<object_id>& objectIdArray)
{
    const size_t arrayLength = cellIdArray.size();
    const int bitshift_step = 8;
    const int numbKeys = 1 << bitshift_step;;
    const int numbThreads = omp_get_max_threads();
    
    // Temporary arrays, resized to match the orignial arrays.
    std::vector<int64_t> cellIdArrayTemp(arrayLength); 
    std::vector<object_id> objectIdArrayTemp(arrayLength);

    for (int bitshift = 0; bitshift <= 40; bitshift += bitshift_step) {
        bitParallelRadix(cellIdArray, objectIdArray, cellIdArrayTemp, objectIdArrayTemp, arrayLength, numbKeys, numbThreads, bitshift);

        cellIdArray.swap(cellIdArrayTemp);
        objectIdArray.swap(objectIdArrayTemp);
    }
}


void SpatialSubdivisionOpenMPAlgorithm::initializeObjectandCellArray(
    const Coordinates& coordinates, 
    double cellLength, 
    double collisionDistSq, 
    std::vector<int64_t>& cellIdArray, 
    std::vector<object_id>& objectIdArray) 
{
    // We will accumulate in local thread buffers to avoid contention.
    // Note: You may need to estimate max size or dynamically use a local std::vector.

    // Make an array of vectors, one per thread:
    #pragma omp parallel
    {
        std::vector<int64_t> cellIdsLocal;
        std::vector<object_id> objIdsLocal;

        #pragma omp for
        for (int object_id_counter = 0; object_id_counter < (int)coordinates.size(); ++object_id_counter) {
            const Coordinate& coordinate = coordinates[object_id_counter];
            
            // 1) Compute home cell:
            int64_t grid_x = static_cast<int64_t>(std::floor(coordinate[0] / cellLength));
            int64_t grid_y = static_cast<int64_t>(std::floor(coordinate[1] / cellLength));
            int64_t grid_z = static_cast<int64_t>(std::floor(coordinate[2] / cellLength));
            
            // 2) Push into local buffers
            cellIdsLocal.push_back(hash_coordinates(grid_x, grid_y, grid_z));
            objIdsLocal.push_back(object_id{object_id_counter, 1});
            
            // 3) Check neighbors
            for (int64_t dx = -1; dx <= 1; ++dx) {
                for (int64_t dy = -1; dy <= 1; ++dy) {
                    for (int64_t dz = -1; dz <= 1; ++dz) {
                        if (dx == 0 && dy == 0 && dz == 0) {
                            continue;
                        }

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

                        if (distSq <= collisionDistSq) {
                            cellIdsLocal.push_back(
                                hash_coordinates(grid_x + dx, grid_y + dy, grid_z + dz));
                            objIdsLocal.push_back(object_id{object_id_counter, 0});
                        }
                    }
                }
            }
        } 

        // Merge each thread's local vectors into the global ones:
        #pragma omp critical
        {
            cellIdArray.insert(cellIdArray.end(), cellIdsLocal.begin(), cellIdsLocal.end());
            objectIdArray.insert(objectIdArray.end(), objIdsLocal.begin(), objIdsLocal.end());
        }
    } // end parallel
}


int64_t SpatialSubdivisionOpenMPAlgorithm::hash_coordinates(
    int64_t grid_x, 
    int64_t grid_y, 
    int64_t grid_z) 
{
    // Encode into a 64-bit hash
    int64_t hash = 0;
    hash |= (grid_x & 0xFFFF);         // Place X in bits 0–15
    hash |= (grid_y & 0xFFFF) << 16;   // Place Y in bits 16–31
    hash |= (grid_z & 0xFFFF) << 32;   // Place Z in bits 32–47
    return hash;
}


void SpatialSubdivisionOpenMPAlgorithm::decode_hash(
    const int64_t& hash, 
    int64_t& grid_x, 
    int64_t& grid_y, 
    int64_t& grid_z) {
    // Extract and reinterpret as signed 16-bit integers
    grid_x = static_cast<int16_t>(hash & 0xFFFF);         // Extract X from bits 0–15
    grid_y = static_cast<int16_t>((hash >> 16) & 0xFFFF); // Extract Y from bits 16–31
    grid_z = static_cast<int16_t>((hash >> 32) & 0xFFFF); // Extract Z from bits 32–47
}

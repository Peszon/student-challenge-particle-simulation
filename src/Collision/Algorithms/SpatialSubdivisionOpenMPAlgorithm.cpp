#include "Collision/Algorithms/SpatialSubdivisionOpenMPAlgorithm.h"

int SpatialSubdivisionOpenMPAlgorithm::run(
    const Coordinates& coordinates, 
    double collisionDistance) 
{
    // The times 1.5 makes it so that the maximum number of phantom cells is 8.
    double const cellLength {collisionDistance * 1.5}; 

    // Calculates the collision distance square once to remove square root operations.
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

    // Each thread takes one particle and compares it to other particles in the same cell.
    #pragma omp parallel for reduction(+:collisionCounter)
    for (size_t i = 0; i < arrayLength - 1; i++) {
        for (size_t j = i + 1; j < arrayLength; j++) {

            // Move the first index if the second index is in another cell.
            if (cellIdArray[i] != cellIdArray[j]) {
                break;
            
            // Preventing double-counting of collisions:
            /* If two particles collide within their shared home cell, the stability of the Radix sort 
            * ensures that the collision is counted only once. If they belong to different home cells, 
            * each particle can be considered a sphere, and counting the collision only when one particle 
            * (particle at element i) is in its home cell is sufficient to avoid double-counting.
            */
            } else if (objectIdArray[j].Home_cell == 0) {
                continue;

            // Do a distance calculation to see if they are close enough to collide.
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
    int numbKeys,            // Number of digits
    int numbThreads,         // Number of threads (parts)   
    int bitshift)            // Which bits the current sort looks at, e.g. 0-7, 8-15, ..., 40-47.
{
    // Using the fact that int division drops the decimals.
    int chunk = (arrayLength + numbThreads - 1) / numbThreads; 

    // Creates a radix counter for each thread.
    std::vector<std::vector<int>> Count(numbThreads, std::vector<int>(numbKeys, 0));

    // First parallel phase: counts the occurance of each digit accros all the entries in cellIdArray.
    #pragma omp parallel for num_threads(numbThreads)
    for(int part = 0; part < numbThreads; part++) {
        int start = part * chunk;
        int end   = std::min(start + chunk, arrayLength);
        for(int i = start; i < end; i++) {
            int digit = (cellIdArray[i] >> bitshift) & 0xFF; 
            Count[part][digit]++;
        }
    }


    // Prefix sequential sums across parts to compute final positions.
    // This prefix sum could be parallelized, if I had more time/theory, using Hillis and Steele's 
    // algorithm or Blelloch's algorithm.  
    int base = 0;
    for(unsigned int digit = 0; digit < numbKeys; digit++) {
        for(unsigned int part = 0; part < numbThreads; part++) {
            Count[part][digit] += base;
            base = Count[part][digit];
        }
    }

    // Second parallel phase: placing elements into the output array
    #pragma omp parallel for num_threads(numbThreads)
    for (int part = 0; part < numbThreads; part++) {
        int start = part * chunk;
        int end   = std::min(start + chunk, arrayLength);

        // for-loop iterating backwards to ensure the stability of the sort.
        for (int i = end - 1; i >= start; i--) {
            unsigned int digit = static_cast<unsigned int>((cellIdArray[i] >> bitshift) & 0xFF);
            int pos = Count[part][digit];
            cellIdArrayTemp[pos - 1] = cellIdArray[i];
            objectIdArrayTemp[pos - 1] = objectIdArray[i];
            Count[part][digit]--;
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
    
    // Temporary arrays.
    std::vector<int64_t> cellIdArrayTemp(arrayLength); 
    std::vector<object_id> objectIdArrayTemp(arrayLength);

    // Doing the bitParallelRadix in 8 bit jump from 0-7 to 40-47 since the 
    // CellIdArray consists of three 16 bit hashes of each coordiante.
    for (int bitshift = 0; bitshift <= 40; bitshift += bitshift_step) {
        bitParallelRadix(cellIdArray, objectIdArray, cellIdArrayTemp, objectIdArrayTemp, arrayLength, numbKeys, numbThreads, bitshift);

        // Swaping so that the original array is the sorted one.
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
    // We will accumulate in local thread buffer vectors to avoid contention.
    #pragma omp parallel
    {
        // Make an array of vectors, one per thread:
        std::vector<int64_t> cellIdsLocal;
        std::vector<object_id> objIdsLocal;

        #pragma omp for
        for (int object_id_counter = 0; object_id_counter < (int)coordinates.size(); ++object_id_counter) {
            const Coordinate& coordinate = coordinates[object_id_counter];
            
            // First compute home cell:
            int64_t grid_x = static_cast<int64_t>(std::floor(coordinate[0] / cellLength));
            int64_t grid_y = static_cast<int64_t>(std::floor(coordinate[1] / cellLength));
            int64_t grid_z = static_cast<int64_t>(std::floor(coordinate[2] / cellLength));
            
            // Push the home cell into local buffers
            cellIdsLocal.push_back(hash_coordinates(grid_x, grid_y, grid_z));
            objIdsLocal.push_back(object_id{object_id_counter, 1});
            
            // Check if particle resides inside any neighbouring cells.
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

                        // Add the neighbouring cell if it is closer than the collision distance.
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

#include "Collision/Algorithms/SpatialSubdivisionAlgorithm.h"


int SpatialSubdivisionAlgorithm::run(const Coordinates& coordinates, double collisionDistance) 
{
    double const cellLength {collisionDistance * 1.5}; // The times 1.5 makes it so that the maximum phantom cells is 8.
    double const collisionDistSq { collisionDistance * collisionDistance };

    std::vector<int64_t> cellIdArray {};
    std::vector<object_id> objectIdArray {};

    initializeObjectandCellArray(coordinates, cellLength, collisionDistSq, cellIdArray, objectIdArray);

    radixSort(cellIdArray, objectIdArray);

    return calculateNumberOfCollisions(coordinates, collisionDistSq, cellIdArray, objectIdArray);
}


int SpatialSubdivisionAlgorithm::calculateNumberOfCollisions(const Coordinates& coordinates, double collisionDistSq, std::vector<int64_t>& cellIdArray,  std::vector<object_id>& objectIdArray) 
{
    const size_t arrayLength = cellIdArray.size();
    int collisionCounter { 0 };

    size_t i = 0;
    // Loop over i as long as there is at least one more element to compare with i.
    while (i < (arrayLength - 1)) {
        
        // Inner loop used to iterate over elements to compare to the element at i. 
        for(size_t j = i + 1; j < arrayLength; j++) {
            
            // If the cell IDs differ, we break from the inner loop because all further elements will also be in a different cell.
            if (cellIdArray[i] != cellIdArray[j]) {
                break;

            // Only checking if there is a collision if in the i:ths particle's home cell.
            // Thereby not 'double counting' the collision that would have occured in the j:th particle's home.
            } else if (objectIdArray[j].Home_cell == 0) {      
                continue; 

            // At this point, both particles share the same cell, and we ignore the 'double counting' cases.
            // If the distance between them is less than the collision distance we count it as a collision.
            } else {
                Coordinate coordParticle1 = coordinates[objectIdArray[i].particle_id];
                Coordinate coordParticle2 = coordinates[objectIdArray[j].particle_id];
                if (calculateDistanceSq(coordParticle1, coordParticle2) < collisionDistSq) {
                    collisionCounter++;
                }
            }
        }

        i++;
    }

    return collisionCounter;
}


double SpatialSubdivisionAlgorithm::calculateDistanceSq(const Coordinate& pos1, const Coordinate& pos2) {
    return pow(pos1[0] - pos2[0], 2) + pow(pos1[1] - pos2[1], 2) + pow(pos1[2] - pos2[2], 2);
}


void SpatialSubdivisionAlgorithm::radixSort(std::vector<int64_t>& cellIdArray,  std::vector<object_id>& objectIdArray) {
    const size_t arrayLength = cellIdArray.size();
    
    // Temporary arrays, resized to match the orignial arrays.
    std::vector<int64_t> replCellIdArray(arrayLength); 
    std::vector<object_id> replObjectIdArray(arrayLength);

    // Since we use three coordinates each represented by 16 bits, we only need to sort for the first 48 bits. 
    for (int bitShift = 0; bitShift <= 40; bitShift += 8) {
        std::vector<int> radixCounter(256,0);

        // Counting the number of occurances of each number in a 8-bit.
        for (const int64_t& cellID : cellIdArray) {
            radixCounter[((cellID >> bitShift) & 0xFF)]++;
        } 

        // Calculating the Radix prefix sum. 
        for (unsigned int i = 1; i < 256; i++) {
            radixCounter[i] += radixCounter[i - 1];
        }

        // Adds the element of the original array to the new array at the place specified by the Radix prefex sum.
        for (int i = static_cast<int>(arrayLength) - 1; i >= 0; i--) {
            unsigned int digit = static_cast<unsigned int>((cellIdArray[i] >> bitShift) & 0xFF);

            replCellIdArray[radixCounter[digit] - 1] = cellIdArray[i];
            replObjectIdArray[radixCounter[digit] - 1] = objectIdArray[i];

            radixCounter[digit]--;
        }

        // Switches the arrays so that the original array is sorted.
        cellIdArray = replCellIdArray;
        objectIdArray = replObjectIdArray;
    }
} 


void SpatialSubdivisionAlgorithm::initializeObjectandCellArray(const Coordinates& coordinates, double cellLength, double collisionDistSq, std::vector<int64_t>& cellIdArray, std::vector<object_id>& objectIdArray) {
    int object_id_counter = 0;
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
                    if (distSq <= collisionDistSq) {
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

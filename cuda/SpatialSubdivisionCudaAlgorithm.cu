#include "Collision/CollisionAlgorithm.h"

#include <cmath>
#include <iostream>
#include <cstdint>
#include <vector>
#include <bitset>

#define MAX_BLOCK_SIZE 384 //Dividable by 32 so each block work with full wraps. 

// CUDA error checking macro.
#define CUDA_CHECK(error)                                                    \
    do {                                                                     \
        cudaError_t err = error;                                             \
        if (err != cudaSuccess) {                                            \
            std::cerr << "CUDA Error: " << cudaGetErrorString(err)           \
                      << " at " << __FILE__ << ":" << __LINE__ << "\n";      \
            exit(EXIT_FAILURE);                                              \
        }                                                                    \
    } while (0)


struct object_id
{
    int particle_id {0}; ///< The unique ID of the particle.
    int home_cell {0};   ///< Flag indicating if the cell is the particle's home cell.
};


double* vecToArr(const std::vector<std::vector<double>>& vec) {
    std::size_t totalsize = 0;

    for (int i=0; i<vec.size(); i++) {
        totalsize += vec[i].size();
    }

    double* newarr=new double[totalsize];
    double* walkarr=newarr;

    for (int i=0; i<vec.size(); i++) {
        std::copy(vec[i].begin(), vec[i].end(), walkarr);
        walkarr += vec[i].size();
    }

    return newarr;
}

__device__ 
int64_t hash_coordinates(int64_t grid_x, int64_t grid_y, int64_t grid_z) {
    int64_t hash = 0;
    hash |= (grid_x & 0xFFFF);         // Place X in bits 0–15
    hash |= (grid_y & 0xFFFF) << 16;   // Place Y in bits 16–31
    hash |= (grid_z & 0xFFFF) << 32;   // Place Z in bits 32–47
    return hash;
}

__global__ 
void initializeObjectandCellArray(
    double *coordinates, 
    int64_t *cellIdArray, 
    object_id *objectIdArray, 
    int coordinatesLength, 
    double collisionDist) 
{
    double const cellLength {collisionDist * 1.5}; 
    double const collisionDistSq {collisionDist * collisionDist};

    int TIdx {threadIdx.x + blockDim.x * blockIdx.x};
    int coordIdx {3 * TIdx};
    int arrayIdx {8 * TIdx};

    // Finds the cell grid coordinate of the home cell.
    if  (TIdx < coordinatesLength) {
        int64_t grid_x = static_cast<int64_t>(floor(coordinates[coordIdx]     / cellLength));
        int64_t grid_y = static_cast<int64_t>(floor(coordinates[coordIdx + 1] / cellLength));
        int64_t grid_z = static_cast<int64_t>(floor(coordinates[coordIdx + 2] / cellLength));

        // Hashes and adds the home cell to cellIdArray list.
        cellIdArray[arrayIdx] = (hash_coordinates(grid_x, grid_y, grid_z));

        // Adds corresponding information regarding the particle to the ObjectIdArray.
        objectIdArray[arrayIdx] = object_id{TIdx, 1};

        // Checks if the particle is present in any neighbourign cells. Can be in maxiumum 7 more cells due to geometry
        int counter {1};
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
                    double nearestX = fmax(minX, fmin(coordinates[coordIdx], maxX));
                    double nearestY = fmax(minY, fmin(coordinates[coordIdx + 1], maxY));
                    double nearestZ = fmax(minZ, fmin(coordinates[coordIdx + 2], maxZ));

                    // Calculates distance to the relevant cell
                    double distSq =   (coordinates[coordIdx]     - nearestX) *     (coordinates[coordIdx] - nearestX)
                                    + (coordinates[coordIdx + 1] - nearestY) * (coordinates[coordIdx + 1] - nearestY)
                                    + (coordinates[coordIdx + 2] - nearestZ) * (coordinates[coordIdx + 2] - nearestZ);

                    // Adds the cell to the cellIdArray array if the distance between the particle and the relevant cell is less than the particle radius.
                    if (distSq <= collisionDistSq) {
                        cellIdArray[arrayIdx + counter] = (hash_coordinates(grid_x + dx, grid_y + dy, grid_z + dz));
                        objectIdArray[arrayIdx + counter] = object_id{TIdx, 0};
                        counter++;
                    }
                }
            }
        }
        for(; counter < 8; counter++) {
            cellIdArray[arrayIdx + counter] = 0xFFFFFFFFFFFF;  //Setting to a default value the particle doesn't cross into 7 other cells. 
            objectIdArray[arrayIdx + counter] = object_id{-1, 0};
        }
    }
}


__global__ 
void placeCellIdFlags(
    int64_t *d_cellIdArray, 
    int64_t *d_cellIdArrayTemp,
    object_id *d_objectIdArray, 
    object_id *d_objectIdArrayTemp, 
    int *d_flagCellIDStart,
    int *d_flagCellIDEnd,
    unsigned int arrayIdLength) 
{
    int TIdx = threadIdx.x + blockDim.x * blockIdx.x;

    // Shift the original d_cellIdArray one step Forwards.
    if (TIdx < (arrayIdLength - 1)) {
        d_cellIdArrayTemp[TIdx + 1] = d_cellIdArray[TIdx];
    }
    
    __syncthreads();

    // Check if the element on index i is not the same as the element on index i-1. If true place a start marker.
    if (TIdx < (arrayIdLength) && TIdx != 0) {
        d_flagCellIDStart[TIdx] = ((d_cellIdArrayTemp[TIdx] != d_cellIdArray[TIdx]) ? 1 : 0);        
    } else if  (TIdx == 0) {
        d_flagCellIDStart[TIdx] = 1;
    }

    __syncthreads();

    // Create a end marker array by shifting the start marker array on step to the left.
    if (TIdx < arrayIdLength && TIdx != 0) {
        d_flagCellIDEnd[TIdx - 1] = d_flagCellIDStart[TIdx];
    } else if (TIdx == arrayIdLength) {
        d_flagCellIDStart[TIdx] = 1;
    }
}


__device__ 
double calculateDistanceSq(
    double x1, 
    double y1, 
    double z1,
    double x2,
    double y2, 
    double z2) 
{
    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;
    return dx * dx + dy * dy + dz * dz;
}

__device__ 
void checkCollision(
    int i, 
    int j,
    const double *coordinates,
    const object_id *objectIdArray,
    double collDistSq,
    int *counter) 
{
    // Only do distance checks if home_id == 1
    if (objectIdArray[i].home_cell == 1) {
        double p1 = objectIdArray[i].particle_id; 
        double p2 = objectIdArray[j].particle_id; 

        double distSq = calculateDistanceSq(
            coordinates[3 * p1], 
            coordinates[3 * p1 + 1], 
            coordinates[3 * p1 + 2], 
            coordinates[3 * p2], 
            coordinates[3 * p2 + 1], 
            coordinates[3 * p2 + 2]);

        if (distSq < collDistSq) {
            counter[i]++; 
        }
    }
}

__global__ 
void calculateNumberOfCollisions(
    const int64_t *d_cellIdArray,
    const object_id *d_objectIdArray,  
    const int *d_flagCellIdStart,
    const int *d_flagCellIdEnd,
    const double *coordinates,
    double collDistSq,
    unsigned int arrayIdLength,
    int *d_collisionCounter)
{
    // Thread index
    int TIdx = blockIdx.x * blockDim.x + threadIdx.x;

    // Bounds check
    if (TIdx >= arrayIdLength) return;

    // Each thread runs only if it sees a start flag of 1 (example condition)
    if (d_flagCellIdStart[TIdx] == 1 && d_flagCellIdEnd[TIdx] != 1) {

        // Outer loop over i
        for (int i = TIdx; i < arrayIdLength; i++) {
            // Check collisions with subsequent j
            for (int j = i + 1; j < arrayIdLength; j++) {
                // Perform collision check
                checkCollision(
                    i,
                    j, 
                    coordinates, 
                    d_objectIdArray, 
                    collDistSq, 
                    d_collisionCounter);

                // If flag says this is the last iteration for j, break afterward
                if (d_flagCellIdEnd[j] == 1) {
                    // One final iteration has just occurred, so exit inner loop
                    break;
                }
            }

            // If flag says this is the last iteration for i, break afterward
            if (d_flagCellIdEnd[i + 1] == 1) {
                break;
            }
        }
    }
}


int run(Coordinates h_vectorCoordinates, double collisionDistance) {
    unsigned int coordinatesLength = h_vectorCoordinates.size();
    unsigned int coordinatesSize = coordinatesLength * 3 * sizeof(double);  // since h_vectorCoordinates contained 3 double per row.
    unsigned int arrayIdLength = coordinatesLength * 8; // For every 3D coordinate get 1 cellId.
    unsigned int arrayIdSize = coordinatesLength* 8 * sizeof(int64_t);

    // Move the coordinate to device global memory.
    double *h_coordinates = vecToArr(h_vectorCoordinates);
    double *d_coordinates; CUDA_CHECK(cudaMalloc(&d_coordinates, coordinatesSize));

    CUDA_CHECK(cudaMemcpy(d_coordinates, h_coordinates, coordinatesSize, cudaMemcpyHostToDevice));

    // Create cellIdArray and objectIdArray.
    int64_t *d_cellIdArray; CUDA_CHECK(cudaMalloc(&d_cellIdArray, arrayIdSize));
    object_id *d_objectIdArray; CUDA_CHECK(cudaMalloc(&d_objectIdArray, arrayIdLength * sizeof(object_id)));

    // Create array on device for out data from the Radix sort.        
    int64_t *d_cellIdArrayOut; CUDA_CHECK(cudaMalloc(&d_cellIdArrayOut, arrayIdSize));
    object_id *d_objectIdArrayOut; CUDA_CHECK(cudaMalloc(&d_objectIdArrayOut, arrayIdLength * sizeof(object_id)));

    // Take advantage of the fact that integer division drops the decimals
    unsigned int gridSize = coordinatesLength / MAX_BLOCK_SIZE;
    if (coordinatesLength % MAX_BLOCK_SIZE != 0) {
        gridSize += 1;
    }

    // Create the cellIdArray and the ObjectIdArray for the obtained coordinates.
    initializeObjectandCellArray<<<gridSize, MAX_BLOCK_SIZE>>>(
        d_coordinates, 
        d_cellIdArray, 
        d_objectIdArray, 
        coordinatesLength, 
        collisionDist);
    
    cudaDeviceSynchronize();

    // Create a variable to hold the temporary storage.
    size_t temp_storage_bytes = 0;
    void *d_temp_storage = nullptr;

    // First call sets temp_storage_bytes
    cub::DeviceRadixSort::SortPairs(nullptr, temp_storage_bytes,
                                    d_cellIdArray, d_cellIdArrayOut,
                                    d_objectIdArray, d_objectIdArrayOut,
                                    arrayIdLength);

    CUDA_CHECK(cudaMalloc(&d_temp_storage, temp_storage_bytes));

    // Second call does the actual sorting
    cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
                                    d_cellIdArray, d_cellIdArrayOut,
                                    d_objectIdArray, d_objectIdArrayOut,
                                    arrayIdLength);

    cudaDeviceSynchronize();

    // Take advantage of the fact that integer division drops the decimals
    unsigned int gridSizeNmbCollisions = arrayIdLength / MAX_BLOCK_SIZE;
    if (coordinatesLength % MAX_BLOCK_SIZE != 0) {
        gridSizeNmbCollisions += 1;
    }

    // Create temporary arrays
    int64_t *d_cellIdArrayTemp = d_cellIdArray;
    object_id *d_objectIdArrayTemp = d_objectIdArray;

    // Create arrays on device that store information regarding when the cell id changes. 
    int *d_flagCellIdStart; CUDA_CHECK(cudaMalloc(&d_flagCellIdStart, arrayIdLength * sizeof(int)));
    int *d_flagCellIdEnd; CUDA_CHECK(cudaMalloc(&d_flagCellIdEnd, arrayIdLength * sizeof(int)));

    // Execute kernal that places start and end flags of the blocks of cell id in cellIdArray.
    placeCellIdFlags<<<gridSizeNmbCollisions, MAX_BLOCK_SIZE>>>(
        d_cellIdArrayOut, 
        d_cellIdArrayTemp, 
        d_objectIdArrayOut, 
        d_objectIdArrayTemp,
        d_flagCellIdStart,
        d_flagCellIdEnd,
        arrayIdLength);

    cudaDeviceSynchronize();

    // Creates an collision counter array on the device, length equals the number of threads.
    int *d_collisionCounter; CUDA_CHECK(cudaMalloc((&d_collisionCounter), arrayIdLength * sizeof(int)));
    CUDA_CHECK(cudaMemset(d_collisionCounter, 0, arrayIdLength * sizeof(int)));
    double collDistSq = collisionDistance * collisionDistance;

    calculateNumberOfCollisions<<<gridSizeNmbCollisions, MAX_BLOCK_SIZE>>>(
        d_cellIdArrayOut,
        d_objectIdArrayOut, 
        d_flagCellIdStart,
        d_flagCellIdEnd,
        d_coordinates,
        collDistSq,
        arrayIdLength,
        d_collisionCounter);
        
    cudaDeviceSynchronize();


    // Reducing the count array with built in cub function.
    int *d_collisionCounterOut = nullptr;
    CUDA_CHECK(cudaMalloc((void**)&d_collisionCounterOut, sizeof(int)));

    // Resetting these so they don't clash with the earlier sort storage
    size_t temp_storage_bytes_reduce = 0;
    void *d_temp_storage_reduce = nullptr;

    // Checking for temporary storage
    cub::DeviceReduce::Sum(
        nullptr, 
        temp_storage_bytes_reduce,
        d_collisionCounter, 
        d_collisionCounterOut, 
        arrayIdLength);

    CUDA_CHECK(cudaMalloc(&d_temp_storage_reduce, temp_storage_bytes_reduce));

    // Perform reduction
    cub::DeviceReduce::Sum(
        d_temp_storage_reduce, 
        temp_storage_bytes_reduce,
        d_collisionCounter,
        d_collisionCounterOut, 
        arrayIdLength);

    // Now copy single int to host
    int hostCollisionSum = 0;
    CUDA_CHECK(cudaMemcpy(&hostCollisionSum, d_collisionCounterOut, 
                        sizeof(int), cudaMemcpyDeviceToHost));

    // Freeing memory after operations.
    delete[] h_coordinates;

    cudaFree(d_coordinates);
    
    cudaFree(d_cellIdArray);
    cudaFree(d_objectIdArray);
    cudaFree(d_cellIdArrayOut);
    cudaFree(d_objectIdArrayOut);

    cudaFree(d_temp_storage);

    cudaFree(d_flagCellIdStart);
    cudaFree(d_flagCellIdEnd);

    cudaFree(d_collisionCounter);
    cudaFree(d_collisionCounterOut);
    cudaFree(d_temp_storage_reduce);

    return hostCollisionSum;
}

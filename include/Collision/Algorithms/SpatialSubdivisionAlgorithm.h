#pragma once
#include "Collision/CollisionAlgorithm.h"

#include <cmath>
#include <iostream>
#include <cstdint>
#include <vector>
#include <bitset>

/**
 * @class SpatialSubdivisionAlgorithm
 * Implements a collision detection algorithm based on spatial subdivision.
 * This algorithm divides the space into a grid of cells and processes potential collisions
 * by considering particles within the same or neighboring cells.
 */
class SpatialSubdivisionAlgorithm : public CollisionAlgorithm
{
    public:
        /**
         * Runs the spatial subdivision collision detection algorithm.
         *
         * @param coordinates A vector of 3D coordinates representing particle positions.
         * @param collision_dist The maximum distance to consider for detecting a collision.
         * @return The number of detected collisions.
         */
        int run(const Coordinates& coordinates, double collision_dist);

    private:
        /**
         * @struct object_id
         * Represents an object with its particle ID and home cell flag.
         * `Home_cell` is 1 for the particle's home cell and 0 for neighboring cells.
         */
        struct object_id
        {
            int particle_id {}; ///< The unique ID of the particle.
            int Home_cell {};   ///< Flag indicating if the cell is the particle's home cell.
        };

        /**
         * Initializes the cell and object arrays for the spatial grid.
         *
         * @param coordinates A vector of 3D particle positions.
         * @param cellLength The size of each cell in the spatial grid.
         * @param particleRadiusSq The squared radius for collision checks.
         * @param cellIdArray The array to store hashed cell IDs.
         * @param objectIdArray The array to store object metadata.
         */
        void initializeObjectandCellArray(
            const Coordinates& coordinates, 
            double cellLength, 
            double particleRadiusSq, 
            std::vector<int64_t>& cellIdArray, 
            std::vector<object_id>& objectIdArray);

        /**
         * Sorts the particles using a radix sort on their hashed cell IDs.
         *
         * @param cellIdArray The array of hashed cell IDs to be sorted.
         * @param objectIdArray The array of object metadata to be rearranged in sync with cellIdArray.
         */
        void radixSort(std::vector<int64_t>& cellIdArray, std::vector<object_id>& objectIdArray);

        /**
         * Calculates the total number of collisions between particles.
         *
         * @param coordinates A vector of 3D particle positions.
         * @param collisionDistSq The squared maximum collision distance.
         * @param cellIdArray The sorted array of hashed cell IDs.
         * @param objectIdArray The sorted array of object metadata.
         * @return The number of detected collisions.
         */
        int calculateNumberOfCollisions(const Coordinates& coordinates, double collisionDistSq, std::vector<int64_t>& cellIdArray, std::vector<object_id>& objectIdArray);

        /**
         * Calculates the squared Euclidean distance between two 3D points.
         *
         * @param p1 The first coordinate.
         * @param p2 The second coordinate.
         * @return The squared distance between the two points.
         */
        double calculateDistanceSq(const Coordinate& p1, const Coordinate& p2);

        /**
         * Decodes a 64-bit hashed cell ID into its 3D grid coordinates.
         *
         * @param hash The hashed cell ID.
         * @param grid_x Output for the X coordinate of the grid cell.
         * @param grid_y Output for the Y coordinate of the grid cell.
         * @param grid_z Output for the Z coordinate of the grid cell.
         */
        void decode_hash(const int64_t& hash, int64_t& grid_x, int64_t& grid_y, int64_t& grid_z);

        /**
         * Hashes the 3D grid cell coordinates into a 64-bit integer.
         *
         * @param x The X coordinate of the grid cell.
         * @param y The Y coordinate of the grid cell.
         * @param z The Z coordinate of the grid cell.
         * @return A 64-bit integer representing the hashed cell ID.
         */
        int64_t hash_coordinates(int64_t x, int64_t y, int64_t z);
};

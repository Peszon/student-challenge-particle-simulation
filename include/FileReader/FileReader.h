#pragma once

#include <string>   
#include <vector>   
#include <fstream>  // Open and read file.
#include <iomanip>  // std::linestream and .getline()

using Coordinate = std::vector<double>; // Represents a 3D coordinate as a vector of doubles
using Coordinates = std::vector<Coordinate>; // Represents a list of 3D coordinates

using std::istringstream;
using std::string;
using std::vector;

namespace DataReader 
{
    /**
     * Reads a file containing 3D coordinate data and returns a list of coordinates.
     *
     * @param FILE_NAME The path to the file containing the coordinates.
     * @return Coordinates A vector of 3D coordinate vectors.
     * @throws std::invalid_argument If the file cannot be opened or contains invalid data.
     */
    Coordinates read_positions_data(const std::string& FILE_NAME);
}

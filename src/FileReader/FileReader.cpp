#include "FileReader/FileReader.h"

namespace DataReader 
{
    /**
     * Reads a file containing 3D coordinate data and returns a list of coordinates.
     *
     * @param FILE_NAME The path to the file containing the coordinates.
     * @return Coordinates A vector of 3D coordinate vectors.
     * @throws std::invalid_argument If the file cannot be opened.
     */
    Coordinates read_positions_data(const string& FILE_NAME)
    {
        std::ifstream file(FILE_NAME); // Open the file for reading
        if (!file.is_open()) {
            throw std::invalid_argument( "Error: Could not open the file with the given filename." );
        }

        Coordinates coordinates;

        string line;
        while (std::getline(file, line)) { // Read the file line by line
            Coordinate coordinate;

            string value;
            istringstream line_stream(line);
            while (std::getline(line_stream, value, ' ')) {
                coordinate.push_back(stod(value));
            }
            coordinates.push_back(coordinate);
        }
        file.close(); 

        return coordinates;
    }
};

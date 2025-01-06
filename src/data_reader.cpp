#include "../include/data_reader.h"


vector<vector<double>> read_positions_data(string FILE_NAME)
{
    std::ifstream file(FILE_NAME); // Open the file for reading
    if (!file.is_open()) {
        throw std::invalid_argument( "Error: Could not open the file with the given filename." );
    }

    vector<vector<double>> coordinates;

    string line;
    while (std::getline(file, line)) { // Read the file line by line
        vector<double> coordinate;

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

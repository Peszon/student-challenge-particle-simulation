#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <string>
#include <iomanip>

#include "../include/data_reader.h"

int main()
{
    string DATA_FILE_NAME { "../data/positions.xyz" };

    vector<vector<double>> coordinates;
    coordinates = read_positions_data(DATA_FILE_NAME);

    // Display all parsed coordinates
    std::cout << "\nParsed Coordinates:\n";
    std::cout << std::setprecision(8);
    for (const auto& row : coordinates) {
        for (const auto& value : row) {
            std::cout << " " << value;
        }
        std::cout << "\n";
    }
    std::cout << "Numb of coordinates: " << coordinates.size() << std::endl;
      
    return 0;
}

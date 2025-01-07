#include "../include/data_reader.h"
#include "../include/brute_force.h"
#include "../include/timer.h"

using std::vector;
using std::string;

int main()
{
    string DATA_FILE_NAME { "../data/positions_large.xyz" };
    vector<vector<double>> coordinates;
    coordinates = read_positions_data(DATA_FILE_NAME);

    Timer t;
    
    for(int i { 0 }; i < 100; ++i) {
        calculate_number_colisions(coordinates, 0.05);
    } 
    std::cout << "Average time taken: " << t.elapsed()/100 << " seconds\n";
    return 0;
}

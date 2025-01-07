#include "../include/brute_force.h"

using std::vector;

double calculate_distance(const vector<double>& pos1, const vector<double>& pos2) {
    return sqrt(pow(pos1[0] - pos2[0], 2) + pow(pos1[1] - pos2[1], 2) + pow(pos1[2] - pos2[2], 2));
}

int calculate_number_colisions(const vector<vector<double>>& coordinates, double collision_distance) {
    int collitions {0};

    unsigned int i {1};
    unsigned int j {0}; 
    while (i < coordinates.size()) {
        while (j < i) {
            if (collision_distance > calculate_distance(coordinates[i], coordinates[j])) {
                ++collitions;
            }
            ++j;
        }
        ++i;
    }

    return collitions;
}

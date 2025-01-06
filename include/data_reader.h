#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <string>

using std::istringstream;
using std::string;
using std::vector;

/* **********************************************
 *
 * Reads the given datafiles and returns a 2D 
 * string vector with the file content.
 *
 * **********************************************/
vector<vector<double>> read_positions_data(string FILE_NAME);

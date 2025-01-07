#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <string>
#include <iomanip>


/* **********************************************
 *
 * Reads the given datafiles and returns a 2D 
 * string vector with the file content.
 * 
 * Input: - FILE_NAME: path of file.
 * Returns: 2D vector with positions of the particles. 
 *
 * **********************************************/
std::vector<std::vector<double>> read_positions_data(std::string& FILE_NAME);

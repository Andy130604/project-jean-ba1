/*
 * Mini-projet 3:
 */

#pragma once

#include <string>
#include <sstream>
#include <vector>

#include "pathsearch.hpp"

//------------- General utilities ----------
Network read_network(std::string network_filename);

//------------- Part 1 -------------
// Various helper function to print the datst structures
std::string to_string(const Reaction &reac, bool verbose = true);
std::string to_string(const Network &network, bool verbose = true);
std::string to_string(const Network &network, const AdjacencyGraph &graph);
std::string to_string(const BFS &bfs_data_structure);
std::string to_string(const Path &path);
std::string to_string(const Paths &paths);

std::string to_string(const Network &network, const Path &path, bool verbose = true);
std::string to_string(const Network &network, const Paths &paths, bool verbose = true);
std::string to_string(const Concentrations &concentrations);

//------------- Part 2 -------------
// helper function to read the concentrations from a file
Concentrations read_initial_concentrations(const Network &network, std::string initial_concentration_filename);

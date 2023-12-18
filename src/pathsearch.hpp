/*
 * Mini-projet 2023 : (C++11)
 */
#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <set>
#include <string>

const double V_IN = 5.0;
const double V_OUT = 1.0;
const double DELTA = 1e-8;

/*-----------------  TYPES AND DATA STRUCTURES    ------------*/

typedef std::string CompoundName;
typedef int ReactionID; // -1 used to mean "no reaction"
typedef int CompoundID; // -1 used  to mean "no compound", helps to stop recursive traversal

struct Reaction
{
    std::pair<CompoundID, CompoundID> compounds; // direction is defined by first->second
    double V_plus;
    double V_minus;
    double K_S;
    double K_P;
};

typedef std::vector<std::map<CompoundID, ReactionID>> b;

struct Network
{
    // the compounds are well ordered
    // vector index == CompoundID
    std::vector<CompoundName> compounds;
    // vector index == ReactionID
    std::vector<Reaction> reactions;
};
typedef std::vector<std::map<CompoundID, ReactionID>> AdjacencyGraph;

struct BFS
{
    CompoundID start;
    // vector index == CompoundID
    std::vector<std::vector<CompoundID>> parents;
    std::vector<int> distances;
};

typedef std::vector<ReactionID> Path;
typedef std::vector<Path> Paths;
typedef std::map<CompoundID, double> Concentrations;

///------------- Part 1 -------------
/*!
 * @brief finds the ID of a compound in a network given its name
 * @return -1 if no compound is found
 */
CompoundID find_compoundID(const Network &network, CompoundName name);

/*!
 * @brief builds the adjacency graph of a network
 */
AdjacencyGraph build_adjacency_graph(const Network &network);

/*!
 * @brief finds the Id of the relation linking a compound to another
 * @return -1 is no reaction is found
 */
ReactionID find_reactionID(const AdjacencyGraph &graph,
                           CompoundID index1,
                           CompoundID index2);

/*!
 * @brief performs breadth-first search on a graph starting fron the compound Id
 * and stores the result of the traversal in a BFS data structure
 */
BFS bfs(const AdjacencyGraph &graph, CompoundID start);

///------------- Part 2 -------------

/*!
 * @brief finds a shortest path between a source to a destination compound using breadth-first search
 * @param srcID  Id of the source compound
 * @param destID Id of the destination compound
 */
Path find_shortest_path(const AdjacencyGraph &graph, CompoundID srcID, CompoundID destID); // ~30 lines

/*!
 * @brief Rercursively finds all the shortest paths from a source to a destination using a BFS result to iterate in the reverse direction (dest -> src)
 * @param graph Adjacency graph of the whole network
 * @param result BFS structure created by BFS algorithm
 * @param src Node ID of the source
 * @param dest Node ID of the destination
 * @param currentPath Path the recursive function is working on
 * @param allPaths Object to store all paths found
 */
void recursive_find_paths(const AdjacencyGraph &graph, BFS &result, CompoundID &src, CompoundID &dest, Path &currentPath, Paths &allPaths);

/*!
 * @brief finds all shortest path between a source to a destination compound using breadth-first search
 * @param srcID  Id of the source compound
 * @param destID Id of the destination compound
 */
Paths find_all_shortest_paths(const AdjacencyGraph &graph, CompoundID srcID, CompoundID destID); // change to adj graph

//------------- Part 3 -------------

/*!
 * @brief michaelis reversible rate
 * @param R the reaction
 * @param S reactive concentration
 * @param P product concentration
 * @return the rate of the reaction
 */
double michaelis_reversible_rate(const Reaction &R, const double S, const double P);

/*!
 * @brief Creates a vector of all the coumpound ID's of the path
 * @param network Whole network - Used to find reactions based on their ID
 * @param path The path to create the compound ID vector
 * @return std::vector<CompoundID> list of compoundID's
 */
std::vector<CompoundID> compute_coumpound_path(const Network &network, const Path &path);

/*!
 * @brief Computes the new concentrations after a time period
 * @param network The whole network
 * @param path The path to compute the new concentrations on
 * @param c_in The old concentrations
 * @param dt The time period
 * @return The new concentrations computed
 */
Concentrations euler_implicite(const Network &network, const Path &path, const Concentrations c_in, double dt);

/*!
 * @brief Checks if a pair of old and new concentrations are stable or not
 * @param c_in The old concentrations
 * @param c_out The new concentrations
 * @return A boolean true if stable false otherwise
 */
bool checkStable(const Concentrations &c_in, const Concentrations &c_out);

/*!
 * @brief computes steady state concentrations in a single path : applies euler_implicite until convergence
 * @param initial_concentrations intitial concentrations in the network's compounds
 * @return the steady state concentrations in the compounds of the given path
 */
Concentrations compute_ss_concentration(const Network &network, const Path &path, const Concentrations &initial_concentrations, double dt = 1e-3);

/*!
 * @brief computes the smallest michaelis_reversible_rate in the given path
 * (the michaelis_reversible_rate is computed between each pair of compounds of the path, and the smallest is returned)
 * @param ss_concentrations steady state concentrations in the compounds of the network
 *
 * @return the smallest (slowest) rate in a given path
 */
double compute_path_rate(const Network &network, const Path &path, const Concentrations &ss_concentrations);

/*!
 * @brief computes fastest path among a set of given paths
 * (the one with the highest compute_path_rate)
 * @param paths the set of paths to consider
 * @param initial_concentrations intitial concentrations in the compounds of the network
 * @param dt the time step used to compute the steady state concentrations
 * @return a fastest path
 */
Path find_fastest_path(const Network &network, const Paths &paths, const Concentrations &initial_concentrations, double dt);

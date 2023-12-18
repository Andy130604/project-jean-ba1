/*
 * Mini-projet 3 : (C++11)
 */

#include <fstream>
#include <iomanip>
#include "utils.hpp"
#include "pathsearch.hpp"
#include <cmath>
#include <array>
#include <queue>
#include <list>
#include <climits>
#include <algorithm>
#include <cassert>
#include <map>
#include <set>

//==================================================================
//                              PART 1
//==================================================================
CompoundID find_compoundID(const Network &network, CompoundName vertex)
{
    CompoundID compoundID = -1;
    size_t i = 0;
    while (compoundID == -1 && i < network.compounds.size())
    {
        if (network.compounds[i] == vertex)
        {
            compoundID = (CompoundID)i;
        }
        i++;
    }

    return compoundID;
}

AdjacencyGraph build_adjacency_graph(const Network &network)
{
    AdjacencyGraph graph(network.compounds.size());
    for (size_t i = 0; i < network.reactions.size(); ++i)
    {
        std::pair<CompoundID, CompoundID> compounds = network.reactions[i].compounds;
        graph[compounds.first].insert({compounds.second, (ReactionID)i});
        graph[compounds.second].insert({compounds.first, (ReactionID)i});
    }

    return graph;
}

BFS bfs(const AdjacencyGraph &adjacency_graph, CompoundID start)
{
    BFS result;
    size_t size = adjacency_graph.size();
    result.start = start;
    result.parents.resize(size);
    result.parents[start] = {-1};
    result.distances.resize(size);
    std::fill(result.distances.begin(), result.distances.end(), INT_MAX);
    result.distances[start] = 0;

    std::queue<CompoundID> queue;
    queue.push(start);
    while (!queue.empty())
    {
        CompoundID currentNode = queue.front();
        queue.pop();
        for (std::pair<CompoundID, ReactionID> pair : adjacency_graph[currentNode])
        {
            if (result.distances[pair.first] > result.distances[currentNode] + 1)
            {
                result.distances[pair.first] = result.distances[currentNode] + 1;
                queue.push(pair.first);
                result.parents[pair.first] = {currentNode};
            }
            else if (result.distances[pair.first] == result.distances[currentNode] + 1)
            {
                result.parents[pair.first].push_back(currentNode);
            }
        }
    }

    return result;
}

ReactionID find_reactionID(const AdjacencyGraph &adjacency_graph, CompoundID index1, CompoundID index2)
{
    ReactionID reactionID = -1;
    std::map<CompoundID, ReactionID> map = adjacency_graph[index1];
    std::map<CompoundID, ReactionID>::iterator it = map.find(index2);
    if (it != map.end())
    {
        reactionID = it->second;
    }

    return reactionID;
}

//==================================================================
//                              PART 2
//==================================================================

Path find_shortest_path(const AdjacencyGraph &graph, CompoundID srcID, CompoundID destID)
{
    BFS result = bfs(graph, srcID);
    Path path;
    CompoundID currentNode = destID;
    CompoundID firstParent = result.parents[currentNode][0];
    while (firstParent != -1)
    {
        path.push_back(find_reactionID(graph, currentNode, firstParent));
        currentNode = firstParent;
        firstParent = result.parents[currentNode][0];
    }
    std::reverse(path.begin(), path.end());
    return path;
}

void recursive_find_paths(const AdjacencyGraph &graph, BFS &result, CompoundID &src, CompoundID &dest, Path &currentPath, Paths &allPaths)
{
    if (src == dest)
    {
        allPaths.push_back(currentPath);
        return;
    }

    for (CompoundID parent : result.parents[dest])
    {
        currentPath.push_back(find_reactionID(graph, dest, parent));
        recursive_find_paths(graph, result, src, parent, currentPath, allPaths);
        currentPath.pop_back();
    }
}

Paths find_all_shortest_paths(const AdjacencyGraph &graph, CompoundID srcID, CompoundID destID)
{
    BFS result = bfs(graph, srcID);
    Paths allPaths;
    Path currentPath;

    recursive_find_paths(graph, result, srcID, destID, currentPath, allPaths);

    for (Path &path : allPaths)
    {
        std::reverse(path.begin(), path.end());
    }

    return allPaths;
}

//==================================================================
//                              PART 3
//==================================================================

// michaelis reversible rate
double michaelis_reversible_rate(const Reaction &R, const double S, const double P)
{
    return (R.V_plus * (S / R.K_S) - R.V_minus * (P / R.K_P)) / (1 + S / R.K_S + P / R.K_P);
}

std::vector<CompoundID> compute_coumpound_path(const Network &network, const Path &path)
{
    std::vector<CompoundID> compound_path;
    if (path.size() == 1)
    {
        Reaction reaction = network.reactions[path[0]];
        compound_path.push_back(reaction.compounds.first);
        compound_path.push_back(reaction.compounds.second);
    }
    else
    {
        for (size_t i = 0; i < path.size() - 1; ++i)
        {
            Reaction reaction = network.reactions[path[i]];
            Reaction next = network.reactions[path[i + 1]];
            CompoundID left = reaction.compounds.first;
            CompoundID right = reaction.compounds.second;
            if (left == next.compounds.first || left == next.compounds.second)
            {
                if (i == 0)
                {
                    compound_path.push_back(right);
                }
                compound_path.push_back(left);
            }
            else
            {
                if (i == 0)
                {
                    compound_path.push_back(left);
                }
                compound_path.push_back(right);
            }
        }

        Reaction reaction = network.reactions[path[path.size() - 1]];
        Reaction prev = network.reactions[path[path.size() - 2]];
        CompoundID left = reaction.compounds.first;
        CompoundID right = reaction.compounds.second;
        if (left == prev.compounds.first || left == prev.compounds.second)
        {
            compound_path.push_back(right);
        }
        else
        {
            compound_path.push_back(left);
        }
    }

    return compound_path;
}

Concentrations euler_implicite(const Network &network, const Path &path, const Concentrations c_in, double dt)
{
    std::vector<double> rates;
    std::vector<CompoundID> compound_path = compute_coumpound_path(network, path);
    for (size_t i = 0; i < path.size(); ++i)
    {
        rates.push_back(michaelis_reversible_rate(network.reactions[path[i]], c_in.find(compound_path[i])->second, c_in.find(compound_path[i + 1])->second));
    }

    Concentrations c_out;
    size_t i = 0;
    for (auto it = c_in.begin(); it != c_in.end(); it++)
    {
        double rateOfChange;
        if (it->first == compound_path[0])
        {
            rateOfChange = V_IN * (1.0 - it->second) - rates[0];
        }
        else if (it->first == compound_path[compound_path.size() - 1])
        {
            rateOfChange = rates[rates.size() - 1] - (it->second * V_OUT);
        }
        else
        {
            rateOfChange = rates[i] - rates[i + 1];
            i++;
        }
        double newConcentration = it->second + dt * rateOfChange;
        if (newConcentration < 0)
        {
            newConcentration = 0.0;
        }
        c_out.insert({it->first, newConcentration});
    }

    return c_out;
}

bool checkStable(const Concentrations &c_in, const Concentrations &c_out)
{
    bool stable = true;
    auto prevIt = c_in.begin();
    auto nextIt = c_out.begin();
    while (prevIt != c_in.end() && stable)
    {
        if (fabs(nextIt->second - prevIt->second) / nextIt->second >= DELTA)
        {
            stable = false;
        }
        prevIt++;
        nextIt++;
    }
    return stable;
}

Concentrations compute_ss_concentration(const Network &network, const Path &path, const Concentrations &initial_concentrations, double dt)
{
    Concentrations c_in, c_out;
    for (size_t i = 0; i < path.size(); ++i)
    {
        Reaction reaction = network.reactions[path[i]];
        CompoundID left = reaction.compounds.first;
        CompoundID right = reaction.compounds.second;
        if (c_in.find(left) == c_in.end())
        {
            c_in.insert({left, initial_concentrations.find(left)->second});
        }
        if (c_in.find(right) == c_in.end())
        {
            c_in.insert({right, initial_concentrations.find(right)->second});
        }
    }

    c_out = euler_implicite(network, path, c_in, dt);
    bool stable = checkStable(c_in, c_out);
    while (!stable)
    {
        c_in = c_out;
        c_out = euler_implicite(network, path, c_in, dt);
        stable = checkStable(c_in, c_out);
    }
    return c_out;
}

double compute_path_rate(const Network &network, const Path &path, const Concentrations &ss_concentrations)
{
    std::vector<CompoundID> compound_path = compute_coumpound_path(network, path);
    double minRate = INT_MAX;
    for (size_t i = 0; i < path.size(); ++i)
    {
        double rate = michaelis_reversible_rate(network.reactions[path[i]], ss_concentrations.find(compound_path[i])->second, ss_concentrations.find(compound_path[i + 1])->second);
        minRate = rate < minRate ? rate : minRate;
    }

    return minRate;
}

Path find_fastest_path(const Network &network, const Paths &paths, const Concentrations &initial_concentrations, double dt)
{
    double maxPathRate = INT_MIN;
    Path bestPath;
    for (Path path : paths)
    {
        Concentrations ss_concentrations = compute_ss_concentration(network, path, initial_concentrations, dt);
        double pathRate = compute_path_rate(network, path, ss_concentrations);
        if (pathRate > maxPathRate)
        {
            maxPathRate = pathRate;
            bestPath = path;
        }
    }

    return bestPath;
}

/*
 * Mini-projet 3
 */
#include "utils.hpp"
#include <fstream>
#include <assert.h>
// #include <filesystem>
#include <iostream>

double get_random_value()
{
    return (rand() % 1000) / 1000.0;
}

Network read_network(std::string network_filename)
{
    std::ifstream file("../" + network_filename);
    // std::ifstream file(network_filename);
    Network network;
    int cId(0);
    int rId(0);

    std::map<std::string, int> compoundIDs;

    if (file)
    {
        std::string line;
        getline(file, line);
        while (line[0] != '-')
        {
            CompoundName c(line);
            network.compounds.push_back(c);
            compoundIDs[c] = cId;
            cId++;
            getline(file, line);
        }
        while (getline(file, line))
        {
            if (!line.empty() and line[0] != '#' and line[0] != '-')
            {
                Reaction reac;
                reac.compounds.first = compoundIDs.at(line);

                getline(file, line);
                reac.compounds.second = compoundIDs.at(line);

                getline(file, line);
                reac.V_plus = stod(line);

                getline(file, line);
                reac.V_minus = stod(line);

                getline(file, line);
                reac.K_S = stod(line);

                getline(file, line);
                reac.K_P = stod(line);

                network.reactions.push_back(reac);
                rId++;
            }
        }
    }
    else
    {
        std::cout << "File not found: " << network_filename << std::endl;
        assert(false);
    }
    return network;
}

std::string to_string(const Reaction &r, bool verbose)
{
    std::stringstream ss;
    if (verbose)
        ss << "Compounds: ";
    ss << r.compounds.first << " -> ";
    ss << r.compounds.second << std::endl;
    if (verbose)
    {
        //  ss << "V_max = " << r.V_max << std::endl;
        //  ss << "K_eq  = " << r.K_eq << std::endl;
        //  ss << "KM_S  = " << r.KM_S << std::endl;
        //  ss << "Reversible = " << r.reversible << std::endl;
        //  ss << "KM_P  = " << r.KM_P << std::endl;

        ss << "V_plus = " << r.V_plus << std::endl;
        ss << "V_minus = " << r.V_minus << std::endl;
        ss << "K_S  = " << r.K_S << std::endl;
        ss << "K_P  = " << r.K_P << std::endl;
    }
    return ss.str();
}

std::string to_string(const BFS &bfs_data_structure)
{
    std::stringstream ss;
    ss << "start:" << bfs_data_structure.start << "\n";
    auto parents = bfs_data_structure.parents;
    for (size_t i(0); i < parents.size(); ++i)
    {
        ss << "Node" << i << ": {";
        for (auto parent : parents[i])
        {
            ss << parent << " ";
        }
        ss << "}, ";
        ss << "distance : " << bfs_data_structure.distances[i] << "\n";
    }
    return ss.str();
}

std::string to_string(const Network &network, const AdjacencyGraph &graph)
{
    std::stringstream ss;

    for (size_t i(0); i < graph.size(); ++i)
    {
        ss << i << ":" << network.compounds[i] << ": {";
        for (auto adj : graph[i])
        {
            auto compound_index(adj.first);
            auto reaction_id(adj.second);
            // ss << compound_index << ":" << network.compounds[compound_index] << ",";
            ss << "{" << compound_index << "," << reaction_id << "}";
        }
        ss << "}\n";
    }
    return ss.str();
}
std::string to_string(const Network &network, bool verbose)
{
    std::stringstream ss;

    for (size_t i(0); i < network.compounds.size(); i++)
    {
        if (verbose)
            ss << "Compound ID = " << i << " Name = ";
        ss << network.compounds.at(i) << std::endl;
    }
    for (size_t i(0); i < network.reactions.size(); i++)
    {
        if (verbose)
            ss << "--- Reaction ID = " << i << " --- " << std::endl;
        ss << to_string(network.reactions.at(i), verbose);
    }
    return ss.str();
}

std::string to_string(const Path &path)
{
    std::stringstream ss;
    ss << "Reactions: ";
    for (auto r : path)
    {
        ss << r << " ";
    }
    ss << "\n";
    return ss.str();
}

std::string to_string(const Paths &paths)
{
    std::stringstream ss;
    for (size_t i(0); i < paths.size(); ++i)
    {
        ss << "Path Number " << i + 1 << " - " << to_string(paths[i]);
    }
    return ss.str();
}

std::string to_string(const Network &network, const Path &path, bool verbose)
{
    std::stringstream ss;
    for (auto r : path)
    {
        if (verbose)
        {
            ss << "Reaction (" << r << "): ";
            CompoundID first = network.reactions.at(r).compounds.first;
            CompoundID second = network.reactions.at(r).compounds.second;
            ss << network.compounds.at(first) << " (" << first << ") -> ";
            ss << network.compounds.at(second) << " (" << second << ")" << std::endl;
        }
        ss << to_string(network.reactions.at(r), verbose);
    }
    return ss.str();
}

std::string to_string(const Network &network, const Paths &paths, bool verbose)
{
    std::stringstream ss;
    size_t k(0);
    for (auto path : paths)
    {
        ss << "Path " << ++k << ":\n";
        ss << to_string(network, path, verbose);
    }
    return ss.str();
}
Concentrations read_initial_concentrations(const Network &network, std::string filename)
{
    std::ifstream file("../" + filename);

    // std::ifstream file(filename);
    Concentrations concentrations;

    if (file)
    {
        std::string line;

        while (getline(file, line))
        {
            if (!line.empty())
            {
                for (size_t i(0); i < line.size(); ++i)
                {
                    if (line[i] == '=')
                    {
                        std::string name = line.substr(1, i - 2); // fonction remove for the '[' char
                        double concentration = stod(line.substr(i + 1, line.size() - 1));
                        int id = find_compoundID(network, name);
                        if (id >= 0)
                            concentrations[id] = concentration;
                        //          if ( id >= 0) concentrations[id] = get_random_value();
                        else
                            std::cerr << "Component " << name << " does not exist." << std::endl;
                        break;
                    }
                }
            }
        }
    }
    else
    {
        std::cout << "File not found: " << filename << std::endl;
        assert(false);
    }
    return concentrations;
}

std::string to_string(const Concentrations &concentrations)
{
    std::stringstream ss;
    for (auto element : concentrations)
    {
        ss << element.first << " = " << element.second << std::endl;
    }
    return ss.str();
}

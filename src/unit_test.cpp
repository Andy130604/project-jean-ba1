#include <algorithm>
#include <cmath> // std::fabs
#include <iomanip>
#include <iostream> // std::cerr, std::endl
#include <limits>   // std::numeric_limits
#include <sstream>
#include <regex>
#include "pathsearch.hpp"
#include "unit_test.hpp"
#include "utils.hpp"

using namespace std;

/**
 * Utility Functions
 */

const AdjacencyGraph SEVEN_PATH_ADJACENCY = {{{{1, 0}, {3, 2}, {6, 8}},
                                              {{0, 0}, {2, 1}, {3, 5}, {5, 6}},
                                              {{1, 1}, {4, 4}},
                                              {{0, 2}, {1, 5}, {4, 3}, {6, 9}},
                                              {{2, 4}, {3, 3}, {5, 7}},
                                              {{1, 6}, {4, 7}},
                                              {{0, 8}, {3, 9}}}};

// Check if two floats are equal within precision tolerance
bool equal(double x, double y, double precision) { return (std::fabs(x - y) < precision); }

// Utility function to print out the header (name)
// Used by test cases to print out the test name
void print_header(std::string const &header)
{
    std::cerr << "---------------" << std::endl;
    std::cerr << header << std::endl;
    std::cerr << "---------------" << std::endl;
}

// Check if two doubles are equal, and print the results out
void check_equal(double expected, double computed)
{
    if (equal(computed, expected))
    {
        std::cerr << "[Passed]" << std::endl
                  << std::endl;
        return;
    }
    std::cerr << "[Failed]" << std::endl;
    std::cerr << "   expected: " << expected << std::endl;
    std::cerr << "   computed: " << computed << std::endl;
}

void check_equal(int expected, int computed)
{
    if (expected == computed)
    {
        std::cerr << "[Passed]" << std::endl
                  << std::endl;
        return;
    }
    std::cerr << "[Failed]" << std::endl;
    std::cerr << "   expected: " << expected << std::endl;
    std::cerr << "   computed: " << computed << std::endl;
}
void check_equal(const Path &expected, const Path &computed)
{
    if (expected == computed)
    {
        std::cerr << "[Passed]" << std::endl
                  << std::endl;
        return;
    }
    std::cerr << "[Failed]" << std::endl;
    std::cerr << "   expected: " << to_string(expected) << std::endl;
    std::cerr << "   computed: " << to_string(computed) << std::endl;
}

void check_equal(const Paths &expected, const Paths &computed)
{
    if (expected == computed)
    {
        std::cerr << "[Passed]" << std::endl
                  << std::endl;
        return;
    }
    std::cerr << "[Failed]" << std::endl;
    std::cerr << "   expected: " << to_string(expected) << std::endl;
    std::cerr << "   computed: " << to_string(computed) << std::endl;
}

void check_equal(const AdjacencyGraph &expected, const AdjacencyGraph &computed, const Network &network)
{
    if (expected == computed)
    {
        std::cerr << "[Passed]" << std::endl
                  << std::endl;
        return;
    }
    std::cerr << "[Failed]" << std::endl;
    std::cerr << "   expected: " << to_string(network, expected) << std::endl;
    std::cerr << "   computed: " << to_string(network, computed) << std::endl;
}

void check_equal(const BFS &expected, const BFS &computed)
{
    if (expected == computed)
    {
        std::cerr << "[Passed]" << std::endl
                  << std::endl;
        return;
    }
    std::cerr << "[Failed]" << std::endl;
    std::cerr << "   expected: " << to_string(expected) << std::endl;
    std::cerr << "   computed: " << to_string(computed) << std::endl;
}

// This `operator==` is a is a technique called "operator overloading"
// This implements the '==' operator for the type BFS so that we can do comparisons to see if
// two BFS results are equal. If you want to know more take a look here:
// https://en.cppreference.com/w/cpp/language/operators
bool operator==(const BFS &a, const BFS &b)
{
    if (a.parents.size() != b.parents.size())
        return false;
    if (a.distances.size() != b.distances.size())
        return false;
    for (size_t i(0); i < a.parents.size(); ++i)
    {
        std::set<int> a_as_set;
        std::set<int> b_as_set;
        copy(a.parents[i].begin(), a.parents[i].end(), inserter(a_as_set, a_as_set.end()));
        copy(a.parents[i].begin(), a.parents[i].end(), inserter(b_as_set, b_as_set.end()));
        if (a_as_set != b_as_set)
            return false;
    }
    return a.start == b.start;
}
bool operator==(const AdjacencyGraph &a, const AdjacencyGraph &b)
{
    // Note: C++ maps are ordered by key.
    if (a.size() != b.size())
        return false;
    for (size_t i(0); i < a.size(); ++i)
    {
        if (a[i] != b[i])
            return false;
    }
    return true;
}

/**
 * Tests
 *
 * These tests are not the tests that will be used to grade
 * your project, but are here to help you debug.
 * You can also add additional tests here.
 */

void test_find_compoundId()
{
    print_header("test_find_compoundId");
    Network network = read_network("data/7paths.txt");
    std::cerr << "Testing with network 7paths.txt " << std::endl;
    check_equal(find_compoundID(network, "C00025"), 0);
    check_equal(find_compoundID(network, "C03912"), 1);
    check_equal(find_compoundID(network, "C00148"), 2);
    check_equal(find_compoundID(network, "C01165"), 3);
    check_equal(find_compoundID(network, "C00077"), 4);
    check_equal(find_compoundID(network, "C00062"), 5);
    check_equal(find_compoundID(network, "C03287"), 6);
}

void test_build_adjacency_graph()
{
    print_header("test_build_adjacency_graph");
    Network network = read_network("data/7paths.txt");
    AdjacencyGraph computed_graph(build_adjacency_graph(network));
    AdjacencyGraph expected_graph(SEVEN_PATH_ADJACENCY);
    check_equal(expected_graph, computed_graph, network);
}
void test_bfs()
{
    print_header("test_bfs");
    AdjacencyGraph graph(SEVEN_PATH_ADJACENCY);
    BFS computed(bfs(graph, 0));
    BFS expected({0,
                  {{-1}, {0}, {1}, {0}, {3}, {1}, {0}},
                  {0, 1, 2, 1, 2, 2, 1}});
    check_equal(expected, computed);
}
void test_find_reactionId()
{
    print_header("test_find_reactionId");
    AdjacencyGraph graph(SEVEN_PATH_ADJACENCY);
    check_equal(find_reactionID(graph, 0, 3), 2);
    check_equal(find_reactionID(graph, 0, 1), 0);
    check_equal(find_reactionID(graph, 0, 6), 8);
    check_equal(find_reactionID(graph, 1, 2), 1);
}
void test_find_shortest_path()
{
    print_header("test_find_shortest_path");
    AdjacencyGraph graph(SEVEN_PATH_ADJACENCY);
    Path computed(find_shortest_path(graph, 0, 2));
    Path expected({0, 1});
    check_equal(expected, computed);
}

void test_find_all_shortest_paths()
{
    print_header("test_find_all_shortest_paths");
    AdjacencyGraph graph(SEVEN_PATH_ADJACENCY);
    Paths computed(find_all_shortest_paths(graph, 3, 2));
    Paths expected({{5, 1}, {3, 4}});
    check_equal(expected, computed);
}

void test_michaelis_reversible_rate()
{
    print_header("test_michaelis_reversible_rate");
    Network network = read_network("data/7paths.txt");
    std::cerr << "Testing with network 7paths.txt " << std::endl;
    check_equal(2.20623, michaelis_reversible_rate(network.reactions[0], 0.5, 0.1));
    check_equal(0.0535988, michaelis_reversible_rate(network.reactions[0], 0.4, 1.0));
    check_equal(-1.14094, michaelis_reversible_rate(network.reactions[0], 0.05, 0.8));

    check_equal(0.700876, michaelis_reversible_rate(network.reactions[1], 0.5, 0.1));
    check_equal(-0.814342, michaelis_reversible_rate(network.reactions[1], 0.4, 1.0));
    check_equal(-1.46621, michaelis_reversible_rate(network.reactions[1], 0.05, 0.8));
}

void test_compute_ss_concentration()
{
    print_header("test_compute_ss_concentration");
    Network network = read_network("data/7paths.txt");
    Concentrations initial = read_initial_concentrations(network, "data/7paths_concentrations.txt");
    std::cerr << "Testing with network 7paths.txt " << std::endl;
    Paths paths({{5, 1}, {3, 4}});

    Path path = paths[0];
    Concentrations cs0 = compute_ss_concentration(network, path, initial, 1e-3);
    check_equal(0.755165, cs0[1]);
    check_equal(0.416494, cs0[2]);
    check_equal(0.916702, cs0[3]);

    // Testing another path
    path = paths[1];
    Concentrations cs1 = compute_ss_concentration(network, path, initial, 1e-3);
    check_equal(0.337187, cs1[2]);
    check_equal(0.932562, cs1[3]);
    check_equal(0.728524, cs1[4]);
}

void test_compute_path_rate()
{
    print_header("test_compute_path_rate");
    Network network = read_network("data/7paths.txt");
    Concentrations initial = read_initial_concentrations(network, "data/7paths_concentrations.txt");
    std::cerr << "Testing with network 7paths.txt " << std::endl;
    Paths paths({{5, 1}, {3, 4}});

    Concentrations ss = {{0, 0.5}, {1, 0.4}, {2, 0.3}, {3, 0.7}, {4, 0.8}, {5, 0.6}, {6, 0.2}};

    Path path = paths[0];
    check_equal(0.163442, compute_path_rate(network, path, ss));

    path = paths[1];
    check_equal(-0.477099, compute_path_rate(network, path, ss));
}

void test_find_fastest_path()
{
    print_header("test_find_fastest_path");
    Network network = read_network("data/7paths.txt");
    Concentrations initial = read_initial_concentrations(network, "data/7paths_concentrations.txt");
    std::cerr << "Testing with network 7paths.txt " << std::endl;
    /*
    AdjacencyGraph graph ( {
                                {{1,0},{3,2},{6,8}},
                                {{0,0}, {2,1},{3,5},{5,6}},
                                {{1,1},{4,4}},
                                {{0,2},{1,8},{4,3},{6,9}},
                                {{2,4},{3,3},{5,7}},
                                {{1,6},{4,7}},
                                {{0,9},{3,9}}
                                });
  */
    Paths paths({{5, 1}, {3, 4}});

    Path fastest_path = find_fastest_path(network, paths, initial, 1e-2);

    check_equal({5, 1}, fastest_path);
}

// Run all of the unit tests
void run_unit_tests(int part)
{
    if (part == 1)
    {
        // TASK1
        test_find_compoundId();
        test_build_adjacency_graph();
        test_find_reactionId();
        test_bfs();
    }
    else if (part == 2)
    {
        // TASK2
        test_find_shortest_path();
        test_find_all_shortest_paths();
    }
    else if (part == 3)
    {
        // TASK 3
        test_michaelis_reversible_rate();
        test_compute_ss_concentration();
        test_compute_path_rate();
        test_find_fastest_path();
    }
    else
    {
        std::cerr << "Part should be either 1, 2, or 3. Provided part = " << part << std::endl;
    }
}

#include <algorithm>
#include <cmath> // std::fabs
#include <iomanip>
#include <iostream> // std::cerr, std::endl
#include <limits>   // std::numeric_limits
#include <sstream>
#include <regex>
#include "pathsearch.hpp"
#include "utils.hpp"

using namespace std;

/**
 * Utility Functions
 */

// Check if two floats are equal within precision tolerance
bool equal(double x, double y, double precision = 1e-5);

// Utility function to print out the header (name)
// Used by test cases to print out the test name
void print_header(std::string const &header);

// Check if two doubles are equal, and print the results out
void check_equal(double expected, double computed);

void check_equal(int expected, int computed);

bool operator==(const BFS &a, const BFS &b);

bool operator==(const AdjacencyGraph &a, const AdjacencyGraph &b);

/**
 * Tests
 *
 * These tests are not the tests that will be used to grade
 * your project, but are here to help you debug.
 * You can also add additional tests here.
 */
void test_find_compoundId();

void test_michaelis_reversible_rate();

void test_compute_ss_concentration();

void test_compute_path_rate();

void run_unit_tests(int part = 1);

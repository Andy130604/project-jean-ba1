all: pathsearch

pathsearch: utils.hpp utils.cpp main.cpp pathsearch.cpp pathsearch.hpp unit_test.hpp unit_test.cpp
	c++ -std=c++14 -Wall main.cpp utils.cpp pathsearch.cpp unit_test.cpp -o pathsearch

run: pathsearch
	./pathsearch

clean:
	rm -f pathsearch

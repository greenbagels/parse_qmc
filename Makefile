CXX=c++
CFLAGS=-std=c++17 -Wall -Wextra -Werror -Wfatal-errors
# Need to link filesystem library
LDFLAGS=-lstdc++fs

default:
	$(CXX) $(CFLAGS) parse.cpp $(LDFLAGS) -o parse

# Copyright (c) 2019, VERSO.
# All rights reserved.

# Variables.
CXX ?= g++
CXXFLAGS = -MMD -MP -I. -std=c++14 -Wextra -Wpedantic -Wall -O3
LDLIBS = -lpthread

# Using all cpp files in current directory.
MAIN = ../bin/posman
SRC = $(wildcard *.cpp)\
			$(wildcard ./structures/*.cpp)\
			$(wildcard ./structures/graph/*.cpp)\

OBJ = $(SRC:.cpp=.o)
DEPS = $(SRC:.cpp=.d)

# Main target.
all : $(MAIN)

$(MAIN) : $(OBJ) main.o
	mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

# Building .o files.
%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

-include ${DEPS}

clean :
	$(RM) $(OBJ) $(DEPS)
	$(RM) $(MAIN)

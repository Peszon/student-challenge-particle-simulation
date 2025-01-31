# Makefile for building the C++ project in production mode

# Automatically use 16 parallel jobs if not specifie

# Compiler to use
CXX := /usr/bin/g++

# Compiler flags for production
CXXFLAGS := -fdiagnostics-color=always \
            -fopenmp \
            -std=c++20 \
            -O3 \
            -Iinclude

# Source files
SRCS := src/main.cpp \
        src/CollisionFramework/CollisionFramework.cpp \
        src/Collision/Algorithms/BruteForceAlgorithm.cpp \
        src/Collision/Algorithms/SpatialSubdivisionAlgorithm.cpp \
        src/Collision/Algorithms/SpatialSubdivisionOpenMPAlgorithm.cpp \
        src/FileReader/FileReader.cpp

# Object files (auto-generated from source files)
OBJS := $(SRCS:.cpp=.o)

# Directory for the executable
BIN_DIR := bin

# Name of the final executable
TARGET := $(BIN_DIR)/output

# Default target
all: $(TARGET)

# Rule to link object files into the final executable
$(TARGET): $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Generic rule to compile .cpp files to .o object files
%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build artifacts
clean:
	rm -f $(OBJS) $(TARGET)

# Phony targets to prevent conflicts with files of the same name
.PHONY: all clean

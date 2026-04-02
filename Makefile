# Paths to Z3
#Z3_DIR := ./z3
#Z3_INCLUDE := $(Z3_DIR)/include
#Z3_BIN := $(Z3_DIR)/bin

# Compiler and flags
CXX := g++
CXXFLAGS := -O2 -fopenmp -Wall -Werror -pedantic #-I$(Z3_INCLUDE) -I. 
#LDFLAGS := -L$(Z3_BIN) -lz3

# Source and target
SRC := HRSA_test.cpp cyclotomic_int9.cpp Z9chi.cpp householder_search.cpp decompose.cpp
TARGET := HRSA_tester

# Build rule
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@ #$(LDFLAGS)

# Run rule
run: $(TARGET)
	./$(TARGET)

# Clean rule
clean:
	rm -f $(TARGET)

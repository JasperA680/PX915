# Compiler
FC = gfortran

# NetCDF flags (Fortran wrapper + underlying C library)
NC_FFLAGS := $(shell nf-config --fflags)
NC_FLIBS   := $(shell nf-config --flibs) $(shell nc-config --libs)

# Compiler flags
FFLAGS = -Wall -O2 $(NC_FFLAGS)

# Directories
SRC_DIR = src/fortran
BUILD_DIR = build

# Source files (existing single-lattice TASEP)
MODEL_SRC = $(SRC_DIR)/tasep.f90
SIM_SRC   = $(SRC_DIR)/simulation.f90
IO_SRC    = $(SRC_DIR)/io.f90
TEST_SRC  = $(SRC_DIR)/test_simulation.f90

# Network-model sources (new)
VEHICLE_SRC  = $(SRC_DIR)/vehicle.f90
NETWORK_SRC  = $(SRC_DIR)/road_network.f90

# Executables
TEST_EXE       = $(BUILD_DIR)/test_simulation
TEST_TYPES_EXE = $(BUILD_DIR)/test_types

# Default target
all: test test_types

# Build existing test program
test: $(TEST_EXE)

$(TEST_EXE): $(MODEL_SRC) $(SIM_SRC) $(IO_SRC) $(TEST_SRC)
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) $^ -o $@ $(NC_FLIBS)

# Phase 1 smoke test
test_types: $(TEST_TYPES_EXE)

$(TEST_TYPES_EXE): $(VEHICLE_SRC) $(NETWORK_SRC) tests/test_types.f90
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) $^ -o $@


# Run the existing simulation test
run: test
	./$(TEST_EXE)

# Run the Phase 1 smoke test
run_test_types: test_types
	./$(TEST_TYPES_EXE)

# Clean build files
clean:
	rm -rf $(BUILD_DIR)
	rm -f $(SRC_DIR)/*.mod *.mod

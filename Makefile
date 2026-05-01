# Compiler
FC = gfortran

# NetCDF flags (Fortran wrapper + underlying C library)
NC_FFLAGS := $(shell nf-config --fflags)
NC_FLIBS   := $(shell nf-config --flibs) $(shell nc-config --libs)

# Compiler flags
FFLAGS = -Wall -O2 $(NC_FFLAGS)

# Directories
SRC_DIR   = src/fortran
BUILD_DIR = build

# Source files
MODEL_SRC = $(SRC_DIR)/tasep.f90
SIM_SRC   = $(SRC_DIR)/simulation.f90
IO_SRC    = $(SRC_DIR)/io.f90
TEST_SRC  = $(SRC_DIR)/test_simulation.f90
FD_SRC    = $(SRC_DIR)/fundamental_diagram.f90

# Executables
TEST_EXE = $(BUILD_DIR)/test_simulation
FUND_EXE = $(BUILD_DIR)/fundamental_diagram

# Default simulation parameters (override: make run L=50 ALPHA=0.3 BETA=0.7)
L       ?= 10
N_STEPS ?= 100
ALPHA   ?= 0.5
BETA    ?= 0.5

# Default target: build both executables
all: $(TEST_EXE) $(FUND_EXE)

# Build simulation driver (needs NetCDF)
$(TEST_EXE): $(MODEL_SRC) $(SIM_SRC) $(IO_SRC) $(TEST_SRC)
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) $^ -o $@ $(NC_FLIBS)

# Build fundamental diagram sweep
$(FUND_EXE): $(MODEL_SRC) $(SIM_SRC) $(IO_SRC) $(FD_SRC)
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) $^ -o $@ $(NC_FLIBS)

# Run simulation (params can be overridden: make run L=50 N_STEPS=500 ALPHA=0.3 BETA=0.7)
run: $(TEST_EXE)
	./$(TEST_EXE) $(L) $(N_STEPS) $(ALPHA) $(BETA)

# Run fundamental diagram sweep
run-fd: $(FUND_EXE)
	./$(FUND_EXE)

# Clean build files
clean:
	rm -rf $(BUILD_DIR)

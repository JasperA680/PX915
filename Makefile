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

# Source files
MODEL_SRC = $(SRC_DIR)/tasep.f90
SIM_SRC   = $(SRC_DIR)/simulation.f90
IO_SRC    = $(SRC_DIR)/io.f90
TEST_SRC  = $(SRC_DIR)/test_simulation.f90

# Executable
TEST_EXE = $(BUILD_DIR)/test_simulation

# Default target
all: test

# Build test program
test: $(TEST_EXE)

$(TEST_EXE): $(MODEL_SRC) $(SIM_SRC) $(IO_SRC) $(TEST_SRC)
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) $^ -o $@ $(NC_FLIBS)


# Run the test
run: test
	./$(TEST_EXE)

# Clean build files
clean:
	rm -rf $(BUILD_DIR)

# Compiler 
FC = gfortran

# Compiler flags
FFLAGS = -Wall -O2

# Directories 
SRC_DIR = src/fortran
BUILD_DIR = build

# Source files
MODEL_SRC = $(SRC_DIR)/tasep.f90
SIM_SRC = $(SRC_DIR)/simulation.f90
TEST_SRC = $(SRC_DIR)/test_simulation.f90

# Executable 
TEST_EXE = $(BUILD_DIR)/test_simulation

# Default target
all: test 

# Build test program 
test: $(TEST_EXE)

$(TEST_EXE): $(MODEL_SRC) $(SIM_SRC) $(TEST_SRC)
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) $^ -o $@


# Run the test
run: test
	./$(TEST_EXE)

# Clean build files
clean:
	rm -rf $(BUILD_DIR)

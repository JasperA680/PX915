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
BENCH_SRC = $(SRC_DIR)/benchmark_tasep.f90

# PDE source files (pde_flux and pde_lanechange must precede pde_module)
PDE_FLUX_SRC = $(SRC_DIR)/pde_flux.f90
PDE_LC_SRC   = $(SRC_DIR)/pde_lanechange.f90
PDE_MOD_SRC  = $(SRC_DIR)/pde_module.f90
PDE_DRV_SRC  = $(SRC_DIR)/pde_driver.f90

# Executables
TEST_EXE  = $(BUILD_DIR)/test_simulation
FUND_EXE  = $(BUILD_DIR)/fundamental_diagram
PDE_EXE   = $(BUILD_DIR)/pde_solver
BENCH_EXE = $(BUILD_DIR)/benchmark_tasep

# Default simulation parameters (override: make run L=50 ALPHA=0.3 BETA=0.7)
L       ?= 10
N_STEPS ?= 100
ALPHA   ?= 0.5
BETA    ?= 0.5

# Default PDE parameters (override: make run-pde PDE_M=400 PDE_IC=riemann)
PDE_M      ?= 200
PDE_STEPS  ?= 500
PDE_VMAX   ?= 1.0
PDE_RHOMAX ?= 1.0
PDE_LEFT   ?= 0.1
PDE_RIGHT  ?= 0.9
PDE_IC     ?= riemann
PDE_FLUX   ?= lf
PDE_BC     ?= open
PDE_OUT    ?= data/output/pde_simulation.nc

# Default target: build all executables
all: $(TEST_EXE) $(FUND_EXE) $(PDE_EXE) $(BENCH_EXE)

# Build simulation driver (needs NetCDF)
$(TEST_EXE): $(MODEL_SRC) $(SIM_SRC) $(IO_SRC) $(TEST_SRC)
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) $^ -o $@ $(NC_FLIBS)

# Build fundamental diagram sweep
$(FUND_EXE): $(MODEL_SRC) $(SIM_SRC) $(IO_SRC) $(FD_SRC)
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) $^ -o $@ $(NC_FLIBS)

# Build PDE solver (pde_flux + pde_lanechange before pde_module due to USE deps)
$(PDE_EXE): $(PDE_FLUX_SRC) $(PDE_LC_SRC) $(PDE_MOD_SRC) $(PDE_DRV_SRC)
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) $^ -o $@ $(NC_FLIBS)

pde: $(PDE_EXE)

# Build benchmark driver (serial baseline for the parallelisation effort)
$(BENCH_EXE): $(MODEL_SRC) $(SIM_SRC) $(IO_SRC) $(BENCH_SRC)
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) $^ -o $@ $(NC_FLIBS)

# Run benchmark sweep (writes data/output/benchmark.nc)
benchmark: $(BENCH_EXE)
	mkdir -p data/output
	./$(BENCH_EXE)

# Run simulation (params can be overridden: make run L=50 N_STEPS=500 ALPHA=0.3 BETA=0.7)
run: $(TEST_EXE)
	./$(TEST_EXE) $(L) $(N_STEPS) $(ALPHA) $(BETA)

# Run fundamental diagram sweep
run-fd: $(FUND_EXE)
	./$(FUND_EXE)

# Run PDE solver (params can be overridden: make run-pde PDE_M=400 PDE_IC=riemann)
run-pde: $(PDE_EXE)
	./$(PDE_EXE) $(PDE_M) $(PDE_STEPS) $(PDE_VMAX) $(PDE_RHOMAX) $(PDE_LEFT) $(PDE_RIGHT) $(PDE_IC) $(PDE_FLUX) $(PDE_BC) $(PDE_OUT)

# Clean build files
clean:
	rm -rf $(BUILD_DIR)

# Compiler
FC = gfortran

# NetCDF flags (Fortran wrapper + underlying C library).
# nf-config --flibs and nc-config --libs both append -lnetcdf, which double-links
# on this Mac.  Strip the trailing -lnetcdf from nf-config and let nc-config provide it.
NC_FFLAGS := $(shell nf-config --fflags)
NC_FLIBS   := $(shell nf-config --flibs | sed -E 's/-lnetcdf$$//') $(shell nc-config --libs)

# Compiler flags
FFLAGS = -Wall -O2 $(NC_FFLAGS)

# Directories
SRC_DIR   = src/fortran
BUILD_DIR = build


# Network-model sources
VEHICLE_SRC      = $(SRC_DIR)/vehicle.f90
NETWORK_SRC      = $(SRC_DIR)/road_network.f90
NETWORK_INIT_SRC = $(SRC_DIR)/network_init.f90
NS_MODEL_SRC     = $(SRC_DIR)/NS_model.f90
JUNCTION_SRC     = $(SRC_DIR)/junction.f90
LANE_CHANGE_SRC  = $(SRC_DIR)/lane_change.f90
NET_SIM_SRC      = $(SRC_DIR)/network_simulation.f90

# Network builder + JSON config + NetCDF writer + driver (Python frontend integration)
NETWORK_BUILDER_SRC = $(SRC_DIR)/network_builder.f90
JSON_CONFIG_SRC     = $(SRC_DIR)/json_config.f90
NETWORK_IO_SRC      = $(SRC_DIR)/network_io.f90
RUN_NETWORK_SRC     = $(SRC_DIR)/run_network.f90

NETWORK_LIB_SRC = $(VEHICLE_SRC) $(NETWORK_SRC) $(NETWORK_INIT_SRC) \
                  $(LANE_CHANGE_SRC) $(JUNCTION_SRC) $(NS_MODEL_SRC) $(NET_SIM_SRC)

# PDE source files (pde_flux and pde_lanechange must precede pde_module due to USE deps)
PDE_FLUX_SRC = $(SRC_DIR)/pde_flux.f90
PDE_LC_SRC   = $(SRC_DIR)/pde_lanechange.f90
PDE_MOD_SRC  = $(SRC_DIR)/pde_module.f90
PDE_DRV_SRC  = $(SRC_DIR)/pde_driver.f90

# Executables
PDE_EXE               = $(BUILD_DIR)/pde_solver
RUN_NETWORK_EXE       = $(BUILD_DIR)/run_network

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
PDE_VLIMIT ?= 1.0
PDE_OUT    ?= data/output/pde_simulation.nc

# Default target: build all executables
all: run_network $(PDE_EXE)

pde: $(PDE_EXE)

# Run targets

run-pde: $(PDE_EXE)
	./$(PDE_EXE) $(PDE_M) $(PDE_STEPS) $(PDE_VMAX) $(PDE_RHOMAX) $(PDE_LEFT) $(PDE_RIGHT) $(PDE_IC) $(PDE_FLUX) $(PDE_BC) $(PDE_VLIMIT) $(PDE_OUT)

# Python-frontend driver: reads JSON config, runs sim, writes NetCDF.
run_network: $(RUN_NETWORK_EXE)

$(RUN_NETWORK_EXE): $(NETWORK_LIB_SRC) $(NETWORK_BUILDER_SRC) $(JSON_CONFIG_SRC) $(NETWORK_IO_SRC) $(RUN_NETWORK_SRC)
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) $^ -o $@ $(NC_FLIBS)

# Clean build files
clean:
	rm -rf $(BUILD_DIR)
	rm -f $(SRC_DIR)/*.mod *.mod

EXE_SINGLE   = conway_single
EXE_MPI_CART = conway_2d_mesh
EXE_MPI_ROW  = conway_1d_mesh
EXE = $(EXE_SINGLE) $(EXE_MPI_CART) $(EXE_MPI_ROW)

CXX = g++
MPICXX = mpicxx

OBJ = main.o   cell_board.o
CPP = main.cpp cell_board.cpp
HDR = cell_board.h  aux.h

OPT = -O2 -g

OPENMPI=BUILD_WITH_OPENMPI
CART=CARTESIAN_PARTITION
ROW=ROWWISE_PARTITION


all : $(EXE_SINGLE) $(EXE_MPI_CART) # $(EXE)

# ============================
#  Single-Process without MPI
# ============================

$(EXE_SINGLE) : $(CPP)
	$(CXX) $(OPT) -o $@ $^

$(EXE_MPI_CART) : $(CPP)
	$(MPICXX) $(OPT) -D$(OPENMPI) -D$(CART) -o $@ $^

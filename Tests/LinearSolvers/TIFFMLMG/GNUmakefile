AMREX_HOME ?= ../../../

#DEBUG	= FALSE
DEBUG	= TRUE

# DIM	= 3
DIM	= 2

# COMP    = gnu

USE_MPI   = TRUE
#USE_OMP   = TRUE
USE_OMP   = FALSE

LIBRARIES += -L/Users/hans/anaconda3/pkgs/libtiff-4.0.8-h8cd0352_9/lib/ -ltiff

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/MLMG/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/C_CellMG/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

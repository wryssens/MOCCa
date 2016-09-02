#-------------------------------------------------------------------------------
# Definition of source files, object location, module location etc...
#-------------------------------------------------------------------------------
OBJDIR :=   obj
SRCDIR :=   src
MODDIR :=   mod

TARGET :=   MOCCa.exe
SRC    :=   CompilationInfo.f90 GenInfo.f90 Force.f90 OptimizedDerivatives.f90 Derivatives.f90 
SRC    +=   CoulombDerivatives.f90 Mesh.f90 Spinor.f90 
SRC    +=   Spwf.f90 SpwfStorage.f90 Damping.f90 Densities.f90 
SRC    +=   Moments.f90 SpecialMoments.f90 MultiGrid.f90 Coulomb.f90 PairingInteraction.f90 
SRC    +=   LipkinNogami.f90 HartreeFock.f90 HFB.f90 GradientHFB.f90 BCS.f90 Pairing.f90 Cranking.f90 
SRC    +=   MeanFields.f90 Energy.f90 ImaginaryTime.f90  DensityMixing.f90 
SRC    +=   Transform.f90 SpwfFactory.f90 Interfaces.f90 InOut.f90 Test.f90 Main.f90 
LIBS   :=   -llapack -lblas

#Make the lists of objects
OBJ :=      $(patsubst %.f90,$(OBJDIR)/%.o,$(SRC))
#-------------------------------------------------------------------------------
# Compilers and some recommended options

#Default compiler is gfortran
CXX :=      gfortran-5
# Default behaviour is not debugging
DEBUG := no

ifeq ($(CXX),gfortran-5)

  ifeq ($(DEBUG),no)
    #Optimal flag
    CXXFLAGS := -O3 
  else
    # Debugging flag
    CXXFLAGS= -Og -fbacktrace -fcheck=all 
    # Alternative
    #CXXFLAGS= -Wno-tabs -fbacktrace -fcheck=all -ggdb -Wall
  endif
  #Move the .mod files to their own directory
  CXXFLAGS += -J$(MODDIR)

else ifeq ($(CXX),pgf90)
  ifeq ($(DEBUG),0)
    #Optimal flag
		CXXFLAGS := -O3
    #	else
		# Debugging flag
		CXXFLAGS=-traceback -Mbounds
  endif
	#Move the .mod files to their own directory
	CXXFLAGS += -module $(MODDIR)
else ifeq ($(CXX),ifort)
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # The `-assume realloc-lhs' is mandatory for correct ifort behaviour.
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ifeq ($(DEBUG),no)
		#Optimal flag
		CXXFLAGS=-O3 -assume realloc-lhs
  else
		# Debugging flag
		CXXFLAGS=-g -traceback -check bounds -warn all -assume realloc-lhs
  endif
	#Move the .mod files to their own directory
	CXXFLAGS += -module $(MODDIR)
endif


#-------------------------------------------------------------------------------

.PHONY: all clean

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

clean:
	rm  -f $(OBJDIR)/*
	rm  -f $(MODDIR)/*

$(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(CXX) $(CXXFLAGS) -c  $< -o $@ 

#-------------------------------------------------------------------------------
# Some observations on my experiences of compiling MOCCa.
#
# 
#    Versions of compilers I've verified that DON'T WORK:
#   COMPILER                  SYSTEM           REASON
#   ----------------------------------------------------------------------------
#  ifort13.0                 HERCULES        COMPILES, BUT SEGFAULTS IN COULOMB
#                                             CHANCES ARE THIS CAN BE CORRECTED
#                                             WITH SOME EFFORT
#  ifort12.1 (update 4)      HYDRA           COMPILES, BUT SEGFAULTS IN COULOMB
#                                             CHANCES ARE THIS CAN BE CORRECTED
#                                             WITH SOME EFFORT
#  gfortran 4.6.1            HYDRA           COMPILES WITH MINOR MODIFICATIONS
#                                             BUT SEGFAULTS HORRIBLY
#  gfortran 4.4.1            HYDRA           CANNOT COMPILE
#  pgf90 13.10               HERCULES        compiles, but throws a (so far)
#                                             untraceable segfault
#  pgf90 11.3                HERCULES        does not compile
#
#     Versions of compilers I've verified that DO WORK:
#   COMPILER                  SYSTEM          REMARKS
#   ----------------------------------------------------------------------------
# - gfortran 4.8.2            HYDRA   
# - gfortran 4.8.1            HYDRA 
# - gfortran 4.7.2            HERCULES 
# - ifort15.01                HYDRA
# - ifort14.02                HERCULES
# - pgf90 11.4-0              HYDRA           Works now, after some effort
#                                             regarding abstract interfaces
#                                             for procedure pointers
#                                             However does not work with 
#                                             optimisation enabled (see below)
#                                                   
#
# This list is not comprehensive, on HERCULES there are more versions to test...
#-------------------------------------------------------------------------------
# Some other notes of interest:
# *) The order of the .f90 files is important, it tells the compiler which 
#    modules depend on one another.
# *) MOCCa needs (almost full) FORTRAN 2003 support, and this is mostly visible 
#    in its use of procedure pointers and polymorphic variables.
#    This means that you should compile with recent compilers whereever possible.
# *) The ifort compiler needs the 
#      -assume realloc-lhs
#    flag. The idea is that the ifort compiler does not automatically 
#    allocate memory correctly upon intrinsic assignment. With this flag,
#    it does enforce this. I do not understand why this is ifort defaults 
#    behaviour, intrinsic allocation is in the FORTRAN standard as far as
#    I know.
# *) When trying to debug with gdb some things to keep in mind:
#	    -  allocatable arrays are difficult to print with it, google to 
#          see the syntax
#       -  Doublecheck your gfortran & gdb versions, since these typically
#          do not like each other when one or the other is older. 
# *) Best results (timewise) I obtain with the following compilation options
#       - gfortran:
#          -Ofast -funsafe-loop-optimizations  
#       - ifort:
#           No testing done yet
#       - PGI:
#           No testing done yet
# *) Do not try to read MOCCa files produced by PGI-compiled MOCCa with non-PGI
#    MOCCa exes. This does not detect symmetries on the file for some reason 
#    I've yet to find out.
# *) It was pretty hard to make MOCCa comply with the PGI compilers. They are 
#    the only compiler that I ever used that insists on having defined an 
#    interface for a procedure pointer BEFORE the declaration of the pointer.
#    Thus it is impossible to use procedure pointers in a module which derive
#    from a procedure in the same module. For this reason, you will sometimes
#    find abstract interfaces in some modules.(Notably derivatives module...)
# *) In addition to this, compiling MOCCa with PGI compilers with any kind of
#    optimisation results in segfaults at runtime. To be investigated
#-------------------------------------------------------------------------------

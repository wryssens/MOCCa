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

CC=gfortran
#CC=ifort
#CC=pgf90

#--------------------------------------------------------------------------------
# Set Compiler flags
ifeq ($(CC),gfortran)

#  DEBUGGING FLAG
#   CFLAGS= -Wno-tabs -fbacktrace -fcheck=all -ggdb -Wall
#   CFLAGS= -Og -fbacktrace -fcheck=all 

#  OPTIMAL FLAG
   CFLAGS= -Ofast

else ifeq ($(CC),pgf90)
  
#  DEBUGGING FLAG
#  CFLAGS=-traceback -Mbounds

#  OPTIMAL FLAG
   CFLAGS=-O3

else ifeq ($(CC),ifort)
# The `-assume realloc-lhs' is mandatory for correct ifort behaviour.
# DEBUGGING FLAG
# CFLAGS=-g -traceback -check bounds -warn all -assume realloc-lhs

# OPTIMAL FLAG
  CFLAGS=-O3 -assume realloc-lhs
endif
#--------------------------------------------------------------------------------
# Library linking
LIB=-lopenblas

#Source files
SRC=CompilationInfo.f90 GenInfo.f90 Force.f90 Derivatives.f90  GnuForInterface.f90  CoulombDerivatives.f90 Mesh.f90 Spinor.f90 Spwf.f90 SpwfStorage.f90 Damping.f90 Densities.f90  Moments.f90  MultiGrid.f90 Coulomb.f90 PairingInteraction.f90 LipkinNogami.f90 HFB.f90 BCS.f90 Pairing.f90 Cranking.f90 MeanFields.f90 ImaginaryTime.f90 Energy.f90 DensityMixing.f90 InOut.f90 Test.f90 Main.f90 SpwfFactory.f90
#Modules to build
OBJ=CompilationInfo.o GenInfo.o Force.o Derivatives.o  GnuForInterface.o  CoulombDerivatives.o  Mesh.o Spinor.o Spwf.o SpwfStorage.o Damping.o Densities.o  Moments.o  MultiGrid.o Coulomb.o PairingInteraction.o LipkinNogami.o HFB.o BCS.o Pairing.o Cranking.o MeanFields.o ImaginaryTime.o Energy.o DensityMixing.o Transform.o SpwfFactory.o InOut.o Test.o Main.o
#Modules to build for the Clusters program
CLUOBJ=CompilationInfo.o GenInfo.o Force.o GnuForInterface.o Derivatives.o CoulombDerivatives.o  Mesh.o Spinor.o Spwf.o SpwfStorage.o Densities.o Moments.o  MultiGrid.o Coulomb.o PairingInteraction.o HFB.o BCS.o Pairing.o Cranking.o MeanFields.o Energy.o InOut.o Transform.o Clusters.o

%.o: %.f90
	$(CC) $(CFLAGS) -c $< $(OPT) $(LIB)

MOCCa.exe: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIB)
	
Clusters.exe : $(CLUOBJ)
	$(CC) $(CFLAGS) -o $@ $^

clean :
	rm  *.mod
	rm  *.o

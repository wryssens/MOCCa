
# MOCCa
> Copyright W. Ryssens, P.H. Heenen & M. Bender
> Tuesday, 16. June 2020 10:10AM 

-------
### File structure

Files in the repository's main folder are

*  `Makefile`: Makefile to compile MOCCa. 

Useful folders/files are

*  `src/` :   directory of the compileable source code.
* `scripts/` : many imperfect scripts by W.R. for various analysis and plotting purposes. 
* `forces.param`: contains a large set of predefined functional values

Auxiliary folders are

* `mod/` : Auxiliary folder where the Makefile places the compiled `.mod` files.
* `obj/` :  Auxiliary folder fo the compiled `.obj` files. 

-------
### How do I compile the code?

The easiest way is to type

`make`

Currently defined options to pass into the Makefile are

* `CXX`: your preferred choice of compiler. Default = `gfortran`.
* `CXXFLAGS`:  options for your preferred compiler.  


__Things do check if the code does not want to compile__

* Are the `src/` and `mod/`  directories present in your directory tree?

-------
### How do I run MOCCa? 

Simply do

	./MOCCa.exe < input.in 
	
and make sure that either you tell MOCCa to initialize from scratch or read a wavefunction file that is in the current directory.

An example input file for a Mg24 calculation from scratch is:

	cat << EOF > MOCCa.data
	&geninfo
	neutrons=12, protons =12
	nx=10 , ny=10 , nz=10
	dx=0.8
	MaxIter=100
	ParameterEstimation=1 
	! Give the code licence to select its own evolution parameters
	/
	&SpwfStorage
	nwt=40
	! 40 total spwfs
	/
	&Densit
	MixingScheme=2
	! Selection of density mixing scheme
	/
	&Pairing
	Type="HFB" , PairingNeutron=1250, AlphaNeutron=1.0
	! Type of pairing and parameters of the pairing interaction
	GuessKappa=.true.
	!Note only put this on to start calculations: it initializes kappa to a
	!guessed value
	/
	&HFConfig
	/
	&Blocking
	/
	&Derivatives
	MaxFDOrder    = -1 , MaxFDLapOrder = -1 , BStack=.false.
	! Selects the kind of derivatives to use on the mesh; maximum precision
	/
	&MomentParam
	MoreConstraints=.true.
	/
	&MomentConstraint
	Constraint=50, l=2, m=0, iteration=10
	! A quadrupole constraint for 10 iterations to kick things off
	/
	&
	&InAndOutput
	OutputFileName="MOCCa.out"
	InputFileName="INIT"
	nwninit=20
	nwpinit=20
	! Initialize the code from scratch (not from a .wf file), take 20 proton
	! and 20 neutron wavefunctions. (Needs to sum to NWT)
	/
	&Force
	afor="SLy4"
	! Name of the functional, read from forces.param
	/
	&Coulomb
	/
	&Cranking
	/
	EOF
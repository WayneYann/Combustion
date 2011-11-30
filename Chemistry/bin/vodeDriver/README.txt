
Compact App: vodeDriver

Author: Marc Day, CCSE-LBNL  (MSDay@lbl.gov, 510-486-5076)

Summary: Evolves ODEs associated with chemical kinetics

Additional info:

The driver reads in a set of cells containing state data, and calls
the function ChemDriver::solveTransient over the set.
ChemDriver::solveTransient is called in LMC to perform the ODE
integration substep in the Strang-split version of the low Mach number
reacting flow integration (Day and Bell, Combust. Theory Modelling
4(4) pp.535-556, 2000).  The input data at each cell consists of
temperature (in K), the mass fractions of all chemical species in the
model, the pressure (in atmospheres) and the time interval.  Data is
passed to the solver over a box of cells (in a FArrayBox, or "FAB),
and the pressure and time interval are the same for all cells.  In
LMC, the state data is stored in an array of FABs (a MultiFab) - the
data is distributed in parallel by FAB, and an "owner computes" rule
applies.  In this demo, a single FAB containing approximately 1000
cells is initialized using point values nonuniformly sampled through a
steady premixed flame solution, which was computed using the CHEMKIN
code PREMIX.  This type of intialization provides a realistic sample
of states expected in production LMC runs.

Background and Operational comments:

(1) The chemical kinetics and thermodynamic relationships are stored
internally in a parameterized format consistent with the CHEMKIN
specification.  The CHEMKIN software library is based on a set of
Fortran routines that read a textfile database describing the species,
thermodynamics and kinetics, and store this information in Fortran
common blocks.  In order to realize substantial performance gains,
CCSE has converted a subset of the CHEMKIN access/evaluation routines
into C/C++ code.  As a result, while CHEMKIN uses a runtime process to
setup a database, CCSE uses a complie-time system.  Therefore, in 
order to change models, the executable must be relinked.

(2) We have provided two example cases for analysis.  One is based on
"CHEMH", and implements a 9-species, 27-reaction model for hydrogen
combustion.  The other is based on "DRM19", and implements a
21-species, 84-reaction model for methane combustion.  Initial data is
specific to each of the two options, and is provided for the tests
using binary FAB files included in this distribution.

(3) In the GNUmakefile, the model is selected by setting the variable
CHEMISTRY_MODEL to DRM19 or CHEMH.  Based on the value of this
variable, a chemistry-specific C++ file is included in the build to
implement the relevant model.  After compilation, the demo is run by
providing the appropriate initial data file at the command line, viz.

       vodeDriver...ex pmf_file=chem-H_0370.fab, or
       vodeDriver...ex pmf_file=drm19_0700.fab

In the driver routine, a cursory check is performed to ensure that a
compatible file is indicated.  A file name must be provided.  Optionally,
the user may set the pressure, time interval, or verbosity.

(4) ChemDriver::solveTransient passes the input data to the Fortran
routine FORT_CONPSOLV (in ChemDriver_{2,3}D.F, depending on the value
of DIM in GNUmakefile at compile time), which loops over the cells and
for each one calls DVODE to integrate a system from t=0 to t=dt.
FORT_CONPSOLV optionally computes a number of diagnostics not relevant
to this demo.  The important piece here is the call to DVODE, which
includes in its arguments the name of a function that returns a vector
of dY/dt given a vector of Y.  For our purposes, the function passed
is "conpFY" (in ChemDriver_F.F).  DVODE also requires a Jacobian of
the discrete equations for time integration, but this is computed by
DVODE internally using a finite difference approach.

(5) Input to conpFY includes the mass fractions, Ym (m=1,Nspec) and
temperature, T.  Output includes the instantateous rate of change of
these quantities, following:

            dYm/dt = Rm . Wm / rho
            dT/dt  = -sum( Hm . dYm/dt ) / Cpmix

Here, chemical Rm is the production rate (moles/cm^3.s), W is the
molecular weight and Hm is the enthalpy (including heat of formation)
of species m.  Cpmix is the specific heat of the mixture.  rho is the
mass density, computed with the current values of Ym and T, and the
(constant) pressure, by the routine CKRHOY.  Cpmix is computed by the
routine CKCPBS, Rm is computed by CKWC, and Hm is computed by CKHMS.
These routines are in chem-H.cpp, or drm19.cpp.

(6) The calls to DVODE are such that each cell is integrated
independently if the compile-time preprocessor variable ALWAYS_NEW_J
is set.  If this variable is unset, DVODE attempts the first
integration step of a new cell using the Jacobian matrix from the
previous cell (if the solver fails to find a solution for that step, a
new Jacobian is automatically computed).  Note that the tolerances are
set such that the final results are essentially independent of this
choice.  If ALWAYS_NEW_J is set, the results should be exactly
independent of the order in which the cells are updated, otherwise
there will be a weak dependency, both in the effort to integrate the
cells and in the final result.



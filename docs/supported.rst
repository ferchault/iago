Features
========

Supported Quantities
--------------------

If present in the output of a calculation, the following quantities are recognised by *iago* for cp2k QS. All quantities here
consist of raw data as calculated by the simulation package itself. None of them is derived (except for unit
conversions).

================== =============== ===========================================================
Column             Unit            Description
================== =============== ===========================================================
frame              --              Frame number. Integer.
steptime           fs              Time of the frame. Float.
a, b, c            Angstrom        Cell lengths. Float.
alpha, beta, gamma Degree          Cell angles. Float.
temperature        Kelvin          Instantaneous temperature. Float.
pressure           bar             Instantaneous pressure. Float.
conserved          eV              Conserved quantity. Float.
coreself           eV              Core Self Energy. Float.
corehamiltonian    eV              Core Hamiltonian. Float.
hartree            eV              Hartree Energy. Float.
xc                 eV              Exchange-Correlation Energy. Float.
hfx                eV              Hartree-Fock-Exchange Energy. Float.
dispersion         eV              Dispersion energy. Float.
total              eV              Total energy. Float.
potential          eV              Potential energy. Float.
kinetic            eV              Kinetic energy. Float.
drift                              Temperature drift. Float.
iasd                               Integrated absolute spin density. Float.
s2                                 Single determinant S**2. Float.
scfcycles          --              Number of SCF steps. Integer.
otcycles           --              Number of outer SCF steps. Integer.
cputime            s               CPU time per MD step. Float.
================== =============== ===========================================================

If present in the output of a calculation, the following quantities are recognised by *iago* for NAMD. All quantities here consist of raw data as calculated by the simulation package itself. None of them is derived (except for unit
conversions). Any keyword starts with NAMD is either an averaged output from NAMD or not clearly defined in NAMD MANUAL.

================== =============== ===========================================================
Column             Unit            Description
================== =============== ===========================================================
frame              --              Frame number. Integer.
conserved          eV              Conserved quantity. Float.
temperature        Kelvin          Instantaneous temperature. Float.
pressure           bar             Instantaneous pressure. Float.
volume		   Angstrom^3	   Instantaneous volume. float.
bond               eV              bond energy. Float.
angle              eV              angle energy. Float.
dihedral           eV		   dihedral energy. Float.
improper	   eV		   improper energy. Float.
electrostatic	   eV		   electrostatic energy. Float.
vdw		   eV		   van der Waals energy. Float.
boundary	   eV		   boundary energy from spherical boundary conditions and harmonic restraints. Float.
external	   eV		   MISC term in NAMD. External electric fields and various steering forces. Float.
kinetic_cm	   eV		   Classical kinetic energy. Float.
total_cm	   eV		   Classical total energy. Float.
potential_cm	   eV		   classical potential energy. Float.
NAMD_total3	   eV		   NAMD TOTAL3. Float.
NAMD_GPRESSURE	   eV		   NAMD GPRESSURE. Float.
NAMD_GPRESSUREAVG  eV		   NAMD GPRESSUREAVG. Float.
NAMD_PRESSUREAVG   eV		   NAMD PRESSUREAVG. Float.
NAMD_TEMPAVG	   Kelvin	   NAMD TEMPAVG. Float.
================== =============== ===========================================================

.. _supported-calculation:

Supported Calculations
----------------------

Based on the trajectory data and the simulation package log file, the following quantities can be calculated.

.. py:currentmodule:: iago.Analyser

.. autosummary::

	Analyser.dynamic_plane
	Analyser.dynamic_distance


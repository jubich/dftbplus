#------------------------------------------------------------------------------#
#  DFTB+: general package for performing fast atomistic simulations            #
#  Copyright (C) 2006 - 2022  DFTB+ developers group                           #
#                                                                              #
#  See the LICENSE file for terms of usage and distribution.                   #
#------------------------------------------------------------------------------#


'''
Interface module for the communication between DFTB+
and Python via cffi.
'''

import os
from .dftb_cffi import ffi, lib
import numpy as np


# DFTB+ conversion factors
# (according to src/dftbp/common/constants.F90)
HARTREE__EV = 27.2113845
BOHR__AA = 0.529177249


class DftbPlus:
    '''cffi interface to DFTB+.

    A shared library of DFTB+ is loaded and calculations of several
    physical quantities can be carried out. After completing the calculations,
    the results of the specified properties can be extracted easily.

    '''


    def __init__(self, hsdpath='./dftb_in.hsd', logfile=None, mpi=None):
        '''Initializes a cffi DFTB+ calculator object.

        Args:

            hsdpath (str): path to DFTB+ input file
            logfile (str): name of log file
            mpi (int): MPI-communicator id

        '''
        self._dftbpluslib = lib

        input_str = ffi.new("char[]", str.encode(hsdpath))

        # pointer to DFTB+ instance
        self._dftb_handler = ffi.new("DftbPlus *")
        # pointer to DFTB+ input
        self._dftb_input = ffi.new("DftbPlusInput *")

        self._refobj_usr = None
        self._calc_extpot_c = None
        self._calc_extpotgrad_c = None
        self._calc_extpot_usr = None
        self._calc_extpotgrad_usr = None

        if mpi is None:
            if logfile is not None:
                output_str = ffi.new("char[]", str.encode(logfile))
                self._dftbpluslib.dftbp_init(self._dftb_handler, output_str)
            else:
                self._dftbpluslib.dftbp_init(self._dftb_handler, ffi.NULL)
        else:
            mpi = ffi.new("int[]", mpi)
            if logfile is not None:
                output_str = ffi.new("char[]", str.encode(logfile))
                self._dftbpluslib.dftbp_init_mpi(self._dftb_handler, output_str, mpi)
            else:
                self._dftbpluslib.dftbp_init_mpi(self._dftb_handler, ffi.NULL, mpi)

        self._dftbpluslib.dftbp_get_input_from_file(self._dftb_handler,
                                                    input_str, self._dftb_input)

        self._dftbpluslib.dftbp_process_input(self._dftb_handler,
                                              self._dftb_input)

        self._natoms = self._dftbpluslib.dftbp_get_nr_atoms(self._dftb_handler)


    def set_geometry(self, coords, latvecs=None, origin=None):
        '''Sets up the desired geometry.

        Args:

            coords (2darray): absolute atomic positions (in atomic units)
            latvecs (2darray): lattice vectors or None for non-periodic
                structures (in atomic units)
            origin (2darray): Coordinate origin (in atomic units)

        '''
        periodic = latvecs is not None
        coords = ffi.from_buffer(
            "double *", np.ascontiguousarray(coords, dtype=np.float64))
        if periodic:
            latvecs = ffi.from_buffer(
                "double *", np.ascontiguousarray(latvecs, dtype=np.float64))
            if origin is None:
                self._dftbpluslib.dftbp_set_coords_and_lattice_vecs(
                    self._dftb_handler, coords, latvecs)
            else:
                origin = ffi.from_buffer(
                    "double *", np.ascontiguousarray(origin, dtype=np.float64))
                self._dftbpluslib.dftbp_set_coords_lattice_origin(
                    self._dftb_handler, coords, latvecs, origin)
        else:
            self._dftbpluslib.dftbp_set_coords(self._dftb_handler, coords)


    def set_external_potential(self, extpot, extpotgrad=None):
        '''Sets up an external potential.

        Args:

            extpot (2darray): external potential at the position of each atom
                (in atomic units). Shape: [natom].
            extpotgrad (2darray): gradient of the external potential at each
                atom (in atomic units). Shape: [natom, 3]. This parameter is
                optional, you can pass None if you did not ask DFTB+ to
                calculate forces.

        '''
        if extpotgrad is not None:
            extpotgrad = ffi.from_buffer(
                "double *", np.ascontiguousarray(extpotgrad, dtype=np.float64))
        else:
            extpotgrad = ffi.NULL
        extpot = ffi.from_buffer(
            "double *", np.ascontiguousarray(extpot, dtype=np.float64))
        self._dftbpluslib.dftbp_set_external_potential(
            self._dftb_handler, extpot, extpotgrad)


    def register_ext_pot_generator(self, refobj, calc_extpot, calc_extpotgrad):
        '''Registers callback functions for population dependent external
            potential calculations.

        Args:

            refobj (pointer): user defined data struct or class which contains
                the necessary data for the potential calculation
            calc_extpot (pointer): pointer to user defined callback function
                which DFTB+ should call, whenever the population dependent
                external potential should be calculated
            calc_extpotgrad (pointer): pointer to user defined callback function
                which DFTB+ should call, whenever the gradient of the population
                dependent external potential should be calculated

        '''
        self._refobj_usr = refobj
        self._self_handle = ffi.new_handle(self)
        self._calc_extpot_usr = calc_extpot
        self._calc_extpotgrad_usr = calc_extpotgrad

        self._calc_extpot_c = self._dftbpluslib._calc_extpot_callback
        self._calc_extpotgrad_c = self._dftbpluslib._calc_extpotgrad_callback

        self._dftbpluslib.dftbp_register_ext_pot_generator(
            self._dftb_handler, self._self_handle,
            self._calc_extpot_c, self._calc_extpotgrad_c)


    @ffi.def_extern()
    def _calc_extpot_callback(refobj, dqatom, extpotatom):
        '''Callback function wrapper to hide the necessary conversions of low
            level types into numpy arrays.

        Args:

            refobj (pointer): user defined data struct or class which contains
                the necessary data for the potential calculation
            dqatom (pointer): population difference with respect to reference
                population (usually the neutral atom). Note: population means
                electrons, so a positive number indicates electron excess
            extpotatom (pointer): potential at the position of each QM-atom.
                Note: it should be the potential as felt by an electron
                (negative potential value means attraction for an electron)

        '''
        self = ffi.from_handle(refobj)
        dqatom_size = ffi.sizeof(dqatom) * self._natoms
        dqatom_array = np.frombuffer(
            ffi.buffer(dqatom, dqatom_size)).reshape(self._natoms,)
        expot_size = ffi.sizeof(extpotatom) * self._natoms
        extpot_array = np.frombuffer(
            ffi.buffer(extpotatom, expot_size)).reshape(self._natoms,)

        self._calc_extpot_usr(self._refobj_usr, dqatom_array, extpot_array)


    @ffi.def_extern()
    def _calc_extpotgrad_callback(refobj, dqatom, extpotatomgrad):
        '''Callback function wrapper to hide the necessary conversions of low
            level types into numpy arrays.

        Args:

            refobj (pointer): user defined data struct or class which contains
                the necessary data for the potential calculation
            dqatom (pointer): population difference with respect to reference
                population (usually the neutral atom). Note: population means
                electrons, so a positive number indicates electron excess
            extpotatomgrad (pointer): potential gradient at the position of each
                QM-atom. Note: it should be the gradient of the potential as
                felt by an electron (negative potential value means attraction
                for an electron)

        '''
        self = ffi.from_handle(refobj)
        dqatom_size = ffi.sizeof(dqatom) * self._natoms
        dqatom_array = np.frombuffer(
            ffi.buffer(dqatom, dqatom_size)).reshape(self._natoms,)
        expotgrad_size = ffi.sizeof(extpotatomgrad) * self._natoms * 3
        extpotgrad_array = np.frombuffer(
            ffi.buffer(extpotatomgrad, expotgrad_size)).reshape(self._natoms, 3)

        self._calc_extpotgrad_usr(
            self._refobj_usr, dqatom_array, extpotgrad_array)


    def get_nr_atoms(self):
        '''Queries the number of atoms.

        Returns:

            self._natoms (int): number of atoms

        '''
        return self._natoms


    def get_energy(self):
        '''Performs the energy (Mermin free) calculation and queries the energy
            of the current geometry.

        Returns:

            energy[0] (float): calculated Mermin free energy (in atomic units)

        '''
        energy = ffi.new("double[1]")

        self._dftbpluslib.dftbp_get_energy(self._dftb_handler, energy)

        return energy[0]


    def get_gradients(self):
        '''Performs the calculation of the atomic gradients and queries the
            gradients of the current geometry.

        Returns:

            gradients (2darray): calculated gradients (in atomic units)

        '''
        gradients = np.ascontiguousarray(
            np.empty((self._natoms, 3)), dtype=np.float64)

        self._dftbpluslib.dftbp_get_gradients(
            self._dftb_handler, ffi.from_buffer("double *", gradients))

        return gradients


    def get_gross_charges(self):
        '''Queries the atomic Gross charges.

            Until state commit: 7c13358d (base: 19.1) of DFTB+:
            For non-trivial results, an energy or force calculation should be
            carried out beforehand. Otherwise, the atom populations are returned.
            More recent versions of DFTB+ adapt the behavior to the subroutines
            for the extraction of energy and gradients, therefore preventing the
            return of trivial charges without a calculation carried out.

            Sign convention: Electron has negative charge, so negative values
            indicate electron excess.

        Returns:

            grosschg (1darray): obtained Gross charges (in atomic units)

        '''
        grosschg = np.ascontiguousarray(
            np.empty(self._natoms), dtype=np.float64)

        self._dftbpluslib.dftbp_get_gross_charges(
            self._dftb_handler, ffi.from_buffer("double *", grosschg))

        return grosschg


    def get_cm5_charges(self):
        '''Queries the atomic CM5 charges.

            Sign convention: Electron has negative charge, so negative values
            indicate electron excess.

        Returns:

            charges (1darray): obtained CM5 charges (in atomic units)

        '''
        charges = np.ascontiguousarray(
            np.empty(self._natoms), dtype=np.float64)

        self._dftbpluslib.dftbp_get_cm5_charges(
            self._dftb_handler, ffi.from_buffer("double *", charges))

        return charges


    def get_api_version(self):
        """Returns current version of the DFTB+ API

        Returns:

            major[0] (int): major version
            minor[0] (int): minor version
            patch[0] (int): patch version

        """
        major = ffi.new("int[1]")
        minor = ffi.new("int[1]")
        patch = ffi.new("int[1]")
        self._dftbpluslib.dftbp_api(major, minor, patch)
        return major[0], minor[0], patch[0]


    def is_instance_safe(self):
        """Returns whether library was built with instance safe components
        only. Instance safe API support the creation of multiple concurrent
        DFTB+ instances within one process.

        Returns:

            (bool): Whether API is instance safe

        """
        return self._dftbpluslib.dftbp_is_instance_safe()


    def get_nr_kpoints(self):
        """Queries the nr. of k-points in the system.

        Returns:

            (int): Nr. of k-points

        """
        return self._dftbpluslib.dftbp_nr_kpoints(self._dftb_handler)


    def get_nr_orbitals(self):
        """Queries the number of basis functions for each atom in current
        geometry.

        Returns:

            (1darray): number of orbitals on each atom. Shape: [natom].

        """
        norbitals = ffi.new(f"int[{self._natoms}]")
        self._dftbpluslib.dftbp_get_nr_orbitals(self._dftb_handler, norbitals)
        orbitals = []
        for atom in range(self._natoms):
            orbitals.append(norbitals[atom])

        return np.asarray(orbitals)


    def get_masses(self):
        """Queries the masses for each atom in current geometry.

        Returns:

            (1darray): mass of each atom. Shape: [natom].

        """
        masses = ffi.new(f"double[{self._natoms}]")
        self._dftbpluslib.dftbp_get_masses(self._dftb_handler, masses)
        mass = []
        for atom in range(self._natoms):
            mass.append(masses[atom])

        return np.asarray(mass)


    def get_stress_tensor(self):
        """Queries the stress tensor of the current periodic box.

        Returns:

            stress_tensor (2darray): Stress Tensor for the periodic box.
                Shape: [3, 3]. Unit: Pascals.
        """
        stress_tensor = np.ascontiguousarray(np.empty((3, 3)), dtype=np.float64)

        self._dftbpluslib.dftbp_get_stress_tensor(
            self._dftb_handler, ffi.from_buffer("double *", stress_tensor))

        return stress_tensor


    def get_elstat_potential(self, locations):
        """Queries electrostatic potential in specific points.

        Args:

            locations (2darray): Coordinates of requested points.
                Shape: [nlocations, 3]. Unit: Bohr.

        Returns:

            (1darray): Values of electrostatic potential. Shape: [nLocations].
        """
        nlocations = ffi.new("int[]", locations.shape[0])
        locations = ffi.from_buffer(
            "double *", np.ascontiguousarray(locations, dtype=np.float64))
        potential = np.ascontiguousarray(
            np.empty((nlocations, )), dtype=np.float64)
        self._dftbpluslib.dftbp_get_elstat_potential(
            self._dftb_handler, nlocations, potential, locations)

        return potential

    def close(self):
        '''Finalizes the DFTB+ calculator.'''
        self._dftbpluslib.dftbp_final(self._dftb_handler)

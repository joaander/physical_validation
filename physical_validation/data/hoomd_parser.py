###########################################################################
#                                                                         #
#    physical_validation,                                                 #
#    a python package to test the physical validity of MD results         #
#                                                                         #
#    Written by Michael R. Shirts <michael.shirts@colorado.edu>           #
#               Pascal T. Merz <pascal.merz@colorado.edu>                 #
#                                                                         #
#    Copyright (C) 2012 University of Virginia                            #
#              (C) 2017 University of Colorado Boulder                    #
#                                                                         #
#    This library is free software; you can redistribute it and/or        #
#    modify it under the terms of the GNU Lesser General Public           #
#    License as published by the Free Software Foundation; either         #
#    version 2.1 of the License, or (at your option) any later version.   #
#                                                                         #
#    This library is distributed in the hope that it will be useful,      #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU    #
#    Lesser General Public License for more details.                      #
#                                                                         #
#    You should have received a copy of the GNU Lesser General Public     #
#    License along with this library; if not, write to the                #
#    Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,     #
#    Boston, MA 02110-1301 USA                                            #
#                                                                         #
###########################################################################
r"""
hoomd_parser.py
"""
import warnings
import numpy
import json

from . import parser
# py2.7 compatibility
from .simulation_data import SimulationData
from .unit_data import UnitData
from .ensemble_data import EnsembleData
from .system_data import SystemData
from .observable_data import ObservableData
from .trajectory_data import TrajectoryData
# replace lines above by this when py2.7 support is dropped:
# from . import SimulationData, UnitData, EnsembleData, SystemData, ObservableData, TrajectoryData
from ..util import error as pv_error

try:
    import gsd
    import gsd.hoomd
except ImportError:
    gsd = None;

class HOOMDParser(parser.Parser):
    """
    HOOMDParser
    """

    @staticmethod
    def units():
        # HOOMD doesn't require that the user use any specific system of units, rather it works in a
        # self-consistent set of units. GROMACS kJ/mol, nm system is reasonable for that
        # Pressure is the only unit requiring conversion
        return UnitData(
            kb=8.314462435405199e-3,
            energy_str='kJ/mol',
            energy_conversion=1.0,
            length_str='nm',
            length_conversion=1.0,
            volume_str='nm^3',
            volume_conversion=1.0,
            temperature_str='K',
            temperature_conversion=1.0,
            pressure_str='kJ/mol/nm^3',
            pressure_conversion=16.605391,
            time_str='ps',
            time_conversion=1.0)

    def __init__(self):
        r"""
        Create a HOOMDParser object
        """
        super(HOOMDParser, self).__init__()

    def get_simulation_data(self,
                            metadata_filename,
                            gsd_filename,
                            log_filename):
        r"""

        Parameters
        ----------
        meta: str, optional
            A string pointing to a .json metadata output file
        gsd: str, optional
            A string pointing to a .gsd trajectory file
        log: str, optional
            A string pointing to a .log thermodynamic quantity output file

        Returns
        -------
        result: SimulationData
            A SimulationData filled with the results of the simulation as described by
            the provided HOOMD files.

        """
        result = SimulationData()
        result.units = self.units()

        # trajectory, error out of the gsd python package is not available
        if gsd is None:
            raise pv_error.InputError('gsd',
                                      'gsd python package not found.')
        position = []
        velocity = []
        traj = gsd.hoomd.open(name=gsd_filename, mode='rb')
        for frame in traj:
            position.append(numpy.array(frame.particles.position));
            velocity.append(numpy.array(frame.particles.velocity));

        result.trajectory = TrajectoryData(
            position,
            velocity)

        with open(metadata_filename) as json_data:
            metadata = json.load(json_data)

        result.dt = float(metadata['hoomd.md.integrate.mode_standard']['dt'])

        system = SystemData()
        system.natoms = traj[0].particles.N;    # assumes N doesn't change
        system.nconstraints = 0;
        system.ndof_reduction_tra = 3
        system.ndof_reduction_rot = 0
        system.mass = traj[0].particles.mass;

        # TODO: handle molecules/constraints/bonds
        # Note: HOOMD has internal knowledge of the number of degrees of freedom, it seems that system.nconstraints
        # is only used to compute the degrees of freedom of the system. We could abuse this field and simply
        # insert a nconstraints value appropriate to set the number of degrees of freedom to match hoomd's.
        # this would also allow for testing 2D and simulations with orientation degrees of freedom (though
        # not with the equipartition analysis engine)

        # read log
        log_data = numpy.genfromtxt(fname=log_filename, names=True);

        volume = None;
        pressure = None;
        temperature = None;
        conserved_quantity = None;

        if 'hoomd.md.integrate.nvt' in metadata:
            ens = 'NVT';
            volume = log_data['volume'][0]
            temperature = metadata['hoomd.md.integrate.nvt'][0]['kT']/result.units.kb
            conserved_quantity = numpy.array(log_data['kinetic_energy'] + log_data['potential_energy'] + log_data['nvt_mtk_reservoir_energy_all'])
        elif 'hoomd.md.integrate.berendsen' in metadata:
            ens = 'NVT';
            volume = log_data['volume'][0]
            temperature = metadata['hoomd.md.integrate.berendsen'][0]['kT']/result.units.kb
        elif 'hoomd.md.integrate.langevin' in metadata:
            ens = 'NVT';
            volume = log_data['volume'][0]
            temperature = metadata['hoomd.md.integrate.langevin'][0]['kT']/result.units.kb
        elif 'hoomd.md.integrate.npt' in metadata:
            ens = 'NPT';
            pressure = metadata['hoomd.md.integrate.npt'][0]['S'][0]
            temperature = metadata['hoomd.md.integrate.npt'][0]['kT']/result.units.kb
            conserved_quantity = numpy.array(log_data['kinetic_energy'] + log_data['potential_energy'] + log_data['npt_thermostat_energy'] + log_data['npt_barostat_energy'])
        elif 'hoomd.md.integrate.nve' in metadata:
            ens = 'NVE';
            volume = log_data['volume'][0]
            conserved_quantity = numpy.array(log_data['kinetic_energy'] + log_data['potential_energy'])
        else:
            raise pv_error.InputError('metadata',
                                      'Integrator not found.')

        result.ensemble = EnsembleData(
            ens,
            natoms=system.natoms,
            volume=volume,
            pressure=pressure,
            temperature=temperature
        )

        result.observables = ObservableData()
        result.observables['kinetic_energy'] = numpy.array(log_data['kinetic_energy'])
        result.observables['potential_energy'] = numpy.array(log_data['potential_energy'])
        result.observables['total_energy'] = numpy.array(log_data['kinetic_energy'] + log_data['potential_energy'])
        if conserved_quantity is not None:
            result.observables['constant_of_motion'] = conserved_quantity;
        result.observables['volume'] = numpy.array(log_data['volume'])
        result.observables['pressure'] = numpy.array(log_data['pressure'])
        result.observables['temperature'] = numpy.array(log_data['temperature']/result.units.kb)

        result.system = system

        return result

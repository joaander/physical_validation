import hoomd, hoomd.md
import math

kT = dict(low=1.2, high=1.27)
P = dict(low=0.163, high=0.232)

for name in kT:

    hoomd.context.initialize()

    # initialize system configuration
    phi=0.4
    unitcell=hoomd.lattice.sc(a=math.pi/6/phi, type_name='A')
    hoomd.init.create_lattice(unitcell=unitcell, n=13)
    all = hoomd.group.all()

    # define interactions
    nl = hoomd.md.nlist.cell()
    lj = hoomd.md.pair.lj(r_cut=2.5, nlist=nl)
    lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0, r_on=2.0)
    lj.set_params(mode='xplor')

    # use langevin integrator to thermalize velocities
    # note: hoomd v2.3 will include a random velocity command that can replace this
    hoomd.md.integrate.mode_standard(dt=0.002)
    langevin = hoomd.md.integrate.langevin(group=all, kT=kT[name], seed=4)
    hoomd.run(100)
    zero = hoomd.md.update.zero_momentum(period=1)
    hoomd.run(1)
    zero.disable()
    langevin.disable()

    # NVT integration
    npt = hoomd.md.integrate.npt(group=all, kT=kT[name], tau=0.5, P=P[name], tauP=0.5)

    # equilibrate
    hoomd.run(100e3)

    # sample
    hoomd.analyze.log(filename="{0}.log".format(name),
                      quantities=['kinetic_energy', 'potential_energy', 'volume', 'pressure', 'temperature', 'npt_thermostat_energy', 'npt_barostat_energy'],
                      period=250,
                      overwrite=True)
    hoomd.dump.gsd(filename="{0}.gsd".format(name),
                   period=250,
                   group=all,
                   phase=0,
                   dynamic=['momentum'],
                   overwrite=True)
    hoomd.run(1e6)
    hoomd.meta.dump_metadata(filename="{0}.json".format(name))


import hoomd, hoomd.md
import math

kT = dict(low=1.0, high=1.03)

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
    lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0)
    lj.set_params(mode='shift')

    # use langevin integrator to thermalize velocities
    # note: hoomd v2.3 will include a random velocity command that can replace this
    hoomd.md.integrate.mode_standard(dt=0.005)
    langevin = hoomd.md.integrate.langevin(group=all, kT=kT[name], seed=4)
    zero = hoomd.md.update.zero_momentum(period=1e6, phase=100)
    hoomd.run(100)
    langevin.disable()

    # NVT integration
    nvt = hoomd.md.integrate.nvt(group=all, kT=kT[name], tau=0.5)

    # equilibrate
    hoomd.run(10e3)

    # sample
    hoomd.analyze.log(filename="{0}.log".format(name),
                      quantities=['kinetic_energy', 'potential_energy', 'volume', 'pressure', 'temperature', 'nvt_mtk_reservoir_energy_all'],
                      period=500,
                      overwrite=True)
    hoomd.dump.gsd(filename="{0}.gsd".format(name),
                   period=500,
                   group=all,
                   phase=0,
                   dynamic=['momentum'],
                   overwrite=True)
    hoomd.run(2.5e6)
    hoomd.meta.dump_metadata(filename="{0}.json".format(name))


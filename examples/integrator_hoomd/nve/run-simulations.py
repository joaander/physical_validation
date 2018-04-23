import hoomd, hoomd.md
import math

dt_list = [0.004, 0.002, 0.001, 0.0005, 0.00025]
steps_base=10000

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
hoomd.md.integrate.mode_standard(dt=dt_list[0])
langevin = hoomd.md.integrate.langevin(group=all, kT=1.2, seed=4)
hoomd.run(1000)
zero = hoomd.md.update.zero_momentum(period=1)
hoomd.run(1)
zero.disable()
langevin.disable()

hoomd.dump.gsd('init.gsd', period=None, group=all, overwrite=True);

for idx,dt in enumerate(dt_list):
    hoomd.context.initialize()

    steps = steps_base * dt_list[0]/dt
    log_period = int(dt_list[0]/dt)
    name = '{0}'.format(idx)

    # read initial condition
    hoomd.init.read_gsd('init.gsd');
    all = hoomd.group.all()

    # define interactions
    nl = hoomd.md.nlist.cell()
    lj = hoomd.md.pair.lj(r_cut=2.5, nlist=nl)
    lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0, r_on=2.0)
    lj.set_params(mode='xplor')

    # NVE integration
    hoomd.md.integrate.mode_standard(dt=dt)
    hoomd.md.integrate.nve(group=all)

    # equilibrate
    hoomd.run(steps)

    # sample
    hoomd.analyze.log(filename="{0}.log".format(name),
                      quantities=['kinetic_energy', 'potential_energy', 'volume', 'pressure', 'temperature'],
                      period=log_period,
                      overwrite=True)
    hoomd.dump.gsd(filename="{0}.gsd".format(name),
                   period=steps,
                   group=all,
                   phase=0,
                   dynamic=['momentum'],
                   overwrite=True)
    hoomd.run(steps)
    hoomd.meta.dump_metadata(filename="{0}.json".format(name))


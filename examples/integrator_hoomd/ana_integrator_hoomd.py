import physical_validation as pv
import os

systems = ['nve', 'nvt']

# change this to fit to your GROMACS installation
parser = pv.data.HOOMDParser()

for sys in systems:
    print('### Analyzing system ' + sys)
    print('## Reading results')
    res = []

    for n in range(0, 5):
        res.append(parser.get_simulation_data(
                   metadata_filename=os.path.join(sys, '{0}.json'.format(n)),
                   gsd_filename=os.path.join(sys, '{0}.gsd'.format(n)),
                   log_filename=os.path.join(sys, '{0}.log'.format(n))
                   ))

    # make plot directory
    if not os.path.exists('ana_integrator_hoomd_plots'):
        os.makedirs('ana_integrator_hoomd_plots')
    sysplot = os.path.join('ana_argon_plots', sys)

    print('## Validating integrator convergence')
    pv.integrator.convergence(res, verbose=True)
    #                         filename=sysplot)
    print()

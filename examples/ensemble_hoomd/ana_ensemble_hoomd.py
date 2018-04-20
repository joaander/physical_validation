import physical_validation as pv
import os

systems = ['be', 'nvt']

parser = pv.data.HOOMDParser()

for sys in systems:
    print('### Analyzing system ' + sys)
    print('## Reading lower temperature result')
    res_low = parser.get_simulation_data(
        metadata_filename=os.path.join(sys, 'low.json'),
        gsd_filename=os.path.join(sys, 'low.gsd'),
        log_filename=os.path.join(sys, 'low.log')
    )
    print('## Reading high temperature result')
    res_high = parser.get_simulation_data(
        metadata_filename=os.path.join(sys, 'high.json'),
        gsd_filename=os.path.join(sys, 'high.gsd'),
        log_filename=os.path.join(sys, 'high.log')
    )

    if not os.path.exists('ana_ensemble_hoomd_plots'):
        os.makedirs('ana_ensemble_hoomd_plots')
    sysplot = os.path.join('ana_ensemble_hoomd_plots', sys)

    print('\n## Validating kinetic energy distribution (strict)')
    print('# Low T:')
    pv.kinetic_energy.distribution(res_low, verbosity=2, strict=True,
                                   filename=sysplot + '_low_mb')
    print('# High T:')
    pv.kinetic_energy.distribution(res_high, verbosity=2, strict=True,
                                   filename=sysplot + '_high_mb')

    print('\n## Validating kinetic energy distribution (non-strict)')
    print('# Low T:')
    pv.kinetic_energy.distribution(res_low, verbosity=2, strict=False)
    print('# High T:')
    pv.kinetic_energy.distribution(res_high, verbosity=2, strict=False)

    print('\n## Validating ensemble')
    if 'pr' in sys:
        # can't plot in 2d...
        quantiles = pv.ensemble.check(res_low, res_high, verbosity=2)
    else:
        quantiles = pv.ensemble.check(res_low, res_high, verbosity=2,
                                      filename=sysplot + '_ensemble')
    if len(quantiles) == 1:
        q_str = '{:.1f}'.format(quantiles[0])
    else:
        q_str = '('
        for q in quantiles:
            q_str += '{:.1f}, '.format(q)
        q_str = q_str[:-2] + ')'
    print('Calculated slope is ' + q_str +
          ' quantiles from the true slope')
    print('\n')

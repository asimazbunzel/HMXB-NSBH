
import sys

import biaswise



def get_value_from_cc_file(fname: str, string: str):
    '''return value associated with string

    Parameters
    ----------
    fname : `string`
       Name of file where to search for information.

    string : `string`
       What to look for in `fname`.

    Returns
    -------
    value : `string/float`
       Formatted returned value.
    '''

    with open(fname, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if string in line:
            if string == 'sn_model':
                return line.strip().split()[1]
            else:
                return float(line.strip().split()[1])


def get_cc_info(star_fname: str, binary_fname: str):
    '''get all important parameters of star & binary at core-collapse

    Parameters
    ----------
    star_fname : `string`
       Name of file with info of collapsing star.

    binary_fname : `string`
       Name of file with info on binary system at core-collapse.

    Returns
    -------
    CollapseInfo : `dictionary`
       Information of binary system and collapsing star at core-collapse.
    '''

    m1_pre_cc = get_value_from_cc_file(star_fname, 'mass_pre_cc')
    m1_c_core_mass = get_value_from_cc_file(star_fname, 'c_core_mass_pre_cc')
    m1_remnant_mass = get_value_from_cc_file(star_fname, 'remnant_mass')
    m1_fallback_fraction = get_value_from_cc_file(star_fname, 'fallback_fraction')
    m2 = get_value_from_cc_file(binary_fname, 'companion_mass')
    P = get_value_from_cc_file(binary_fname, 'period_pre_cc')

    # store everything in a dictionary
    CollapseInfo = {'m1': m1_pre_cc, 'm1_c_core_mass': m1_c_core_mass,
            'm1_remnant_mass': m1_remnant_mass, 'm1_fallback_fraction': m1_fallback_fraction,
            'm2': m2, 'P': P}

    return CollapseInfo


if __name__ == '__main__':

    # get filename where kick info on star is located
    star_info_fname = sys.argv[1]
    # same but for binary info
    binary_info_fname = sys.argv[2]
    # folder where grid will be saved
    natal_kicks_folder = sys.argv[3]

    # collect important info at core-collapse
    Info = get_cc_info(star_info_fname, binary_info_fname)

    # compute natal kicks & post-kick binary parameters
    binary = biaswise.binary.BinarySystem(m1=Info['m1'], m1_core_mass=Info['m1_c_core_mass'],
            m1_remnant_mass=Info['m1_remnant_mass'],
            m1_fallback_fraction=Info['m1_fallback_fraction'], m2=Info['m2'], P=Info['P'])

    # define type of distribution of natal kicks & compute them
    binary.set_natal_kick_distribution(n_trials=50000, distribution_id='Maxwell',
            kick_scaling=lambda x: (1-binary.m1_fallback_fraction)*x)
    binary.get_natal_kick_distribution()

    # compute orbital parameters & make a grid out of them
    binary.get_orbital_distribution(verbose=True)
    binary.get_post_kick_grid(use_unbounded_for_norm=False, verbose=True)

    # save grid to a file
    fname = f'{natal_kicks_folder}/grid.data'
    binary.save_target_grid(fname=fname)

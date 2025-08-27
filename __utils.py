import numpy as np

from .__basics import *


def default_parameters():
    param = {}

    #### [STAR] ####
    param['Rs'] = None
    param['Loggs'] = None
    param['meta'] = None
    param['Ts'] = None
    param['Ts_err'] = None
    param['Ts_het'] = None
    param['het_frac'] = None
    param['Ts_spot'] = None
    param['spot_frac'] = None
    param['Ts_fac'] = None
    param['fac_frac'] = None

    #### [PLANET] ####
    param['name_p'] = None
    param['major-a'] = None
    param['Rp'] = None
    param['Mp'] = None
    param['Tp'] = None

    #### [STAR FITTING PARAMETERS] ####
    param['incl_star_activity'] = False
    param['stellar_activity_parameters'] = 5
    param['light_star_mods'] = False
    param['stellar_spec_dir'] = None

    #### [PLANETARY FITTING PARAMETERS] ####
    param['fit_Rp'] = True
    param['fit_Mp'] = False
    param['fit_T'] = False
    param['TP_profile'] = None
    param['P_upper'] = 8
    param['T_lower'] = 100
    param['T_upper'] = 2000

    #### [ATMOSPHERIC FITTING PARAMETERS] ####
    param['bare_rock'] = False
    param['clr_prior'] = 'uniform'
    param['gas_par_space'] = 'clr'
    param['vmr_upper'] = -0.3
    param['fit_H2O'] = False
    param['fit_CH4'] = False
    param['fit_C2H2'] = False
    param['fit_C2H4'] = False
    param['fit_C2H6'] = False
    param['fit_C4H2'] = False
    param['fit_CH3Cl'] = False
    param['fit_CH3SH'] = False
    param['fit_NH3'] = False
    param['fit_HCN'] = False
    param['fit_H2S'] = False
    param['fit_SO2'] = False
    param['fit_SO3'] = False
    param['fit_CS2'] = False
    param['fit_CO'] = False
    param['fit_CO2'] = False
    param['fit_N2O'] = False
    param['fit_OCS'] = False
    param['fit_DMS'] = False
    param['fit_PH3'] = False
    param['fit_N2'] = False
    param['fit_H2'] = False
    param['gas_fill'] = None
    param['incl_clouds'] = False
    param['cloud_type'] = None
    param['incl_haze'] = False
    param['haze_type'] = None

    #### [DATASET FITTING PARAMETERS] ####
    param['fit_offset'] = False
    param['offset_range'] = 100

    #### [ExoTR PARAMETER] ####
    param['optimizer'] = None
    param['opac_data'] = None
    param['use_float32'] = False
    param['opac_dir'] = None

    #### [MULTINEST PARAMETERS] ####
    param['multimodal'] = True
    param['max_modes'] = 100
    param['nlive_p'] = 600
    param['ev_tolerance'] = 0.5
    param['multinest_resume'] = True
    param['multinest_verbose'] = False

    #### [POST-PROCESSING PARAMETERS] ####
    param['calc_likelihood_data'] = False
    param['n_likelihood_data'] = 10240
    param['extended_wl_plot'] = False
    param['plot_models'] = False
    param['plot_posterior'] = False
    param['plot_elpd_stats'] = False
    param['elpd_reference'] = None
    param['truths'] = None

    #### [SYNTHESIZING SPECTRUM PARAMETERS] ####
    param['add_noise'] = False
    param['gaussian_noise'] = True
    param['use_noise_file'] = False
    param['noise_file'] = None
    param['snr'] = 20
    param['return_bins'] = False

    param['wkg_dir'] = os.getcwd()
    param['supported_molecules'] = ['H2O', 'CH4', 'C2H2', 'C2H4', 'C2H6', 'C4H2', 'CH3Cl', 'CH3SH', 'NH3', 'HCN', 'H2S', 'SO2', 'SO3', 'CS2', 'CO', 'CO2', 'N2O', 'OCS', 'DMS', 'PH3', 'N2', 'H2']

    for mol in param['supported_molecules']:
        param[mol + '_contribution'] = True
    param['cld_contribution'] = True
    param['CIA_contribution'] = True
    param['Rayleigh_contribution'] = True
    param['star_act_contribution'] = True
    param['haze_contribution'] = True

    mm = {'H': 1.00784, 'He': 4.002602, 'C': 12.0107, 'N': 14.0067, 'O': 15.9994, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453}
    mm['H2'] = 2. * mm['H']
    mm['H2O'] = (2. * mm['H']) + mm['O']
    mm['CH4'] = mm['C'] + (4. * mm['H'])
    mm['C2H2'] = (2. * mm['C']) + (2. * mm['H'])
    mm['C2H4'] = (2. * mm['C']) + (4. * mm['H'])
    mm['C2H6'] = (2. * mm['C']) + (6. * mm['H'])
    mm['C4H2'] = (4. * mm['C']) + (2. * mm['H'])
    mm['CH3Cl'] = (mm['C'] + (3. * mm['H']) + mm['Cl'])
    mm['CH3SH'] = (mm['C'] + (4. * mm['H']) + mm['S'])
    mm['NH3'] = mm['N'] + (3. * mm['H'])
    mm['HCN'] = mm['H'] + mm['C'] + mm['N']
    mm['H2S'] = (2. * mm['H']) + mm['S']
    mm['SO2'] = mm['S'] + (2. * mm['O'])
    mm['SO3'] = mm['S'] + (3. * mm['O'])
    mm['CS2'] = mm['C'] + (2. * mm['S'])
    mm['CO'] = mm['C'] + mm['O']
    mm['CO2'] = mm['C'] + (2. * mm['O'])
    mm['N2O'] = (2. * mm['N']) + mm['O']
    mm['OCS'] = mm['O'] + mm['C'] + mm['S']
    mm['DMS'] = (2. * mm['C']) + (6. * mm['H']) + mm['S']
    mm['PH3'] = mm['P'] + (3. * mm['H'])
    mm['O2'] = 2. * mm['O']
    mm['O3'] = 3. * mm['O']
    mm['N2'] = 2. * mm['N']

    param['mm'] = mm

    param['formatted_labels'] = {}
    for mol in param['supported_molecules']:
        if mol == 'H2O':
            param['formatted_labels'][mol] = "Log(H$_2$O)"
        if mol == 'CH4':
            param['formatted_labels'][mol] = "Log(CH$_4$)"
        if mol == 'C2H2':
            param['formatted_labels'][mol] = "Log(C$_2$H$_2$)"
        if mol == 'C2H4':
            param['formatted_labels'][mol] = "Log(C$_2$H$_4$)"
        if mol == 'C2H6':
            param['formatted_labels'][mol] = "Log(C$_2$H$_6$)"
        if mol == 'C4H2':
            param['formatted_labels'][mol] = "Log(C$_4$H$_2$)"
        if mol == 'CH3Cl':
            param['formatted_labels'][mol] = "Log(CH$_3$Cl)"
        if mol == 'CH3SH':
            param['formatted_labels'][mol] = "Log(CH$_3$SH)"
        if mol == 'NH3':
            param['formatted_labels'][mol] = "Log(NH$_3$)"
        if mol == 'HCN':
            param['formatted_labels'][mol] = "Log(HCN)"
        if mol == 'H2S':
            param['formatted_labels'][mol] = "Log(H$_2$S)"
        if mol == 'SO2':
            param['formatted_labels'][mol] = "Log(SO$_2$)"
        if mol == 'SO3':
            param['formatted_labels'][mol] = "Log(SO$_3$)"
        if mol == 'CS2':
            param['formatted_labels'][mol] = "Log(CS$_2$)"
        if mol == 'CO':
            param['formatted_labels'][mol] = "Log(CO)"
        if mol == 'CO2':
            param['formatted_labels'][mol] = "Log(CO$_2$)"
        if mol == 'N2O':
            param['formatted_labels'][mol] = "Log(N$_2$O)"
        if mol == 'OCS':
            param['formatted_labels'][mol] = "Log(OCS)"
        if mol == 'DMS':
            param['formatted_labels'][mol] = "Log(DMS)"
        if mol == 'PH3':
            param['formatted_labels'][mol] = "Log(PH$_3$)"
        if mol == 'N2':
            param['formatted_labels'][mol] = "Log(N$_2$)"
        if mol == 'H2':
            param['formatted_labels'][mol] = "Log(H$_2$)"

    return param


def read_parfile(param, parfile=None):
    """
        This function reads a parameter file, processes its contents, and updates a given dictionary with the
        obtained parameters. It's similar to `take_standard_parameters`, but the parameter file is provided by the user.

        The function requires the parameter file name to be explicitly provided, otherwise it will exit the
        program with an error message.

        It skips lines that start with '%' or are empty and interprets the remaining lines as key-value pairs separated
        by tabs. It handles 'mol', 'mol_vmr', and 'range_mol' keys similarly to the 'take_standard_parameters' function.

        The function also makes specific modifications based on certain keys. For instance, it adjusts parameters
        related to the gas_fill type, cloud inclusion, and optimization method. If 'incl_clouds' is True, it sets
        the fitting flags based on the cloud type. It also creates the output directory if it doesn't exist.

        Args:
            param (dict): The dictionary to update with the parameters from the file.
            parfile (str, optional): The name of the parameter file. If not provided, the program will exit.

        Returns:
            dict: The updated dictionary with the parameters from the file, including 'wkg_dir' and 'out_dir' entries.
    """
    cwd = os.getcwd()
    if parfile is None:
        print('No parameter file provided. Please provide one.')
        sys.exit()
    else:
        with open(cwd + '/' + parfile, 'r') as file:
            paramfile = file.readlines()
        for i in paramfile:
            if i[0] == '%' or i[0] == '\n':
                pass
            else:
                paramline = list(i.split('\t'))
                paramline[-1] = paramline[-1][:-1]
                if len(paramline) >= 2:
                    try:
                        param[paramline[0]] = float(paramline[-1])
                    except ValueError:
                        if str(paramline[1]) == str(True):
                            param[paramline[0]] = bool(paramline[1])
                        elif str(paramline[1]) == str(False):
                            param[paramline[0]] = bool("")
                        elif str(paramline[1]) == str(None):
                            param[paramline[0]] = None
                        else:
                            param[paramline[0]] = str(paramline[1])

                    if paramline[0] == 'file_output_name':
                        param[paramline[0]] = str(paramline[1])
                else:
                    paramline = str(paramline[0]).split()
                    if paramline[0] == 'mol':
                        param[paramline[0]] = paramline[1].split(',')
                    elif paramline[0] == 'mol_vmr' or paramline[0] == 'range_mol':
                        param[paramline[0]] = paramline[1].split(',')
                        for ob in range(0, len(param[paramline[0]])):
                            param[paramline[0]][ob] = float(param[paramline[0]][ob])
                        if paramline[0] == 'mol_vmr':
                            for num, mol in enumerate(param['mol']):
                                param['vmr_' + mol] = param['mol_vmr'][num]
                        else:
                            pass
                    else:
                        try:
                            param[paramline[0]] = float(paramline[-1])
                        except ValueError:
                            if str(paramline[1]) == str(True):
                                param[paramline[0]] = bool(paramline[1])
                            elif str(paramline[1]) == str(False):
                                param[paramline[0]] = bool("")
                            elif str(paramline[1]) == str(None):
                                param[paramline[0]] = None
                            else:
                                param[paramline[0]] = str(paramline[1])

    param['wkg_dir'] = cwd + '/'
    param['out_dir'] = param['wkg_dir'] + param['output_directory']
    try:
        os.mkdir(param['out_dir'])
    except FileExistsError:
        pass
    del param['output_directory']

    os.system('cp ' + cwd + '/' + parfile + ' ' + param['out_dir'])

    if not param['fit_Mp']:
        try:
            param['Mp'] += 0.0
        except KeyError:
            print('ERROR - Please specify the planetary mass by using "Mp" key in the parameter file.')
            sys.exit()
    else:
        try:
            param['Mp'] += 0.0
        except (KeyError, TypeError):
            if param['Mp'] is None:
                pass
            else:
                print('ERROR - Please specify the planetary mass or set it to None in the parameter file; use "Mp" as key.')
                sys.exit()

        if param['Mp'] is None:
            param['Mp_prior'] = 'uniform'
        else:
            try:
                param['Mp'] += 0.0
                param['Mp_err'] += 0.0
                param['Mp_orig'] = param['Mp'] + 0.0
                param['Mp_prior'] = 'gaussian'
            except KeyError:
                print('ERROR - Please specify the planetary mass ("Mp") and the associated errorbar ("Mp_err") to use a gaussian distribution. Otherwise, set "Mp" to "None" to use a uniform distribution.')
                sys.exit()

    if param['TP_profile'] is not None:
        param['TP_profile'] = np.loadtxt(param['TP_profile'])

    if param['gas_fill'] is None:
        param['gas_fill'] = 'H2'

    for mol in param['supported_molecules']:
        if param['gas_fill'] == mol:
            param['fit_' + mol] = False

    if param['bare_rock']:
        param['gas_fill'] = None
        param['incl_clouds'] = False
        param['incl_haze'] = False
        for mol in param['supported_molecules']:
            param['fit_' + mol] = False

    if param['incl_clouds']:
        if param['cloud_type'] == 'water':
            param['fit_wtr_cld'] = True
            param['fit_amn_cld'] = False
            param['fit_gen_cld'] = False
        elif param['cloud_type'] == 'ammonia':
            param['fit_wtr_cld'] = False
            param['fit_amn_cld'] = True
            param['fit_gen_cld'] = False
        else:
            param['fit_wtr_cld'] = False
            param['fit_amn_cld'] = False
            param['fit_gen_cld'] = True
    else:
        param['fit_wtr_cld'] = False
        param['fit_amn_cld'] = False
        param['fit_gen_cld'] = False
        param['CR_H2O'] = 1.0
        param['CR_NH3'] = 1.0

    if param['incl_clouds'] and param['fit_wtr_cld']:
        try:
            param['vmr_H2O'] += 0.0
        except KeyError:
            if param['fit_H2O']:
                pass
            else:
                if param['fit_wtr_cld']:
                    print('You have selected to use the water cloud model, but the volume mixing ratio of water is not provided nor it is fitted. Switching to general cloud model.')
                    param['fit_wtr_cld'] = False
                    param['fit_gen_cld'] = True
                    param['P_top'] = param['pH2O'] + 0.0
                else:
                    param['fit_gen_cld'] = True
                param['vmr_H2O'] = 10. ** (-20.)

    if param['incl_haze']:
        if param['haze_type'] == 'tholin':
            param['fit_tholin'] = True
            param['fit_soot'] = False
        elif param['haze_type'] == 'soot':
            param['fit_tholin'] = False
            param['fit_soot'] = True
    else:
        param['fit_tholin'] = False
        param['fit_soot'] = False

    if param['optimizer'] == 'multinest':
        param['nlive_p'] = int(param['nlive_p'])
        param['max_modes'] = int(param['max_modes'])
    elif param['optimizer'] == 'dynesty':
        param['nlive_init'] = int(param['nlive_init'])
        param['nlive_batch'] = int(param['nlive_batch'])
    else:
        pass

    if param['optimizer'] is not None:
        param['fit_molecules'] = []
        for mol in param['supported_molecules']:
            if param['gas_par_space'] == 'vmr':
                param['vmr_' + mol] = 0.0
            if param['fit_' + mol]:
                param['fit_molecules'].append(mol)

    if param['incl_star_activity']:
        if param['optimizer'] is None:
            try:
                param['Ts_phot'] += 0.0
            except KeyError:
                param['Ts_phot'] = param['Ts'] + 0.0

        param['stellar_activity_parameters'] = int(param['stellar_activity_parameters'])

        if param['meta'] is None:
            param['meta'] = 0.0

        if not os.path.isdir(param['pkg_dir'] + 'PHOENIX_models/'):
            if os.path.isdir(param['pkg_dir'] + 'PHOENIX_models_light/'):
                param['light_star_mods'] = True

    return param


def par_and_calc(param):
    """
        This function modifies and adds entries to a parameter dictionary based on existing values. It is intended
        to be used after reading a parameter file with `read_parfile`.

        If 'fit_T' is False and 'Tp' is not already defined, it calculates 'Tp' using a formula that includes
        the stellar radius ('Rs'), the semi-major axis of the planet's orbit ('major-a'), and the stellar temperature ('Ts').

        It also defines 'P_standard' and 'P' arrays which are logarithmically spaced pressure arrays.

        Args:
            param (dict): The dictionary containing the parameters to be updated.

        Returns:
            dict: The updated dictionary with modified and added parameters.
    """

    def add_boundary_knots(spline):
        """
        Add knots infinitesimally to the left and right.

        Additional intervals are added to have zero 2nd and 3rd derivatives,
        and to maintain the first derivative from whatever boundary condition
        was selected. The spline is modified in place.
        """
        # determine the slope at the left edge
        leftx = spline.x[0]
        lefty = spline(leftx)
        leftslope = spline(leftx, nu=1)

        # add a new breakpoint just to the left and use the
        # known slope to construct the PPoly coefficients.
        leftxnext = np.nextafter(leftx, leftx - 1)
        leftynext = lefty + leftslope * (leftxnext - leftx)
        leftcoeffs = np.array([0, 0, leftslope, leftynext])
        spline.extend(leftcoeffs[..., None], np.r_[leftxnext])

        # repeat with additional knots to the right
        rightx = spline.x[-1]
        righty = spline(rightx)
        rightslope = spline(rightx, nu=1)
        rightxnext = np.nextafter(rightx, rightx + 1)
        rightynext = righty + rightslope * (rightxnext - rightx)
        rightcoeffs = np.array([0, 0, rightslope, rightynext])
        spline.extend(rightcoeffs[..., None], np.r_[rightxnext])

    # planet
    if not param['bare_rock']:
        if not param['fit_T'] and param['TP_profile'] is None:
            try:
                param['Tp'] += 0.0
            except TypeError:
                t1 = ((param['Rs'] * const.R_sun) / (2. * param['major-a'] * const.au)) ** 0.5
                param['Tp'] = t1 * ((1 - 0.3) ** 0.25) * param['Ts']

        if param['TP_profile'] is not None:
            tp_prof = param['TP_profile'] + 0.0
            param['TP_profile'] = CubicSpline(tp_prof[:, 1], tp_prof[:, 0], bc_type='natural')
            add_boundary_knots(param['TP_profile'])

    param['P_standard'] = 10. ** np.arange(-1.0, 12.1, step=0.1)
    param['P'] = 10. ** np.arange(-1.0, param['P_upper']+0.1, step=0.1)
    return param


def calc_mean_mol_mass(param):
    """
        This function calculates the mean molecular weight of a planet's atmosphere. The composition of the atmosphere
        is determined by the volume mixing ratios (vmr) of different molecules in the parameter dictionary. The molecular
        masses of individual molecules are defined within the function.

        Based on the fill gas type ('H2', 'N2', or 'CO2'), the function adjusts the corresponding vmr and calculates
        the mean molecular weight of the atmosphere considering other molecules too. The 'ret_mode' key is used to
        decide whether to print the vmr of the fill gas and the mean molecular weight or not.

        Args:
            param (dict): The dictionary containing the parameters including molecular volume mixing ratios (vmr),
                          fill gas type, H2O, and ret_mode flag. The dictionary is expected to be in a particular
                          format defined by previous processing steps.

        Returns:
            dict: The updated parameter dictionary including mean molecular weight and updated vmr values for the
                  corresponding fill gas.
    """
    param['vmr_' + param['gas_fill']] = np.ones(len(param['P']))
    param['mean_mol_weight'] = np.zeros(len(param['P']))
    for mol in param['fit_molecules']:
        param['vmr_' + mol] = np.ones(len(param['P'])) * param['vmr_' + mol]
        param['vmr_' + param['gas_fill']] -= param['vmr_' + mol]
        if len(np.where(param['vmr_' + param['gas_fill']] < 0)[0]) > 0:
            idx = np.where(param['vmr_' + param['gas_fill']] < 0)[0]
            param['vmr_' + param['gas_fill']][idx] = (10.**(-12.))
        param['mean_mol_weight'] += (param['vmr_' + mol] * param['mm'][mol])
    if param['gas_fill'] == 'H2':
        param['mean_mol_weight'] += (param['vmr_' + param['gas_fill']] * (0.8 * param['mm']['H2'] + 0.2 * param['mm']['He']))
    else:
        param['mean_mol_weight'] += (param['vmr_' + param['gas_fill']] * param['mm'][param['gas_fill']])

    if not param['ret_mode']:
        print('VMR ' + param['gas_fill'] + ' \t\t = \t' + str(param['vmr_' + param['gas_fill']][-1]))
        print('mu \t\t = \t' + str(param['mean_mol_weight'][-1]))

    return param


def load_input_spectrum(param):
    """
        This function loads an input spectrum from a file and adds it to the parameter dictionary. The input spectrum
        file's format and location are specified in the parameter dictionary.

        Depending on the 'ret_mode' key in the dictionary, it either loads a spectrum containing wavelength, transit
        depth, and error on transit depth or loads wavelength bins only.

        If the function cannot find a required input file, it will print a descriptive error message and then
        terminate the program.

        Args:
            param (dict): A dictionary that contains the parameters necessary to load the input spectrum.
                          The dictionary is expected to be in a particular format defined by previous processing steps.

        Returns:
            dict: The updated parameter dictionary with the added spectrum or wavelength bins.
    """
    if param['ret_mode']:
        try:
            spectrum = np.loadtxt(param['wkg_dir'] + param['spectrum'])
            param['spectrum'] = {}
            if len(spectrum[0, :]) == 3:
                wl_idx = 0
                param['spectrum']['bins'] = False
                param['min_wl'] = min(spectrum[:, 0])
                param['max_wl'] = max(spectrum[:, 0])
            else:
                wl_idx = 2
                param['spectrum']['wl_bins'] = spectrum[:, 0:3]  # wavelength bin_low in micron
                param['spectrum']['bins'] = True
                param['min_wl'] = min(param['spectrum']['wl_bins'][:, 0])
                param['max_wl'] = max(param['spectrum']['wl_bins'][:, 1])
                param['spectrum']['wl_low'] = param['spectrum']['wl_bins'][:, 0]
                param['spectrum']['wl_high'] = param['spectrum']['wl_bins'][:, 1]
            param['spectrum']['wl'] = spectrum[:, wl_idx]  # wavelength in micron
            param['spectrum']['T_depth'] = spectrum[:, wl_idx + 1]  # Transit depth
            param['spectrum']['error_T'] = spectrum[:, wl_idx + 2]  # error on Transit depth
            param['sorted_data_idx'] = np.argsort(param['spectrum']['wl'])
            multi_spec_mode = False
        except KeyError:
            try:
                _ = np.loadtxt(param['wkg_dir'] + param['spec1'])
                multi_spec_mode = True
            except KeyError:
                print(
                    'An input spectrum is required, in the parameter file, use the "spectrum" keyword followed by the path of the file.')
                sys.exit()

        if multi_spec_mode:
            n_spec = 1
            while True:
                try:
                    _ = np.loadtxt(param['wkg_dir'] + param['spec' + str(n_spec + 1)])
                    n_spec += 1
                except KeyError:
                    break

            spec1 = np.loadtxt(param['wkg_dir'] + param['spec1'])
            param['spectrum'] = {}
            if len(spec1[0, :]) == 3:
                wl_idx = 0
                param['spectrum']['bins'] = False
            else:
                wl_idx = 2
                param['spectrum']['wl_bins'] = spec1[:, 0:3]  # wavelength in micron
                param['spectrum']['bins'] = True
            param['spectrum']['wl'] = spec1[:, wl_idx]  # wavelength in micron
            param['spectrum']['error_T'] = spec1[:, wl_idx + 2]

            for i in range(1, n_spec):
                spec1 = np.loadtxt(param['wkg_dir'] + param['spec' + str(i + 1)])
                param['spectrum']['wl'] = np.append(param['spectrum']['wl'], spec1[:, wl_idx])  # wavelength in micron
                param['spectrum']['error_T'] = np.append(param['spectrum']['error_T'], spec1[:, wl_idx + 2])  # error on Transit depth
                if len(spec1[0, :]) > 3:
                    param['spectrum']['wl_bins'] = np.append(param['spectrum']['wl_bins'], spec1[:, 0:3], axis=0)
            if len(spec1[0, :]) == 3:
                param['min_wl'] = min(param['spectrum']['wl'])
                param['max_wl'] = max(param['spectrum']['wl'])
            else:
                param['min_wl'] = min(param['spectrum']['wl_bins'][:, 0])
                param['max_wl'] = max(param['spectrum']['wl_bins'][:, 1])

                param['spectrum']['wl_low'] = param['spectrum']['wl_bins'][:, 0]
                param['spectrum']['wl_high'] = param['spectrum']['wl_bins'][:, 1]

                del param['spectrum']['wl_bins']

            param['sorted_data_idx'] = np.argsort(param['spectrum']['wl'])

            if not param['fit_offset']:
                spec1 = np.loadtxt(param['wkg_dir'] + param['spec1'])
                param['spectrum']['T_depth'] = spec1[:, wl_idx + 1]  # Transit depth
                for i in range(1, n_spec):
                    spec1 = np.loadtxt(param['wkg_dir'] + param['spec' + str(i + 1)])
                    param['spectrum']['T_depth'] = np.append(param['spectrum']['T_depth'], spec1[:, wl_idx + 1])  # Transit depth
            else:
                param['n_offsets'] = n_spec - 1
                for i in range(0, param['n_offsets'] + 1):
                    spectrum = np.loadtxt(param['wkg_dir'] + param['spec' + str(i + 1)])
                    param['spectrum' + str(i + 1)] = {}
                    param['spectrum' + str(i + 1)]['T_depth'] = spectrum[:, wl_idx + 1]  # Transit depth
    else:
        try:
            spectrum = np.loadtxt(param['pkg_dir'] + 'Data/wl_bins/' + param['wave_file'] + '.dat')
            param['spectrum'] = {}
            try:
                param['spectrum']['wl_bins'] = spectrum + 0.0  # wavelength bins in micron
                param['spectrum']['wl'] = np.mean(param['spectrum']['wl_bins'], axis=1)  # wavelength in micron
                param['min_wl'] = min(param['spectrum']['wl_bins'][:, 0])
                param['max_wl'] = max(param['spectrum']['wl_bins'][:, 1])
                param['spectrum']['wl_low'] = param['spectrum']['wl_bins'][:, 0]
                param['spectrum']['wl_high'] = param['spectrum']['wl_bins'][:, 1]
                param['spectrum']['bins'] = True
                param['sorted_data_idx'] = np.argsort(param['spectrum']['wl'])
            except IndexError:
                param['spectrum']['wl'] = spectrum
                param['sorted_data_idx'] = np.argsort(param['spectrum']['wl'])
                param['spectrum']['bins'] = False
                param['min_wl'] = min(param['spectrum']['wl'])
                param['max_wl'] = max(param['spectrum']['wl'])
        except KeyError:
            print('No wavelength bins provided, using the native spectral resolution of the opacities.')
            param['wave_file'] = None
            param['spectrum']['wl'] = None
            param['spectrum']['bins'] = False

    return param


def find_nearest(array, value):
    """
        This function returns the index of the element in 'array' that is closest to 'value'.
        It handles arrays with NaN values by ignoring them and only considering the other elements of the array.

        Args:
            array (numpy.array): A numpy array in which the nearest value to 'value' is to be found.
            value (float): A numeric value. The function will find the index of the element in 'array' closest to this value.

        Returns:
            int: The index of the element in 'array' closest to 'value'.
    """
    idx = np.nanargmin(np.absolute(array - value))
    return idx


def cloud_pos(param):
    """
        Updates the atmospheric profile with a water cloud layer.

        Parameters:
        param (dict): A dictionary containing the model parameters.

        Returns:
        param (dict): The same dictionary with updated water mixing ratio profile.

        The function takes an atmospheric model parameter dictionary. If the model includes a water cloud layer,
        it calculates the position of the cloud in the atmosphere based on the cloud top pressure. The function then
        updates the water mixing ratio accordingly. If the model doesn't include a cloud layer, it assumes a constant
        water mixing ratio throughout the atmosphere.
    """
    if param['fit_wtr_cld']:
        pos_cldw = int(find_nearest(param['P_standard'], param['Pw_top']))

        if (param['cldw_depth'] + param['P_standard'][pos_cldw]) > param['P_standard'][-1]:
            param['cldw_depth'] = param['P_standard'][-1] - param['P_standard'][pos_cldw]

        pwbot = int(find_nearest(param['P_standard'], (param['cldw_depth'] + param['P_standard'][pos_cldw])))

        depth_w = pwbot - pos_cldw
        if depth_w == 0:
            depth_w = 1

        if not hasattr(param['vmr_H2O'], "__len__"):
            pass
        else:
            param['vmr_H2O'] = param['vmr_H2O'][-1] + 0.0
        watermix = np.ones((len(param['P_standard']))) * (param['CR_H2O'] * param['vmr_H2O'])
        dw = (np.log10(param['vmr_H2O']) - np.log10(param['CR_H2O'] * param['vmr_H2O'])) / depth_w
        for i in range(0, len(watermix)):
            if i <= pos_cldw:
                pass
            elif pos_cldw < i <= pos_cldw + depth_w:
                watermix[i] = 10. ** (np.log10(watermix[i - 1]) + dw)
            elif i > pos_cldw + depth_w:
                watermix[i] = watermix[i - 1]
        watermix = watermix[:len(param['P'])]
    else:
        watermix = np.ones((len(param['P']))) * param['vmr_H2O']

    param['vmr_H2O'] = watermix + 0.0
    return param


def ranges(param):
    """
        Sets the ranges for the fit parameters in the model.

        Parameters:
        param (dict): A dictionary containing the model parameters.

        Returns:
        param (dict): The same dictionary with added range parameters.

        This function takes an atmospheric model parameter dictionary.
        Depending on which molecules and parameters are set to be fitted,
        it sets the corresponding range parameters in the dictionary.
    """
    if param['fit_wtr_cld']:
        param['ptopw_range'] = [0.0, 7.5]  # Top pressure H2O
        param['dcldw_range'] = [0.0, 8.0]  # Depth H2O cloud
        param['crh2o_range'] = [-12.0, 0.0]  # Condensation Ratio H2O
    if param['fit_amn_cld']:
        param['ptopa_range'] = [0.0, 7.5]  # Top pressure NH3
        param['dclda_range'] = [0.0, 8.0]  # Depth HN3 cloud
        param['crnh3_range'] = [-12.0, 0.0]  # Condensation Ratio NH3
    if param['fit_gen_cld']:
        param['ptop_range'] = [0.0, 9.0]  # Top pressure
    if param['fit_tholin']:
        param['dtholin_range'] = [-3, 2]  # diameter of tholin haze particle
        param['vmrtholin_range'] = [-10, -1]
    if param['fit_soot']:
        param['dsoot_range'] = [-3, 2]  # diameter of soot haze particle
        param['vmrsoot_range'] = [-10, -1]
    if param['fit_T'] and param['Tp'] is None and param['TP_profile'] is None:
        param['tp_range'] = [param['T_lower'], param['T_upper']]  # Atmospheric equilibrium temperature
    elif param['fit_T'] and param['Tp'] is None and param['TP_profile'] is not None:
        param['tp_range'] = [-300.0, 300.0]
    elif param['fit_T'] and param['Tp'] is not None:
        param['tp_range'] = [max([100.0, param['Tp'] - 500.0]), param['Tp'] + 500.0]  # Atmospheric equilibrium temperature
    if len(param['fit_molecules']) > 0:
        if param['gas_par_space'] == 'clr':
            if param['clr_prior'] == 'uniform':
                param['gas_clr_range'] = [-25, 25]  # Centered-log-ratio standard range
        elif param['gas_par_space'] == 'vmr':
            param['gas_vmr_range'] = [-12, param['vmr_upper']]
    if param['fit_Rp']:
        param['rp_range'] = [param['Rp'] * 0.5, param['Rp'] * 2.0]  # Planetary Radius
    if param['fit_Mp']:
        if param['Mp_prior'] == 'gaussian':
            param['mp_range'] = [param['Mp_orig'] - (5.0 * param['Mp_err']), param['Mp_orig'] + (5.0 * param['Mp_err'])]
            if param['Mp_orig'] - (5.0 * param['Mp_err']) < 0.01:
                param['mp_range'][0] = 0.01
            if param['Mp_orig'] + (5.0 * param['Mp_err']) > 20.0:
                param['mp_range'][1] = 20.0
        else:
            param['mp_range'] = [0.01, 20.0]  # Planetary Mass
    if param['fit_offset']:
        for i in range(0, param['n_offsets']):
            param['off' + str(i + 1) + '_range'] = [-param['offset_range'] / 1e6, param['offset_range'] / 1e6]
    if param['incl_star_activity']:
        param['delta_range'] = [0.0, 0.5]  # Fraction of the star surface affected by activity
        param['Ts_phot_range'] = [param['Ts'] - (5.0 * param['Ts_err']), param['Ts'] + (5.0 * param['Ts_err'])]
        if param['stellar_activity_parameters'] == int(3):
            param['Ts_het_range'] = [1500, 1.2 * param['Ts']]

    return param


def define_modified_clr_prior(ngas, clr_prior_type):
    def linear_to_clr(x):
        """Transformation to the center log ratio space"""
        n = x.size
        gx = (np.prod(x)) ** (1 / n)
        clr = np.log(x / gx)
        return clr

    def uniform_distribution_to_clr(ng, fmin=1e-12, ns=int(1e6)):
        """generate a center log ratio distrubtion, corresponding to the log-uniform distribution of [fmin, 1]"""
        # generate uniform distributions in log space between 1e-12 and 1
        uniform_sample = 10.0 ** np.random.uniform(low=np.log10(fmin), high=0, size=(ns, ng))
        # only use samples which sum to less than 1
        inds = np.where(np.sum(uniform_sample, axis=1) <= 1)[0]

        filler = (10. ** 0) - np.sum(uniform_sample[inds, :], axis=1)
        filler = np.reshape(filler, (len(filler), 1))
        gases = np.append(uniform_sample[inds, :], filler, axis=1)

        # gases = uniform_sample[inds,:]
        nsg = gases.shape[0]

        # apply CLR transformation to all sets of gases
        clr = np.empty(gases.shape)
        for i in range(nsg):
            clr[i, :] = linear_to_clr(gases[i, :])

        # we return only one prior, as all should be identical
        return clr[:, 0]

    clr = uniform_distribution_to_clr(ngas)
    counts, bins = np.histogram(clr, bins=int(1e4), density=True)
    bns = np.concatenate((np.array([bins[:-1]]).T, np.array([bins[1:]]).T), axis=1)

    if clr_prior_type == 'hybrid' and ngas > 1:
        avg = np.median(counts[find_nearest(bins, -5):find_nearest(bins, 0)])
        std = np.std(counts[find_nearest(bins, -5):find_nearest(bins, 0)])
        uniform = np.random.uniform(-2 * std, 2 * std, size=(len(bins[:-1])))
        idxs = np.where(bins < -5)[0]
        counts[idxs] = uniform[idxs] + avg
        counts /= np.sum(counts * (bns[:, 1] - bns[:, 0]))

    cdf = np.cumsum(counts * (bns[:, 1] - bns[:, 0]))
    cdf[0] = 0.0
    cdf[-1] = 1.0
    bnss = np.mean(bns, axis=1)
    ppf = interp1d(cdf, bnss)

    return ppf


def pre_load_variables(param, for_plotting=False):
    """
        Loads opacity data from a MATLAB (.mat) file and updates the parameter dictionary.

        Parameters:
        param (dict): A dictionary containing the model parameters.

        Returns:
        param (dict): The updated parameter dictionary with opacity data.

        This function loads the opacity data from a MATLAB file into the parameter dictionary.
        It specifically excludes certain keys from the MATLAB file that are not related to the opacity data.
        It also initializes the opacity without clouds to zero.
    """
    if not param['bare_rock']:
        data = scipy.io.loadmat(param['pkg_dir'] + 'Data/opac/opac_082025.mat')
        opac_data_keys = []
        for i in data.keys():
            if i != '__header__' and i != '__globals__' and i != '__version__':
                param[i] = np.array(data[i]) + 0.0
                opac_data_keys.append(i)

        del data

        # delete the opacities of the molecules we don't need
        for mol in param['supported_molecules']:
            if mol not in param['fit_molecules'] + [param['gas_fill']]:
                del param['opac' + mol.lower()]

        if param['gas_fill'] == 'N2' or param['gas_fill'] == 'H2':
            param['opac' + param['gas_fill'].lower()] = np.zeros_like(param['opach2o'])
        if 'N2' in param['fit_molecules']:
            param['opacn2'] = np.zeros_like(param['opach2o'])
        if 'H2' in param['fit_molecules']:
            param['opach2'] = np.zeros_like(param['opach2o'])
        param['opac_data_keys'] = opac_data_keys
        param['opacaer_no_cloud'] = param['opacaerh2o'] * 0.0

        if param['opac_data'] is not None:
            for mol in param['fit_molecules'] + [param['gas_fill']]:
                # todo: we dont have to calculate opact, opacp, opacw with every loop, just once.
                if param['opac_dir'] is not None:
                  param['opact'], param['opacp'], param['opacw'], param['opac' + mol.lower()] = readcross(param['opac_dir'] + param['opac_data'] + '/opac' + mol + '.dat')
                else:
                  param['opact'], param['opacp'], param['opacw'], param['opac' + mol.lower()] = readcross(param['pkg_dir'] + 'Data/opac/' + param['opac_data'] + '/opac' + mol + '.dat')
                # making opacity files float32 to save space (cut memory usage by half)
                if param['use_float32']:
                    param['opac' + mol.lower()] = np.array(param['opac' + mol.lower()], dtype=np.float32)

        if not for_plotting:
            strt = find_nearest(param['opacw'][0] * 1e6, param['min_wl'] - 0.05) - 20
            end = find_nearest(param['opacw'][0] * 1e6, param['max_wl'] + 0.05) + 20
            param['opacw'] = (param['opacw'][0][strt:end]).reshape(1, -1)
            ## cut temperature range for opacity loading
            strtT = find_nearest(param['opact'][0],param['T_lower'])
            endT = find_nearest(param['opact'][0],param['T_upper'])
            param['opact'] = (param['opact'][0][strtT:endT+1]).reshape(1,-1)
            ## cut pressure range for opacity loading (just the upper bound)
            endP = find_nearest(param['opacp'][0],10**param['P_upper'])   # cut pressure list at value specified in retpar
            param['opacp'] = (param['opacp'][0][:endP+1]).reshape(1,-1)
            for mol in param['fit_molecules'] + [param['gas_fill']]:
                param['opac' + mol.lower()] = param['opac' + mol.lower()][:endP + 1, strtT:endT + 1, strt:end]


        if not param['fit_T']:
            P = np.array([param['P'][::-1]]).T
            n = len(P[:, 0])

            if param['TP_profile'] is None:
                T = np.ones_like(P) * param['Tp']
            else:
                T = param['TP_profile'](P)
            for i in range(0, len(P)):
                T[i] = min(max(100.0, T[i]), 2000.0)

            param['forward'] = {}
            param['forward']['PL'] = np.array([np.sqrt(P[0:n - 1, 0] * P[1:n, 0])]).T
            param['forward']['TL'] = np.array([(T[0:n - 1, 0] + T[1:n, 0]) / 2.0]).T

            I1 = np.ones(len(param['forward']['PL'])).reshape(len(param['forward']['PL']), 1)
            I2 = np.ones(len(param['opacw'][0]))
            param['forward']['S'] = {}

            # Calculate the molecular cross-sections
            for mol in param['fit_molecules'] + [param['gas_fill']]:
                param['forward']['S'][mol] = (interpn((param['opacp'][0], param['opact'][0], param['opacw'][0]), param['opac' + mol.lower()], np.array([param['forward']['PL'] * I2, param['forward']['TL'] * I2, I1 * param['opacw'][0]]).T)).T
                if param['use_float32']:
                    param['forward']['S'][mol] = np.array(param['forward']['S'][mol], dtype=np.float32)
                del param['opac' + mol.lower()]

    return param


def reso_range(start, finish, res, bins=False):
    wl_low = [start]
    res = 1. / res
    wl_high = [start + (start * res)]
    while wl_high[-1] < finish:
        wl_low.append(wl_high[-1])
        wl_high.append(wl_low[-1] + (wl_low[-1] * res))

    bns = np.array([wl_low, wl_high]).T

    if not bins:
        return np.mean(bns, axis=1)
    else:
        return bns


def retrieval_par_and_npar(param):
    """
        This function determines the parameters that will be retrieved based on the settings in the input parameter dictionary.
        For each parameter set to be fitted, it appends a corresponding label to the list of parameters and increments the
        gas counter for gas species to be fitted.

        Args:
            param (dict): A dictionary containing the parameter settings. Each entry corresponds to a boolean value
            indicating whether the respective parameter should be fitted.

        Returns:
            parameters (list): A list of string labels for the parameters to be retrieved.
            n_parameters (int): The number of parameters to be retrieved.
    """
    parameters = []
    if param['fit_offset']:
        for i in range(0, param['n_offsets']):
            parameters.append("offset$_" + str(i + 1) + "$")
    if param['fit_Rp']:
        parameters.append("R$_p$")
    if param['fit_Mp']:
        parameters.append("M$_p$")
    if param['fit_T'] and param['TP_profile'] is None:
        parameters.append("T$_p$")
    elif param['fit_T'] and param['TP_profile'] is not None:
        parameters.append("$\Delta$T$_p$")
    if param['fit_wtr_cld']:
        parameters.append("Log(P$_{top, H_2O}$)")
        parameters.append("Log(D$_{H_2O}$)")
        parameters.append("Log(CR$_{H_2O}$)")
    if param['fit_amn_cld']:
        parameters.append("Log(P$_{top, NH_3}$)")
        parameters.append("Log(D$_{NH_3}$)")
        parameters.append("Log(CR$_{NH_3}$)")
    if param['fit_gen_cld']:
        parameters.append("Log(P$_{top}$)")
    if param['fit_tholin']:
        parameters.append("Log(diam$_{tholin}$)")
        parameters.append("Log(vmr$_{tholin}$)")
    if param['fit_soot']:
        parameters.append("Log(diam$_{soot}$)")
        parameters.append("Log(vmr$_{soot}$)")
    if param['fit_H2O']:
        parameters.append("clr(H$_2$O)")
    if param['fit_CH4']:
        parameters.append("clr(CH$_4$)")
    if param['fit_C2H2']:
        parameters.append("clr(C$_2$H$_2$)")
    if param['fit_C2H4']:
        parameters.append("clr(C$_2$H$_4$)")
    if param['fit_C2H6']:
        parameters.append("clr(C$_2$H$_6$)")
    if param['fit_C4H2']:
        parameters.append("clr(C$_4$H$_2$)")
    if param['fit_CH3Cl']:
        parameters.append("clr(CH$_3$Cl)")
    if param['fit_CH3SH']:
        parameters.append("clr(CH$_3$SH)")
    if param['fit_NH3']:
        parameters.append("clr(NH$_3$)")
    if param['fit_HCN']:
        parameters.append("clr(HCN)")
    if param['fit_H2S']:
        parameters.append("clr(H$_2$S)")
    if param['fit_SO2']:
        parameters.append("clr(SO$_2$)")
    if param['fit_SO3']:
        parameters.append("clr(SO$_3$)")
    if param['fit_CS2']:
        parameters.append("clr(CS$_2$)")
    if param['fit_CO']:
        parameters.append("clr(CO)")
    if param['fit_CO2']:
        parameters.append("clr(CO$_2$)")
    if param['fit_N2O']:
        parameters.append("clr(N$_2$O)")
    if param['fit_OCS']:
        parameters.append("clr(OCS)")
    if param['fit_DMS']:
        parameters.append("clr(DMS)")
    if param['fit_PH3']:
        parameters.append("clr(PH$_3$)")
    if param['fit_N2']:
        parameters.append("clr(N$_2$)")
    if param['fit_H2']:
        parameters.append("clr(H$_2$)")
    if param['incl_star_activity']:
        if param['stellar_activity_parameters'] == int(3):
            parameters.append("$\delta$" + "$_{het}$")
            parameters.append("T$_{het}$")
            parameters.append("T$_{phot}$")
        elif param['stellar_activity_parameters'] == int(5):
            parameters.append("$\delta$" + "$_{spot}$")
            parameters.append("$\delta$" + "$_{fac}$")
            parameters.append("T$_{spot}$")
            parameters.append("T$_{fac}$")
            parameters.append("T$_{phot}$")

    return parameters, len(parameters)


def clr_to_vmr(param, clogr):
    """
        This function converts clr (centered log ratio) to vmr (volume mixing ratios) for all molecules
        considered in the fit and in the gas filling the rest of the atmosphere.
        The conversion is made using the equation exp(clr) / sum(exp(clr))
        for each molecule, ensuring that the sum of all volume mixing ratios is 1.

        Args:
            param (dict): A dictionary containing various parameters, including 'gas_fill' and 'fit_molecules'.
                'gas_fill' is the key for the molecule filling the rest of the atmosphere.
                'fit_molecules' is a list of the keys for the molecules considered in the fit.
                The dictionary will be updated in-place to include the calculated vmr's.
            clogr (dict): A dictionary with the clr's for each molecule (keys being molecule names).

        Returns:
            param (dict): The input dictionary updated with the calculated vmr's. Each vmr is added
                with a key of the form 'vmr_' + molecule name.
    """
    clogr[param['gas_fill']] = 0.0
    clr_arr = []
    for mol in param['fit_molecules']:
        clogr[param['gas_fill']] -= clogr[mol]
        clr_arr.append(clogr[mol])
    clr_arr.append(clogr[param['gas_fill']])

    vmr_arr = clr_inv(clr_arr)

    for i, mol in enumerate(param['fit_molecules']):
        param['vmr_' + mol] = vmr_arr[i]

    return param


def add_noise(param, data):
    """
        This function adds Gaussian noise to a given set of spectral data.
        It generates noise based on the provided signal-to-noise ratio (SNR)
        and applies it to the spectral data using the nested function gaussian_noise.

        Args:
            param (dict): A dictionary containing various parameters, including 'spectrum' and 'snr'.
                'spectrum' is a dictionary that contains 'wl', the wavelengths at which the spectrum is defined.
                'snr' is the signal-to-noise ratio for the data.
            data (numpy.ndarray): A 2D numpy array with the wavelengths in the first column and the corresponding
                spectral data in the second column.

        Returns:
            spec (numpy.ndarray): A 2D numpy array similar to 'data', but with Gaussian noise added to the spectral data.
                It has wavelengths in the first column, noisy spectral data in the second column, and the noise (standard deviation)
                in the third column.
    """

    def gaussian_noise(spec, no_less_zero=False):
        '''
        Adds gaussian noise with sigma=err to spectrum

        Parameters
        ----------
        spectrum : array-like
            planet spectrum or contrast.
        err : array-like
            error bars on spectrum for each point.

        Returns
        -------
        spec_with_error : array
            spectrum with gaussian noise added.
        '''

        spec = spec + 0.0
        for i in range(0, len(spec[:, 1])):
            point = np.random.normal(spec[i, 1], spec[i, 2])
            if no_less_zero:
                while point < 0.0:
                    point = np.random.normal(spec[i, 1], spec[i, 2])
            spec[i, 1] = point + 0.0
        return spec

    # if we are using noise file (like pandexo) then use these as the error
    if param['use_noise_file']:
        try:
            error_file = np.loadtxt(param['wkg_dir']+param['noise_file'],skiprows=1)
        except KeyError:
            print("With use_noise_file on, you must provide a noise file")
            sys.exit()
        wl_e = error_file[:,0]
        sp_e = error_file[:,1]
        err_e = error_file[:,3]
        # use spectres to get new errors for our wavelength grid
        _, err = spectres(param['spectrum']['wl'],wl_e,sp_e,err_e)
    # otherwise, calculate noise based on snr like usual
    else:
        err = np.full(len(param['spectrum']['wl']), (max(data[:, 1] / param['snr'])))
    spec = np.array([data[:, 0], data[:, 1], err]).T
    if param['gaussian_noise']:
        spec = gaussian_noise(spec, no_less_zero=True)

    return spec


def elapsed(t):
    """
        This function converts a given time duration in microseconds into a more readable format.
        It takes as input the duration in microseconds and outputs the time in days, hours, minutes,
        seconds, and milliseconds as appropriate. The function rounds the output at each level of conversion
        to the nearest integer.

        Args:
            t (float): A time duration in nanoseconds.

        Prints:
            A string that represents the time duration in a human-readable format. Depending on the length
            of the duration, this could be in the format of days, hours, minutes, seconds, and/or milliseconds.
    """
    milliseconds = round(t / 10 ** 6., 0)  # in milliseconds
    if milliseconds > 10 ** 3:
        seconds = int(milliseconds / 10 ** 3.)  # in seconds
        milliseconds = milliseconds - (seconds * (10 ** 3.))
        if seconds / 60. > 1:
            minutes = int(seconds / 60.)
            seconds = int(seconds - (minutes * 60.))
            if minutes / 60. > 1:
                hours = int(minutes / 60.)
                minutes = int(minutes - (hours * 60.))
                if hours / 24. > 1:
                    days = int(hours / 24.)
                    hours = int(hours - (days * 24.))
                    print('ExoTR runtime : ' + str(days) + ' days, ' + str(hours) + ' hours, ' + str(
                        minutes) + ' minutes, and ' + str(seconds) + ' seconds')
                else:
                    print('ExoTR runtime : ' + str(hours) + ' hours, ' + str(minutes) + ' minutes, and ' + str(
                        seconds) + ' seconds')
            else:
                print('ExoTR runtime : ' + str(minutes) + ' minutes and ' + str(seconds) + ' seconds')
        else:
            print('ExoTR runtime : ' + str(seconds) + ' seconds and ' + str(milliseconds) + ' milliseconds')
    else:
        print('ExoTR runtime : ' + str(milliseconds) + ' milliseconds')


def Ts_prior(param, Ts_phot_cube, Ts_het_cube=None, Ts_spot_cube=None, Ts_fac_cube=None):
    if param['Ts_err'] is None:
        Ts_phot_value = Ts_phot_cube * (param['Ts_phot_range'][1] - param['Ts_phot_range'][0]) + param['Ts_phot_range'][0]  # flat prior
    else:
        Ts_range = np.linspace(param['Ts_phot_range'][0], param['Ts_phot_range'][1], num=10000, endpoint=True)
        Ts_cdf = sp.stats.norm.cdf(Ts_range, param['Ts'], param['Ts_err'])
        Ts_cdf = np.array([0.0] + list(Ts_cdf) + [1.0])
        Ts_range = np.array([Ts_range[0]] + list(Ts_range) + [Ts_range[-1]])
        Ts_pri = interp1d(Ts_cdf, Ts_range)
        Ts_phot_value = Ts_pri(Ts_phot_cube)

    if param['stellar_activity_parameters'] == int(3):
        Ts_het_value = Ts_het_cube * (param['Ts_het_range'][1] - param['Ts_het_range'][0]) + param['Ts_het_range'][0]

        return Ts_het_value, Ts_phot_value

    elif param['stellar_activity_parameters'] == int(5):
        param['Ts_spot_range'] = [1500, Ts_phot_value]
        param['Ts_fac_range'] = [Ts_phot_value, 1.2 * param['Ts']]
        Ts_spot_value = Ts_spot_cube * (param['Ts_spot_range'][1] - param['Ts_spot_range'][0]) + param['Ts_spot_range'][0]
        Ts_fac_value = Ts_fac_cube * (param['Ts_fac_range'][1] - param['Ts_fac_range'][0]) + param['Ts_fac_range'][0]

        return Ts_spot_value, Ts_fac_value, Ts_phot_value


def readcross(fname):
    NTEMP = 20  # Number of temperature values (from Row 1)
    NPRESSURE = 10  # Number of pressure values (from Row 2)

    wave = []
    opac = []

    with open(fname, 'r') as fim:
        lines = fim.readlines()

    line_idx = 0  # Start from the first line

    # Read the temperature values from the first line
    temp_line = lines[line_idx].strip()
    temp_values = [float(x) for x in temp_line.split()]
    if len(temp_values) != NTEMP:
        raise ValueError(f"Expected {NTEMP} temperature values, got {len(temp_values)}")
    temp = temp_values
    line_idx += 1  # Move to the next line

    # Read the pressure values from the second line
    pres_line = lines[line_idx].strip()
    pres_values = [float(x) for x in pres_line.split()]
    if len(pres_values) != NPRESSURE:
        raise ValueError(f"Expected {NPRESSURE} pressure values, got {len(pres_values)}")
    pres = pres_values
    line_idx += 1  # Move to the next line

    # Now read the data blocks
    while line_idx < len(lines):
        # Read wave number (should be a line with 1 value)
        wave_line = lines[line_idx].strip()
        if not wave_line:
            line_idx += 1
            continue  # Skip empty lines
        wave_number = float(wave_line)
        wave.append(wave_number)
        line_idx += 1

        # Read NPRESSURE blocks of data
        opac_block = []
        for _ in range(NPRESSURE):
            if line_idx >= len(lines):
                raise ValueError("Unexpected end of file while reading data block")

            # Read initial opacity values
            data_line = lines[line_idx].strip()
            data_values = data_line.split()
            opacities = [float(val) for val in data_values[1:]]
            line_idx += 1

            # If opacities are split across multiple lines
            while len(opacities) < NTEMP:
                if line_idx >= len(lines):
                    raise ValueError("Unexpected end of file while reading opacities")
                data_line = lines[line_idx].strip()
                data_values = data_line.split()
                opacities.extend([float(val) for val in data_values])
                line_idx += 1

            if len(opacities) > NTEMP:
                opacities = opacities[:NTEMP]  # Ensure exact size
            opac_block.append(opacities)

        opac.append(opac_block)

    del lines

    # Convert lists to NumPy arrays for easier handling
    temp = np.array(temp)  # Shape: (NTEMP,)
    pres = np.array(pres)  # Shape: (NPRESSURE,)
    wave = np.array(wave)  # Shape: (number of wave numbers,)
    opac = np.array(opac)  # Shape: (number of wave numbers, NPRESSURE, NTEMP)

    temp = temp.reshape(1, -1)
    pres = pres.reshape(1, -1)
    wave = wave.reshape(1, -1)
    opac = opac.transpose(1, 2, 0)

    return temp, pres, wave, opac

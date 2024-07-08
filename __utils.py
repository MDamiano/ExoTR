from .__basics import *


def take_standard_parameters(pkg_dir):
    """
        This function reads the 'standard_parameters.dat' file located in the provided directory and processes
        its content into a Python dictionary. It ignores lines that start with '%' or are empty. It interprets
        the rest of the lines as key-value pairs separated by tabs.

        Special handling is included for 'mol', 'mol_vmr', and 'range_mol' keys. For these keys, the values are split
        into a list of items, with 'mol_vmr' and 'range_mol' values further converted into floats. If 'mol_vmr' is the key,
        additional keys in the format 'vmr_' + molecule are added to the dictionary.

        Boolean, None, and numeric values are automatically converted to their appropriate Python types.
        Strings are kept as strings.

        After processing the file, the function adds 'wkg_dir' with the current working directory and
        'supported_molecules' with a list of specific molecules to the dictionary.

        Args:
            pkg_dir (str): The directory in which the 'standard_parameters.dat' file is located.

        Returns:
            dict: A dictionary representing the parameters from the file, along with 'wkg_dir' and 'supported_molecules' entries.
    """
    parfile = 'standard_parameters.dat'
    with open(pkg_dir + parfile, 'r') as file:
        paramfile = file.readlines()
    param = {}
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

    param['wkg_dir'] = os.getcwd()
    param['supported_molecules'] = ['H2O', 'CH4', 'C2H2', 'C2H4', 'C2H6', 'NH3', 'HCN', 'H2S', 'SO2', 'CO', 'CO2', 'N2']

    for mol in param['supported_molecules']:
        param[mol + '_contribution'] = True
    param['cld_contribution'] = True
    param['CIA_contribution'] = True
    param['Rayleigh_contribution'] = True
    param['star_act_contribution'] = True
    param['haze_contribution'] = True

    mm = {'H': 1.00784, 'He': 4.002602, 'C': 12.0107, 'N': 14.0067, 'O': 15.9994, 'S': 32.065}
    mm['H2'] = 2. * mm['H']
    mm['H2O'] = (2. * mm['H']) + mm['O']
    mm['CH4'] = mm['C'] + (4. * mm['H'])
    mm['C2H2'] = (2. * mm['C']) + (2. * mm['H'])
    mm['C2H4'] = (2. * mm['C']) + (4. * mm['H'])
    mm['C2H6'] = (2. * mm['C']) + (6. * mm['H'])
    mm['NH3'] = mm['N'] + (3. * mm['H'])
    mm['HCN'] = mm['H'] + mm['C'] + mm['N']
    mm['H2S'] = (2. * mm['H']) + mm['S']
    mm['SO2'] = mm['S'] + (2. * mm['O'])
    mm['CO'] = mm['C'] + mm['O']
    mm['CO2'] = mm['C'] + (2. * mm['O'])
    mm['O2'] = 2. * mm['O']
    mm['O3'] = 3. * mm['O']
    mm['N2'] = 2. * mm['N']

    param['mm'] = mm

    param['formatted_labels'] = {}
    param['formatted_labels']['H2'] = "Log(H$_2$)"
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
        if mol == 'NH3':
            param['formatted_labels'][mol] = "Log(NH$_3$)"
        if mol == 'HCN':
            param['formatted_labels'][mol] = "Log(HCN)"
        if mol == 'H2S':
            param['formatted_labels'][mol] = "Log(H$_2$S)"
        if mol == 'SO2':
            param['formatted_labels'][mol] = "Log(SO$_2$)"
        if mol == 'CO':
            param['formatted_labels'][mol] = "Log(CO)"
        if mol == 'CO2':
            param['formatted_labels'][mol] = "Log(CO$_2$)"
        if mol == 'N2':
            param['formatted_labels'][mol] = "Log(N$_2$)"

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
        #print('Reading parfile: "' + parfile + '"')
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
            if param['fit_' + mol]:
                param['fit_molecules'].append(mol)
        if len(param['fit_molecules']) < 1 or 'H2O' not in param['fit_molecules']:
            param['vmr_H2O'] = 10. ** (-20.)

    if param['incl_star_activity']:
        if param['optimizer'] is None:
            try:
                param['Ts_phot'] += 0.0
            except KeyError:
                param['Ts_phot'] = param['Ts']

        if not param['fit_Ts_phot']:
            param['Ts_phot'] = param['Ts']

        if param['meta'] is None:
            param['meta'] = 0.0

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
    # planet
    if not param['fit_T']:
        try:
            param['Tp'] += 0.0
        except KeyError:
            t1 = ((param['Rs'] * const.R_sun) / (2. * param['major-a'] * const.au)) ** 0.5
            param['Tp'] = t1 * ((1 - 0.3) ** 0.25) * param['Ts']

    param['P_standard'] = 10. ** np.arange(-1.0, 12.1, step=0.1)
    param['P'] = 10. ** np.arange(-1.0, 8.1, step=0.1)

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
    vmr_M, MMM = 0.0, 0.0
    for mol in param['fit_molecules']:
        if mol != 'H2O':
            vmr_M += param['vmr_' + mol]
            MMM += param['vmr_' + mol] * param['mm'][mol]
    param['vmr_' + param['gas_fill']] = np.zeros(len(param['P']))
    param['mean_mol_weight'] = MMM * np.ones(len(param['P']))

    for i in range(0, len(param['P'])):
        param['vmr_' + param['gas_fill']][i] = (10. ** 0.0) - (vmr_M + param['vmr_H2O'][i])
        param['mean_mol_weight'][i] += ((param['vmr_H2O'][i] * param['mm']['H2O']) + (param['vmr_' + param['gas_fill']][i] * param['mm'][param['gas_fill']]))

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
            except IndexError:
                param['spectrum']['wl'] = spectrum
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


def particlesizef(g, T, P, M, MM, KE, deltaP):
    """
        Calculate particle size in exoplanet atmospheres using various input parameters.

        Parameters:
        g (float): Gravity in SI units (m/s^2).
        T (float): Temperature in Kelvin.
        P (float): Pressure in Pascal.
        M (float): Mean molecular mass of the atmosphere in g/mol.
        MM (float): Molecular mass of the condensable species in g/mol.
        KE (float): Eddy diffusion coefficient in m^2/s.
        deltaP (float): Difference between partial pressure and saturation vapor pressure in Pa.

        Returns:
        r0 (float): Mode radius of the droplet size distribution.
        r1 (float): Radius of the droplet for the lower quartile of the distribution.
        r2 (float): Radius of the droplet for the upper quartile of the distribution.
        VP (float): Volume of the droplet in cm^3.

        The function uses a number of hard-coded constants and assumptions, such as
        the density of the condensed material of water (1000 kg/m^3), the accomodation
        factor (1), and the width of the log-normal size distribution (2). The output
        particle size is in microns, and the volume is in cm^3.
    """
    # Calculate particle size in exoplanet atmospheres

    # input
    # g in SI
    # T in K
    # P in Pa
    # M: mean molecular mass of the atmosphere; in g / mol
    # MM: molecular mass of the condensable species; in g / mol
    # KE: Eddy diffusion coefficient; in m2 s - 1
    # deltaP: difference between partial pressure and saturation vapor
    # pressure, in Pa

    # assume
    # density of condensed material of water 1000 kg / m3
    # accomodation factor of unity
    # sig = 2

    # output particle size in micron, and volumn in cm ^ 3

    # Derived parameters
    H = (const.k_B.value * T) / M / const.u.value / g
    u = KE / H
    mu = ((8.76E-6 * (293.85 + 72)) / (293.85 + 72)) * ((T / 293.85) ** 1.5)  # SI
    lamb = (2. * mu) / P / ((8 * M * 1.0E-3 / math.pi / 8.314472 / T) ** 0.5)  # m
    # KK = 4 * KB * T / 3. / mu
    deltan = deltaP / const.k_B.value / T

    # droplet
    rho = 1.0E+3  # kg m-3
    acc = 1.0

    # mass diffusion coefficient
    D = 0.12E-4

    # Particle Size and Number
    Cc, fa = 1, 1
    Cc1, fa1 = 2, 2
    sig = 2

    check = 0
    while (abs(Cc1 - Cc) + abs(fa1 - fa)) > 0.001:
        Cc = Cc1
        fa = fa1
        cc = -((48. * math.pi ** 2.) ** (1. / 3.)) * D * MM * const.u.value * fa * deltan / rho * np.exp(- np.log(sig) ** 2.)  # effective condensation coefficient D
        aa = rho * g / mu / ((162. * math.pi ** 2.) ** (1. / 3.)) / H * Cc * np.exp(- np.log(sig) ** 2.)
        bb = -u / H

        V = ((-bb + np.sqrt((bb ** 2.) - (4. * aa * cc))) / 2. / aa) ** (3. / 2.)
        d1 = ((6. * V / math.pi) ** (1. / 3.)) * np.exp(- np.log(sig) ** 2.)

        kn = lamb / d1
        Cc1 = 1. + kn * (1.257 + 0.4 * np.exp(- 1.1 / kn))
        fa1 = (1. + kn) / (1. + 2. * kn * (1. + kn) / acc)

        Vs = V + 0.0
        check += 1
        if check > 1e4:
            break

    r0 = (3. * Vs / 4. / math.pi) ** (1. / 3.) * np.exp(- 1.5 * np.log(sig) ** 2.) * 1.0E+6
    r1 = (3. * Vs / 4. / math.pi) ** (1. / 3.) * np.exp(- 1.0 * np.log(sig) ** 2.) * 1.0E+6
    r2 = (3. * Vs / 4. / math.pi) ** (1. / 3.) * np.exp(- 0.5 * np.log(sig) ** 2.) * 1.0E+6
    VP = Vs * 1.0E+6

    return r0, r1, r2, VP


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

    param['vmr_H2O'] = gaussian_filter1d(watermix, 1, mode='nearest')
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
        param['ptopw_range'] = [0.0, 8.0]  # Top pressure H2O
        param['dcldw_range'] = [0.0, 8.5]  # Depth H2O cloud
        param['crh2o_range'] = [-12.0, 0.0]  # Condensation Ratio H2O
    if param['fit_amn_cld']:
        param['ptopa_range'] = [0.0, 8.0]  # Top pressure NH3
        param['dclda_range'] = [0.0, 8.5]  # Depth HN3 cloud
        param['crnh3_range'] = [-12.0, 0.0]  # Condensation Ratio NH3
    if param['fit_gen_cld']:
        param['ptop_range'] = [0.0, 9.0]  # Top pressure
    if param['incl_haze']:
        param['dhaze_range'] = [-3, 2] # diameter of haze particle
        param['vmrhaze_range'] = [-10, -1]
    if param['fit_T'] and param['Tp'] is None:
        param['tp_range'] = [100.0, 2000.0]  # Atmospheric equilibrium temperature
    elif param['fit_T'] and param['Tp'] is not None:
        param['tp_range'] = [min([100.0, param['Tp'] - 500.0]), param['Tp'] + 500.0]  # Atmospheric equilibrium temperature
    else:
        pass
    if param['fit_Rp']:
        param['rp_range'] = [param['Rp'] * 0.5, param['Rp'] * 2.0]  # Planetary Radius
    if param['fit_offset']:
        for i in range(0, param['n_offsets']):
            param['off' + str(i + 1) + '_range'] = [-param['offset_range'] / 1e6, param['offset_range'] / 1e6]
    if param['incl_star_activity']:
        if param['fit_het_frac']:
            param['delta_range'] = [0.0, 0.5]  # Fraction of the star surface affected by activity
        if param['fit_Ts_phot']:
            param['Ts_phot_range'] = [1000, 10000]

    return param


def take_star_spectrum(param, Ts, meta=None):
    # PHOENIX STELLAR SPECTRA FROM : http://svo2.cab.inta-csic.es/theory/newov2/index.php

    def resize(fl):
        idx_start = find_nearest(fl[:, 0] * (10. ** -4.), param['min_wl'])
        idx_stop = find_nearest(fl[:, 0] * (10. ** -4.), param['max_wl'])
        if param['spectrum']['bins']:
            return custom_spectral_binning(np.array([param['spectrum']['wl_low'], param['spectrum']['wl_high'], param['spectrum']['wl']]).T, fl[idx_start:idx_stop, 0] * (10. ** -4.), fl[idx_start:idx_stop, 1], bins=True)
        else:
            return custom_spectral_binning(param['spectrum']['wl'], fl[idx_start:idx_stop, 0] * (10. ** -4.), fl[idx_start:idx_stop, 1])

    directory = param['pkg_dir'] + 'PHOENIX_models/'
    t_star = round(float(Ts), -2)
    interp_Ts = True
    if t_star - Ts < 0:
        t_star1 = t_star + 0.0
        t_star2 = t_star + 100
    elif t_star - Ts > 0:
        t_star1 = t_star - 100
        t_star2 = t_star + 0.0
    else:
        t_star = str(t_star / 100)
        interp_Ts = False

    if interp_Ts:
        t_star2_ratio = 1.0 - ((t_star2 - Ts) / (t_star2 - t_star1))
        t_star1_ratio = 1.0 - t_star2_ratio
        if t_star1 < 1000:
            t_star1 = '00' + str(int(t_star1 / 100))
        elif t_star1 < 10000:
            t_star1 = '0' + str(int(t_star1 / 100))
        else:
            t_star1 = str(int(t_star1 / 100))
        if t_star2 < 1000:
            t_star2 = '00' + str(int(t_star2 / 100))
        elif t_star2 < 10000:
            t_star2 = '0' + str(int(t_star2 / 100))
        else:
            t_star2 = str(int(t_star2 / 100))

    try:
        param['Loggs'] += 0.0
    except KeyError:
        param['Loggs'] = np.log10(const.G.value * (param['Ms'] * const.M_sun.value) / ((param['Rs'] * const.R_sun.value) ** 2.)) + 2

    loggs = round(param['Loggs'] * 2.) / 2.
    interp_loggs = True
    if loggs - param['Loggs'] < 0:
        loggs1 = loggs + 0.0
        loggs2 = loggs + 0.5
    elif loggs - param['Loggs'] > 0:
        loggs1 = loggs - 0.5
        loggs2 = loggs + 0.0
    else:
        loggs = '-' + str(loggs)
        interp_loggs = False

    if interp_loggs:
        loggs2_ratio = 1.0 - ((loggs2 - param['Loggs']) / (loggs2 - loggs1))
        loggs1_ratio = 1.0 - loggs2_ratio
        loggs1 = '-' + str(loggs1)
        loggs2 = '-' + str(loggs2)

    if meta is None:
        s_meta = '-' + str(0.0)
        interp_meta = False
    else:
        s_meta = round(meta * 2.) / 2.
        interp_meta = True
        if s_meta - meta < 0:
            s_meta1 = s_meta + 0.0
            s_meta2 = s_meta + 0.5
        elif s_meta - meta > 0:
            s_meta1 = s_meta - 0.5
            s_meta2 = s_meta + 0.0
        else:
            interp_meta = False
            s_meta = '-' + str(0.0)

    if interp_meta:
        s_meta2_ratio = 1.0 - ((s_meta2 - meta) / (s_meta2 - s_meta1))
        s_meta1_ratio = 1.0 - s_meta2_ratio
        if s_meta1 <= 0:
            s_meta1 = '-' + str(abs(s_meta1))
        else:
            s_meta1 = '+' + str(s_meta1)
        if s_meta2 <= 0:
            s_meta2 = '-' + str(abs(s_meta2))
        else:
            s_meta2 = '+' + str(s_meta2)

    s_files = {}
    if interp_Ts and interp_loggs and interp_meta:
        if float(t_star1) < 26:
            s_meta1 = '-' + str(0.0)
        if float(t_star2) < 26:
            s_meta2 = '-' + str(0.0)
        s_files['A'] = np.genfromtxt(directory + 'lte' + t_star1 + loggs1 + s_meta1 + '.txt', skip_header=6)
        s_files['B'] = np.genfromtxt(directory + 'lte' + t_star2 + loggs1 + s_meta1 + '.txt', skip_header=6)
        s_files['C'] = np.genfromtxt(directory + 'lte' + t_star1 + loggs2 + s_meta1 + '.txt', skip_header=6)
        s_files['D'] = np.genfromtxt(directory + 'lte' + t_star2 + loggs2 + s_meta1 + '.txt', skip_header=6)
        s_files['E'] = np.genfromtxt(directory + 'lte' + t_star1 + loggs1 + s_meta2 + '.txt', skip_header=6)
        s_files['F'] = np.genfromtxt(directory + 'lte' + t_star2 + loggs1 + s_meta2 + '.txt', skip_header=6)
        s_files['G'] = np.genfromtxt(directory + 'lte' + t_star1 + loggs2 + s_meta2 + '.txt', skip_header=6)
        s_files['H'] = np.genfromtxt(directory + 'lte' + t_star2 + loggs2 + s_meta2 + '.txt', skip_header=6)
    elif interp_Ts and interp_loggs and not interp_meta:
        s_files['A'] = np.genfromtxt(directory + 'lte' + t_star1 + loggs1 + s_meta + '.txt', skip_header=6)
        s_files['B'] = np.genfromtxt(directory + 'lte' + t_star2 + loggs1 + s_meta + '.txt', skip_header=6)
        s_files['C'] = np.genfromtxt(directory + 'lte' + t_star1 + loggs2 + s_meta + '.txt', skip_header=6)
        s_files['D'] = np.genfromtxt(directory + 'lte' + t_star2 + loggs2 + s_meta + '.txt', skip_header=6)
    elif interp_Ts and interp_meta and not interp_loggs:
        if float(t_star1) < 26:
            s_meta1 = '-' + str(0.0)
        if float(t_star2) < 26:
            s_meta2 = '-' + str(0.0)
        s_files['A'] = np.genfromtxt(directory + 'lte' + t_star1 + loggs + s_meta1 + '.txt', skip_header=6)
        s_files['B'] = np.genfromtxt(directory + 'lte' + t_star2 + loggs + s_meta1 + '.txt', skip_header=6)
        s_files['C'] = np.genfromtxt(directory + 'lte' + t_star1 + loggs + s_meta2 + '.txt', skip_header=6)
        s_files['D'] = np.genfromtxt(directory + 'lte' + t_star2 + loggs + s_meta2 + '.txt', skip_header=6)
    elif interp_loggs and interp_meta and not interp_Ts:
        s_files['A'] = np.genfromtxt(directory + 'lte' + t_star + loggs1 + s_meta1 + '.txt', skip_header=6)
        s_files['B'] = np.genfromtxt(directory + 'lte' + t_star + loggs2 + s_meta1 + '.txt', skip_header=6)
        s_files['C'] = np.genfromtxt(directory + 'lte' + t_star + loggs1 + s_meta2 + '.txt', skip_header=6)
        s_files['D'] = np.genfromtxt(directory + 'lte' + t_star + loggs2 + s_meta2 + '.txt', skip_header=6)
    elif interp_Ts and not interp_loggs and not interp_meta:
        s_files['A'] = np.genfromtxt(directory + 'lte' + t_star1 + loggs + s_meta + '.txt', skip_header=6)
        s_files['B'] = np.genfromtxt(directory + 'lte' + t_star2 + loggs + s_meta + '.txt', skip_header=6)
    elif interp_loggs and not interp_Ts and not interp_meta:
        s_files['A'] = np.genfromtxt(directory + 'lte' + t_star + loggs1 + s_meta + '.txt', skip_header=6)
        s_files['B'] = np.genfromtxt(directory + 'lte' + t_star + loggs2 + s_meta + '.txt', skip_header=6)
    elif interp_meta and not interp_Ts and not interp_loggs:
        s_files['A'] = np.genfromtxt(directory + 'lte' + t_star + loggs + s_meta1 + '.txt', skip_header=6)
        s_files['B'] = np.genfromtxt(directory + 'lte' + t_star + loggs + s_meta2 + '.txt', skip_header=6)
    else:
        s_files['A'] = np.genfromtxt(directory + 'lte' + t_star + loggs + s_meta + '.txt', skip_header=6)

    # plt.plot(s_files['D'][:, 0] * (10. ** -4.), s_files['D'][:, 1])
    # plt.plot(s_files['B'][:, 0] * (10. ** -4.), s_files['B'][:, 1])
    # plt.plot(s_files['C'][:, 0] * (10. ** -4.), s_files['C'][:, 1])
    # plt.plot(s_files['A'][:, 0] * (10. ** -4.), s_files['A'][:, 1])

    for s_keys in s_files.keys():
        s_files[s_keys] = resize(s_files[s_keys])

    if interp_Ts and interp_loggs and interp_meta:
        sp = s_meta1_ratio * (loggs1_ratio * (t_star1_ratio * s_files['A'] + t_star2_ratio * s_files['B']) + loggs2_ratio * (t_star1_ratio * s_files['C'] + t_star2_ratio * s_files['D'])) + \
             s_meta2_ratio * (loggs1_ratio * (t_star1_ratio * s_files['E'] + t_star2_ratio * s_files['F']) + loggs2_ratio * (t_star1_ratio * s_files['G'] + t_star2_ratio * s_files['H']))
    elif interp_Ts and interp_loggs and not interp_meta:
        sp = loggs1_ratio * (t_star1_ratio * s_files['A'] + t_star2_ratio * s_files['B']) + loggs2_ratio * (t_star1_ratio * s_files['C'] + t_star2_ratio * s_files['D'])
    elif interp_Ts and interp_meta and not interp_loggs:
        sp = s_meta1_ratio * (t_star1_ratio * s_files['A'] + t_star2_ratio * s_files['B']) + s_meta2_ratio * (t_star1_ratio * s_files['C'] + t_star2_ratio * s_files['D'])
    elif interp_loggs and interp_meta and not interp_Ts:
        sp = s_meta1_ratio * (loggs1_ratio * s_files['A'] + loggs2_ratio * s_files['B']) + s_meta2_ratio * (loggs1_ratio * s_files['C'] + loggs2_ratio * s_files['D'])
    elif interp_Ts and not interp_loggs and not interp_meta:
        sp = t_star1_ratio * s_files['A'] + t_star2_ratio * s_files['B']
    elif interp_loggs and not interp_Ts and not interp_meta:
        sp = loggs1_ratio * s_files['A'] + loggs2_ratio * s_files['B']
    elif interp_meta and not interp_Ts and not interp_loggs:
        sp = s_meta1_ratio * s_files['A'] + s_meta2_ratio * s_files['B']
    else:
        sp = s_files['A'] + 0.0

    # plt.plot(param['spectrum']['wl'], s_files['D'])
    # plt.plot(param['spectrum']['wl'], s_files['B'])
    # plt.plot(param['spectrum']['wl'], s_files['C'])
    # plt.plot(param['spectrum']['wl'], s_files['A'])
    #
    # plt.plot(param['spectrum']['wl'], sp, linestyle='--')
    # plt.xlim([0.5, 5.5])
    # plt.ylim([0.0, 500000])
    # plt.show()

    return sp


def pre_load_variables(param):
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
    data = scipy.io.loadmat(param['pkg_dir'] + 'Data/opac/opac_052024.mat')
    opac_data_keys = []
    for i in data.keys():
        if i != '__header__' and i != '__globals__' and i != '__version__':
            param[i] = np.array(data[i]) + 0.0
            opac_data_keys.append(i)

    del data

    param['opac_data_keys'] = opac_data_keys
    param['opacaer_no_cloud'] = param['opacaerh2o'] * 0.0

    return param


def custom_spectral_binning(x, wl, model, err=None, bins=False):
    """
        Applies custom spectral binning to a model spectrum based on the given wavelength grid.

        Parameters:
        x (np.array): An array of wavelengths at which the model is to be binned.
        wl (np.array): An array of the original wavelengths of the model.
        model (np.array): An array representing the spectrum of the model.
        err (np.array, optional): An array representing the error associated with each point in the spectrum.
                                  Default is None.

        Returns:
        np.array: The binned model spectrum.
                  If `err` is not None, also returns an array of the binned errors.

        The function applies a custom binning scheme to the input model spectrum based on the input wavelength grid (`x`).
        If the error array (`err`) is provided, the function also applies the same binning scheme to the errors and returns
        the binned error array.
    """
    binned_mod = []
    if err is not None:
        binned_er = []

    if not bins:
        y = np.roll(x, 1) + 0.0
        dx = (x - y)[1:]
        limits = []

        i, intermed = 0, 0
        while i in range(0, len(dx)):
            if dx[i] == dx[0]:
                lim = (dx[i] / 2., dx[i] / 2.)
                limits.append(lim)
            elif dx[i] > 2 * np.median(dx[intermed:i]):
                lim = (dx[i - 1] / 2., dx[i - 1] / 2.)
                limits.append(lim)
                i += 1
                if i != len(dx):
                    lim = (dx[i] / 2., dx[i] / 2.)
                    limits.append(lim)
                    intermed = i + 1
                else:
                    break
            else:
                lim = (dx[i - 1] / 2., dx[i] / 2.)
                limits.append(lim)
            i += 1
        limits.append((dx[i - 1] / 2., dx[i - 1] / 2.))

        for i in range(0, len(x)):
            yy = np.array(model[np.where((wl > x[i] - limits[i][0]) & (wl < x[i] + limits[i][1]))[0]])
            binned_mod.append(np.mean(yy))
            if err is not None:
                er = np.array(err[np.where((wl > x[i] - limits[i][0]) & (wl < x[i] + limits[i][1]))[0]])
                binned_er.append(np.sqrt(np.sum(er ** 2.)) / len(er))
    else:
        for i in range(0, len(x[:, 0])):
            yy = np.array(model[np.where((wl > x[i, 0]) & (wl < x[i, 1]))[0]])
            binned_mod.append(np.mean(yy))
            if err is not None:
                er = np.array(err[np.where((wl > x[i, 0]) & (wl < x[i, 1]))[0]])
                binned_er.append(np.sqrt(np.sum(er ** 2.)) / len(er))

    if err is None:
        return np.array(binned_mod)
    else:
        return np.array(binned_mod), np.array(binned_er)


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
    if param['fit_T']:
        parameters.append("Log(T$_p$)")
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
    if param['incl_haze']:
        parameters.append("Log(diam$_{haze}$)")
        parameters.append("Log(vmr$_{haze}$)")
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
    if param['fit_NH3']:
        parameters.append("clr(NH$_3$)")
    if param['fit_HCN']:
        parameters.append("clr(HCN)")
    if param['fit_H2S']:
        parameters.append("clr(H$_2$S)")
    if param['fit_SO2']:
        parameters.append("clr(SO$_2$)")
    if param['fit_CO']:
        parameters.append("clr(CO)")
    if param['fit_CO2']:
        parameters.append("clr(CO$_2$)")
    if param['fit_N2']:
        parameters.append("clr(N$_2$)")
    if param['incl_star_activity']:
        if param['fit_het_frac']:
            parameters.append("$\delta$")
        if param['fit_Ts_het']:
            parameters.append("T$_{het}$")
        if param['fit_Ts_phot']:
            parameters.append("T$_{phot}$")

    return parameters, len(parameters)


def clr_to_vmr(param, clr):
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
            clr (dict): A dictionary with the clr's for each molecule (keys being molecule names).

        Returns:
            param (dict): The input dictionary updated with the calculated vmr's. Each vmr is added
                with a key of the form 'vmr_' + molecule name.
    """
    clr[param['gas_fill']] = 0.0
    sumb = 0.0
    for mol in param['fit_molecules']:
        clr[param['gas_fill']] -= clr[mol]
        sumb += np.exp(clr[mol])

    sumb += np.exp(clr[param['gas_fill']])

    for mol in param['fit_molecules']:
        param['vmr_' + mol] = np.exp(clr[mol])/sumb
    param['vmr_' + param['gas_fill']] = np.exp(clr[param['gas_fill']]) / sumb

    return param


def plot_chemistry(param, cube, solutions=None):
    """
        Plots the vertical distribution of molecular volume mixing ratios and mean molecular weight.

        Parameters:
        param (dict): Dictionary of parameters containing information about the model's state and parameters.
        cube (np.array): An array with model parameters, to be updated in `param`.
        solutions (int, optional): Identifier for a specific solution to be labeled in the plot's file name.
                                   Default is None.

        This function generates two plots:
        1. Vertical distribution of molecular volume mixing ratios (VMRs) for different species. The y-axis represents pressure,
           and the x-axis represents VMR. Water, ammonia, and optically thick clouds are highlighted with dashed lines, if present.
        2. Vertical distribution of mean molecular weight (MMM). The y-axis represents pressure,
           and the x-axis represents MMM. Water, ammonia, and optically thick clouds are highlighted with dashed lines, if present.

        Both plots are saved to the directory specified in `param['out_dir']`, with filenames based on whether `solutions` is provided.
    """
    clr = {}
    par = 0
    if param['fit_offset']:
        for i in range(0, param['n_offsets']):
            param['offset' + str(i + 1)] = cube[par]  # offset between datasets
            par += 1
    if param['fit_Rp']:
        param['Rp'] = cube[par]
        par += 1
    if param['fit_T']:
        param['Tp'] = cube[par]
        par += 1
    if param['fit_wtr_cld']:
        param['Pw_top'] = (10. ** cube[par])
        param['cldw_depth'] = (10. ** cube[par + 1])
        param['CR_H2O'] = (10. ** cube[par + 2])
        par += 3
    if param['fit_amn_cld']:
        param['Pa_top'] = (10. ** cube[par])
        param['clda_depth'] = (10. ** cube[par + 1])
        param['CR_NH3'] = (10. ** cube[par + 2])
        par += 3
    if param['fit_gen_cld']:
        param['P_top'] = (10. ** cube[par])
        par += 1

    for mol in param['fit_molecules']:
        clr[mol] = cube[par]  # Molecules
        par += 1

    param = clr_to_vmr(param, clr)
    param = cloud_pos(param)
    param = calc_mean_mol_mass(param)

    fig, ax = plt.subplots()

    for mol in param['fit_molecules']:
        if mol != 'H2O':
            print(str(mol) + ' -> Top: ' + str((param['vmr_' + mol] * np.ones(len(param['P'])))[0]) + ', Bottom: ' + str((param['vmr_' + mol] * np.ones(len(param['P'])))[-1]))
            ax.loglog(param['vmr_' + mol] * np.ones(len(param['P'])), param['P'], label=mol)
        else:
            print(str(mol) + ' -> Top: ' + str(param['vmr_' + mol][0]) + ', Bottom: ' + str(param['vmr_' + mol][-1]))
            ax.loglog(param['vmr_' + mol], param['P'], label=mol)
    print(str(param['gas_fill']) + ' -> Top: ' + str((param['vmr_' + param['gas_fill']] * np.ones(len(param['P'])))[0]) + ', Bottom: ' + str((param['vmr_' + param['gas_fill']] * np.ones(len(param['P'])))[-1]))
    ax.loglog(param['vmr_' + param['gas_fill']] * np.ones(len(param['P'])), param['P'], label=param['gas_fill'])

    ax.set_ylim((param['P'][0], 10. ** 8.))
    ax.set_xlim((1e-18, 1.5))
    plt.gca().invert_yaxis()
    ax.set_xlabel('Molecular VMR')
    ax.set_ylabel('Pressure [Pa]')

    def pa_to_bar(y):
        return y / (10. ** 5.)

    def bar_to_pa(y):
        return y * (10. ** 5.)

    if param['fit_wtr_cld']:
        wtr2 = param['vmr_H2O'] - np.roll(param['vmr_H2O'], 1)
        wtr2[0] = 0.0
        plt.hlines(param['P'][min(np.where(wtr2 != 0.0)[0]) - 1], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black', label='H$_2$O cloud')
        plt.hlines(param['P'][max(np.where(wtr2 != 0.0)[0])], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black')
    if param['fit_amn_cld']:
        wtr2 = param['vmr_NH3'] - np.roll(param['vmr_NH3'], 1)
        wtr2[0] = 0.0
        plt.hlines(param['P'][min(np.where(wtr2 != 0.0)[0]) - 1], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black', label='NH$_3$ cloud')
        plt.hlines(param['P'][max(np.where(wtr2 != 0.0)[0])], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black')
    if param['fit_gen_cld']:
        indx = find_nearest(param['P'], param['P_top'])
        plt.hlines(param['P'][indx], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black', label='Optically thick cloud')

    ax.yaxis.set_ticks(10. ** np.arange(np.log10(param['P'][0]), 9, 1))
    secax_y = ax.secondary_yaxis('right', functions=(pa_to_bar, bar_to_pa))
    secax_y.set_ylabel('Pressure [bar]')

    ax.legend(loc='lower left')
    if solutions is None:
        plt.savefig(param['out_dir'] + 'Chemistry.pdf')
        plt.savefig(param['out_dir'] + 'Chemistry.png')
    else:
        plt.savefig(param['out_dir'] + 'Chemistry (solution ' + str(solutions) + ').pdf')
        plt.savefig(param['out_dir'] + 'Chemistry (solution ' + str(solutions) + ').png')
    plt.close()

    fig, ax = plt.subplots()
    ax.semilogy(param['mean_mol_weight'], param['P'])
    ax.set_ylim((param['P'][0], 10. ** 8.))
    plt.gca().invert_yaxis()
    ax.set_xlabel('Mean molecular weight')
    ax.set_ylabel('Pressure [Pa]')

    ax.set_xlim((ax.get_xlim()[0], ax.get_xlim()[1]))

    if param['fit_wtr_cld']:
        wtr2 = param['vmr_H2O'] - np.roll(param['vmr_H2O'], 1)
        wtr2[0] = 0.0
        plt.hlines(param['P'][min(np.where(wtr2 != 0.0)[0]) - 1], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black', label='H$_2$O cloud')
        plt.hlines(param['P'][max(np.where(wtr2 != 0.0)[0])], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black')
    if param['fit_amn_cld']:
        wtr2 = param['vmr_NH3'] - np.roll(param['vmr_NH3'], 1)
        wtr2[0] = 0.0
        plt.hlines(param['P'][min(np.where(wtr2 != 0.0)[0]) - 1], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black', label='NH$_3$ cloud')
        plt.hlines(param['P'][max(np.where(wtr2 != 0.0)[0])], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black')
    if param['fit_gen_cld']:
        indx = find_nearest(param['P'], param['P_top'])
        plt.hlines(param['P'][indx], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black', label='Optically thick cloud')

    ax.yaxis.set_ticks(10. ** np.arange(np.log10(param['P'][0]), 9, 1))
    secax_y = ax.secondary_yaxis('right', functions=(pa_to_bar, bar_to_pa))
    secax_y.set_ylabel('Pressure [bar]')
    if solutions is None:
        plt.savefig(param['out_dir'] + 'MMM.pdf')
        plt.savefig(param['out_dir'] + 'MMM.png')
    else:
        plt.savefig(param['out_dir'] + 'MMM (solution ' + str(solutions) + ').pdf')
        plt.savefig(param['out_dir'] + 'MMM (solution ' + str(solutions) + ').png')
    plt.close()


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

    err = np.full(len(param['spectrum']['wl']), (max(data[:, 1] / param['snr'])))

    spec = np.array([data[:, 0], data[:, 1], err]).T
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
                    print('Runtime : ' + str(days) + ' days, ' + str(hours) + ' hours, ' + str(
                        minutes) + ' minutes, and ' + str(seconds) + ' seconds')
                else:
                    print('Runtime : ' + str(hours) + ' hours, ' + str(minutes) + ' minutes, and ' + str(
                        seconds) + ' seconds')
            else:
                print('Runtime : ' + str(minutes) + ' minutes and ' + str(seconds) + ' seconds')
        else:
            print('Runtime : ' + str(seconds) + ' seconds and ' + str(milliseconds) + ' milliseconds')
    else:
        print('Runtime : ' + str(milliseconds) + ' milliseconds')


def Ts_prior(param, Ts_het_cube, Ts_phot_cube=None):
    if param['fit_Ts_phot']:
        if param['Ts_err'] is None:
            Ts_phot_value = Ts_phot_cube * (param['Ts_phot_range'][1] - param['Ts_phot_range'][0]) + param['Ts_phot_range'][0]  # ignorant prior
        else:
            Ts_range = np.linspace(param['Ts_phot_range'][0], param['Ts_phot_range'][1], num=10000)
            Ts_cdf = sp.stats.norm.cdf(Ts_range, param['Ts'], param['Ts_err'])
            Ts_pri = interp1d(Ts_cdf, Ts_range)
            Ts_phot_value = Ts_pri(Ts_phot_cube)
    else:
        Ts_phot_value = param['Ts_phot'] + 0.0

    param['Ts_het_range'] = [0.5 * Ts_phot_value, 1.2 * Ts_phot_value]
    Ts_het_value = Ts_het_cube * (param['Ts_het_range'][1] - param['Ts_het_range'][0]) + param['Ts_het_range'][0]  # ignorant prior

    return Ts_het_value, Ts_phot_value

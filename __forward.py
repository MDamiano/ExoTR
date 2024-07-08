from .__utils import *
import copy


class FORWARD_MODEL:
    """
        This class is used to generate a forward model for atmospheric studies of a given planet.

        Attributes:
            param (dict): The dictionary object containing all the necessary parameters for the atmospheric modeling.
            package_dir (str): The directory path to the required data files for the model.
            working_dir (str): The working directory path where the computations will be performed. If not provided,
                               it defaults to the current working directory.
    """
    def __init__(self, param):
        self.param = param
        self.package_dir = param['pkg_dir']
        try:
            self.working_dir = param['wkg_dir']
        except KeyError:
            self.working_dir = os.getcwd()

    def atmospheric_structure(self):
        """
            This function calculates and returns the atmospheric structure of a given planet by using the parameters provided
            in the 'param' object of the class. The function considers factors like planet radius, planet mass, star radius,
            atmospheric profile, molecular opacities, continuum absorption and scattering (CIA), Rayleigh scattering, cloud
            contributions and others.

            It uses these parameters to calculate the wavelength-dependent atmospheric opacity, the apparent size of the planet,
            and the transit depth, which refers to the fraction of the star's light blocked by the planet.

            Returns:
                tuple: A tuple containing the wavelengths of light (in micrometers) and the corresponding transit depth values.

            Raises:
                KeyError: If 'Mp' (Mass of Planet) is not provided in self.param.
                SystemExit: If KeyError is raised, the function stops execution.
        """
        R0 = self.param['Rp'] * const.R_earth.value  # Radius of Planet in meters
        try:
            M0 = self.param['Mp'] * const.M_earth.value  # Mass of Planet in kilograms
        except KeyError:
            print('ERROR - This version of the code does not support unknown planetary mass, please specify.')
            sys.exit()
        RS = self.param['Rs'] * const.R_sun.value  # Radius of Star in meters

        # Atmospheric Profile
        P = np.array([self.param['P'][::-1]]).T
        T = np.ones_like(P) * self.param['Tp']
        fH2O = np.array([self.param['vmr_H2O'][::-1]]).T
        mu = np.array([self.param['mean_mol_weight'][::-1]]).T
        n = len(P[:, 0])
        Z, g = np.zeros_like(P), np.zeros_like(P)

        for i in range(1, n):
            g[i - 1, 0] = 6.674E-11 * M0 / ((R0 + Z[i - 1, 0]) ** 2.0)
            H = 1.38E-23 * (T[i - 1, 0] + T[i, 0]) / (mu[i - 1, 0] + mu[i, 0]) / g[i - 1, 0] / 1.66E-27
            Z[i, 0] = Z[i - 1, 0] + np.log(P[i - 1, 0] / P[i, 0]) * H

        g[-1] = g[-2] - (g[-3] - g[-2])

        g = np.array([g[:-1, 0]]).T
        PL = np.array([np.sqrt(P[0:n - 1, 0] * P[1:n, 0])]).T
        TL = np.array([(T[0:n - 1, 0] + T[1:n, 0]) / 2.0]).T
        mul = np.array([(mu[0:n - 1, 0] + mu[1:n, 0]) / 2.0]).T
        NL = PL / 1.38E-23 / TL
        N = {'H2O': NL * np.sqrt(fH2O[0:n - 1] * fH2O[1:n])}
        for mol in self.param['fit_molecules']:
            if mol != 'H2O':
                N[mol] = NL * self.param['vmr_' + mol]
        ffill = np.array([self.param['vmr_' + self.param['gas_fill']][::-1]]).T
        N[self.param['gas_fill']] = NL * np.sqrt(ffill[0:n - 1] * ffill[1:n])

        I1 = np.array([np.ones(len(TL))]).T
        I2 = np.ones(len(self.param['opacw'][0]))

        # Get molecular opacities
        S = {}
        if len(self.param['fit_molecules']) < 1 or 'H2O' not in self.param['fit_molecules']:
            S['H2O'] = (interpn((self.param['opacp'][0], self.param['opact'][0], self.param['opacw'][0]), self.param['opach2o'], np.array([PL * I2, TL * I2, I1 * self.param['opacw'][0]]).T)).T
        for mol in self.param['fit_molecules'] + [self.param['gas_fill']]:
            if mol == 'N2' or mol == 'H2':
                S[mol] = np.zeros_like(S['H2O'])
            else:
                S[mol] = (interpn((self.param['opacp'][0], self.param['opact'][0], self.param['opacw'][0]), self.param['opac' + mol.lower()], np.array([PL * I2, TL * I2, I1 * self.param['opacw'][0]]).T)).T

        op = np.zeros_like(S['H2O'] * (N['H2O'] * I2))

        # CIA opacity
        (self.param['opaclc']).T[0][0] = 1e-7
        SCIA = (interpn(((self.param['opaclc']).T[0], self.param['opactc'][0]), self.param['opaccia'], np.array([I1 * self.param['opacw'][0], TL * I2]).T) * 1e-10).T
        if self.param['gas_fill'] == 'H2' and self.param['CIA_contribution']:
            op += SCIA * (N['H2'] * I2) * (N['H2'] * I2) * 0.8

        # Rayleigh opacity
        SR = I1 * 8.49e-45 * (self.param['opacw'][0] * 100.0) ** (-4.0) * 1e-4
        if self.param['Rayleigh_contribution']:
            op += SR * (NL * I2)

        # Opacity of Atmosphere(m^-1)
        for mol in self.param['fit_molecules'] + [self.param['gas_fill']]:
            if self.param[mol + '_contribution']:
                op += S[mol] * (N[mol] * I2)
            else:
                pass

        # Cloud
        if self.param['incl_clouds'] and self.param['cld_contribution']:
            if self.param['fit_gen_cld']:
                indx = find_nearest(PL, self.param['P_top'])
                op[:indx + 1, :] = 10. ** 5.
            else:
                # Particle Size
                RA, VP, DA = np.zeros(len(PL)), np.zeros(len(PL)), np.zeros(len(PL))
                for i in range(0, len(PL)):
                    _, _, RA[i], VP[i] = particlesizef(g[i], TL[i], PL[i], mul[i], 18.0, 10.0, 0.01 * N['H2O'][i] / NL[i] * PL[i])
                    DA[i] = max(min(2.0 * RA[i], 100), 0.1)

                # Opacity
                if self.param['fit_wtr_cld']:
                    tck = interp2d(self.param['opacwa'][0], self.param['opacda'][0], self.param['opacaerh2o'])
                elif self.param['fit_amn_cld']:
                    tck = interp2d(self.param['opacwa'][0], self.param['opacda'][0], self.param['opacaernh3'])
                SAER = tck(self.param['opacw'][0], DA) * 1e-4

                # Mass
                VI = np.array([VP * 1e-3]).T
                NA = (fH2O[0:n - 1] - fH2O[1:n]) * 0.018 * PL / 8.3144621 / TL / VI

                op += (SAER * (NA * I2))
        else:
            pass

        # Tholin Haze
        if self.param['incl_haze'] and self.param['haze_contribution']:
            tck_haze = interp2d(self.param['opacwa'][0], self.param['opacda'][0], self.param['opacaertholin'])
            d_haze = np.full(np.shape(NL),self.param['diam_haze'])
            f_haze = self.param['vmr_haze']
            S_haze = tck_haze(self.param['opacw'][0], np.ndarray.flatten(d_haze)) * 1e-4
            #parameterize by fraction f_haze and mass (mp*mmm/volume/density) // density assumed: 800 kg/m^3
            N_haze = f_haze * NL * const.u.value * mul / (4/3 * math.pi * (d_haze/2/1e6)**3*800)
            op += (S_haze * (N_haze * I2))

        w = self.param['opacw']
        # Total Opacity of Transit Path
        ta = np.zeros((n, len(w[0])))
        for i in range(0, n - 1):
            ta[i, :] = 2.0 * op[i, :] * np.sqrt(2.0 * (Z[i + 1] - Z[i]) * (R0 + Z[i + 1]))
            for ip in range(i + 1, n - 1):
                dx = (Z[ip + 1] - Z[ip]) * (R0 + Z[ip]) / np.sqrt((Z[ip] - Z[i]) * (2 * R0 + Z[i] + Z[ip]))
                ta[i, :] = ta[i, :] + (2.0 * op[i, :] * dx)

        # Apparent size of Planet
        sp = (math.pi * (R0 ** 2.0)) + np.zeros_like(w)
        for i in range(1, n - 1):
            # sp += 2 * math.pi * (Z[i + 1] - Z[i]) * (R0 + Z[i]) * (1 - (np.exp(-ta[i]) + np.exp(-ta[i + 1])) / 2)
            sp = sp + (2.0 * math.pi * (Z[i + 1] - Z[i]) * (R0 + Z[i]) * (1.0 - (np.exp(-ta[i, :]) + np.exp(-ta[i + 1, :])) / 2.0))

        # Transit Depth
        de = sp / math.pi / (RS ** 2.0)

        return w[0] * 1000000., de[0]

    def stellar_activity(self):
        st_het = take_star_spectrum(self.param, self.param['Ts_het'], meta=self.param['meta'])
        st_phot = take_star_spectrum(self.param, self.param['Ts_phot'], meta=self.param['meta'])
        stellar_contr = 1.0 / (1.0 - (self.param['st_frac'] * (1.0 - (st_het / st_phot))))

        return stellar_contr


def forward(param, evaluation=None, retrieval_mode=True):
    """
        Runs the forward model to generate the model spectrum.

        Parameters:
        param (dict): Dictionary of parameters containing information about the model's state and parameters.
        evaluation (dict, optional): A dictionary that contains the values of the parameters to be evaluated. Default is None.
        retrieval_mode (bool, optional): If True, only the model spectrum is returned. If False, the wavelength array and the model spectrum
                                         are returned. Default is True.

        Returns:
        np.array: The model spectrum if retrieval_mode is True.
        (np.array, np.array): Wavelength array and the model spectrum if retrieval_mode is False.

        The function performs several steps:
        1. If evaluation is not None, the function updates the parameter values based on the evaluation dictionary.
        2. It calls `cloud_pos` and `calc_mean_mol_mass` functions to update the cloud and mean molecular weight parameters.
        3. It calls the `FORWARD_MODEL` function to generate the model spectrum.
        4. If the 'wl' field of the 'spectrum' parameter is not None, the function performs a custom spectral binning using `custom_spectral_binning`.
        5. Depending on the value of retrieval_mode, the function returns either just the model spectrum or both the wavelength array and the model spectrum.
    """
    param = copy.deepcopy(param)

    if evaluation is not None:
        if param['fit_Rp']:
            param['Rp'] = evaluation['Rp']
        if param['fit_T']:
            param['Tp'] = evaluation['Tp']
        if param['fit_wtr_cld']:
            param['Pw_top'] = (10. ** evaluation['pH2O'])
            param['cldw_depth'] = (10. ** evaluation['dH2O'])
            param['CR_H2O'] = (10. ** evaluation['crH2O'])
        if param['fit_amn_cld']:
            param['Pa_top'] = (10. ** evaluation['pNH3'])
            param['clda_depth'] = (10. ** evaluation['dNH3'])
            param['CR_NH3'] = (10. ** evaluation['crNH3'])
        if param['fit_gen_cld']:
            param['P_top'] = (10. ** evaluation['ptop'])
        if param['incl_haze']:
            param['diam_haze'] = (10. ** evaluation['dhaze'])
            param['vmr_haze'] = (10. ** evaluation['vmrhaze'])
        if not param['bare_rock']:
            clr = {}
            for mol in param['fit_molecules']:
                clr[mol] = evaluation[mol]
            param = clr_to_vmr(param, clr)
        if param['fit_het_frac']:
            param['st_frac'] = evaluation['st_frac']
        if param['fit_Ts_het']:
            param['Ts_het'] = evaluation['Ts_het']
        if param['fit_Ts_phot']:
            param['Ts_phot'] = evaluation['Ts_phot']

    if not param['bare_rock']:
        param = cloud_pos(param)
        param = calc_mean_mol_mass(param)
        mod = FORWARD_MODEL(param)
        wl, trans = mod.atmospheric_structure()

        if param['spectrum']['bins']:
            model = custom_spectral_binning(np.array([param['spectrum']['wl_low'], param['spectrum']['wl_high'], param['spectrum']['wl']]).T, wl, trans, bins=param['spectrum']['bins'])
        else:
            model = custom_spectral_binning(param['spectrum']['wl'], wl, trans, bins=param['spectrum']['bins'])
    else:
        model = np.ones(len(param['spectrum']['wl'])) * (((param['Rp'] * const.R_earth.value) / (param['Rs'] * const.R_sun.value)) ** 2.)

    if param['incl_star_activity'] and param['star_act_contribution']:
        if param['bare_rock']:
            mod = FORWARD_MODEL(param)
        star_act = mod.stellar_activity()
        model = model * star_act

    if retrieval_mode:
        return model
    else:
        return param['spectrum']['wl'], model

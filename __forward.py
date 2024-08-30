from .__utils import *


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
        self.param['Ray_scat_mol'] = ['H2O', 'CO2', 'N2', 'H2']
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

        def Ray_cross(gas):
            DenS = 101325.0 / const.k_B.value / 273.0
            # wl must be in micron!!!
            if gas == 'H2':
                ray_crs = ((8.14e-13 * np.power(self.param['opacw'][0] * 1e10, -4)) + (1.28e-6 * np.power(self.param['opacw'][0] * 1e10, -6)) + (1.61 * np.power(self.param['opacw'][0] * 1e10, -8))) * 1e-4  # Dalgarno 1962, in m^2
            else:
                if gas == 'N2':
                    refidx = 1 + 6.8552e-5 + 3.243157e-2 / (144.0 - ((self.param['opacw'][0] * 1e6) ** (-2.0)))
                elif gas == 'H2O':
                    refidx = np.full(len(self.param['opacw'][0] * 1e6), 1.000261)
                elif gas == 'CO2':
                    w = 1.0 / ((self.param['opacw'][0] * 1e6) ** 2.)
                    refidx = 1 + 1e-5 * (0.154489 / (0.0584738 - w) + 8309.1927 / (210.92417 - w) + 287.64190 / (60.122959 - w))

                refidx[np.where(refidx < 1)[0]] = 1.0
                ray_crs = 1.061 * 8.0 * np.pi ** 3 * ((refidx ** 2) - 1) ** 2 / 3.0 / (self.param['opacw'][0] ** 4) / DenS / DenS

            return ray_crs

        R0 = self.param['Rp'] * const.R_earth.value  # Radius of Planet in meters
        M0 = self.param['Mp'] * const.M_earth.value  # Mass of Planet in kilograms
        RS = self.param['Rs'] * const.R_sun.value  # Radius of Star in meters

        # Atmospheric Profile
        P = np.array([self.param['P'][::-1]]).T
        if self.param['TP_profile'] is None:
            T = np.ones_like(P) * self.param['Tp']
        else:
            if self.param['fit_T']:
                T = self.param['TP_profile'](P) + self.param['Tp']
            else:
                T = self.param['TP_profile'](P)
        for i in range(0, len(P)):
            T[i] = min(max(100.0, T[i]), 2000.0)

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
        NL = PL / const.k_B.value / TL

        I1 = np.array([np.ones(len(PL))]).T
        I2 = np.ones(len(self.param['opacw'][0]))

        N = {'H2O': NL * np.sqrt(fH2O[0:n - 1] * fH2O[1:n])}
        for mol in self.param['fit_molecules']:
            if mol != 'H2O':
                N[mol] = NL * self.param['vmr_' + mol][:-1].reshape(np.shape(I1))
        ffill = np.array([self.param['vmr_' + self.param['gas_fill']][::-1]]).T
        N[self.param['gas_fill']] = NL * np.sqrt(ffill[0:n - 1] * ffill[1:n])

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
        if self.param['Rayleigh_contribution']:
            SR = np.zeros(np.shape(I1 * self.param['opacw'][0]))
            tot_vmr_cr = np.zeros(np.shape(SR))
            for mol in self.param['fit_molecules'] + [self.param['gas_fill']]:
                if mol in self.param['Ray_scat_mol']:
                    SR += (self.param['vmr_' + mol][:-1].reshape(np.shape(I1)) * Ray_cross(mol))
                    for j in range(len(self.param['opacw'][0])):
                        tot_vmr_cr[:, j] += self.param['vmr_' + mol][:-1]
                else:
                    pass

            SR /= tot_vmr_cr

            # SR = I1 * 8.49e-45 * (self.param['opacw'][0] * 100.0) ** (-4.0) * 1e-4
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
            d_haze = np.full(np.shape(NL), self.param['diam_haze'])
            f_haze = self.param['vmr_haze']
            S_haze = tck_haze(self.param['opacw'][0], np.ndarray.flatten(d_haze)) * 1e-4
            # parameterize by fraction f_haze and mass (mp*mmm/volume/density) // density assumed: 800 kg/m^3
            N_haze = f_haze * NL * const.u.value * mul / (4 / 3 * math.pi * (d_haze / 2 / 1e6) ** 3 * 800)
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
        def take_star_spectrum(param, Ts, meta=None):
            # PHOENIX STELLAR SPECTRA FROM : http://svo2.cab.inta-csic.es/theory/newov2/index.php

            def resize(fl, final=False):
                wl_temp = reso_range(0.85 * param['min_wl'], 1.1 * param['max_wl'], res=20000, bins=False)
                if not final:
                    idx_start = find_nearest(fl[:, 0] * (10. ** -4.), 0.8 * param['min_wl'])
                    idx_stop = find_nearest(fl[:, 0] * (10. ** -4.), 1.15 * param['max_wl'])
                    return spectres(wl_temp, fl[idx_start:idx_stop, 0] * (10. ** -4.), fl[idx_start:idx_stop, 1], fill=False)
                else:
                    if param['light_star_mods']:
                        wl_temp = reso_range(0.5, 5.5, res=20000, bins=False)
                    return spectres(param['spectrum']['wl'][param['sorted_data_idx']], wl_temp, fl, fill=False)

            if not param['light_star_mods']:
                directory = param['pkg_dir'] + 'PHOENIX_models/'
                skp_hdr = 6
            else:
                directory = param['pkg_dir'] + 'PHOENIX_models_light/'
                skp_hdr = 0

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
                if float(t_star1) < 26 or float(t_star2) < 26:
                    s_meta1 = '-' + str(0.0)
                    s_meta2 = '-' + str(0.0)
                s_files['A'] = np.genfromtxt(directory + 'lte' + t_star1 + loggs1 + s_meta1 + '.txt', skip_header=skp_hdr)
                s_files['B'] = np.genfromtxt(directory + 'lte' + t_star2 + loggs1 + s_meta1 + '.txt', skip_header=skp_hdr)
                s_files['C'] = np.genfromtxt(directory + 'lte' + t_star1 + loggs2 + s_meta1 + '.txt', skip_header=skp_hdr)
                s_files['D'] = np.genfromtxt(directory + 'lte' + t_star2 + loggs2 + s_meta1 + '.txt', skip_header=skp_hdr)
                s_files['E'] = np.genfromtxt(directory + 'lte' + t_star1 + loggs1 + s_meta2 + '.txt', skip_header=skp_hdr)
                s_files['F'] = np.genfromtxt(directory + 'lte' + t_star2 + loggs1 + s_meta2 + '.txt', skip_header=skp_hdr)
                s_files['G'] = np.genfromtxt(directory + 'lte' + t_star1 + loggs2 + s_meta2 + '.txt', skip_header=skp_hdr)
                s_files['H'] = np.genfromtxt(directory + 'lte' + t_star2 + loggs2 + s_meta2 + '.txt', skip_header=skp_hdr)
            elif interp_Ts and interp_loggs and not interp_meta:
                s_files['A'] = np.genfromtxt(directory + 'lte' + t_star1 + loggs1 + s_meta + '.txt', skip_header=skp_hdr)
                s_files['B'] = np.genfromtxt(directory + 'lte' + t_star2 + loggs1 + s_meta + '.txt', skip_header=skp_hdr)
                s_files['C'] = np.genfromtxt(directory + 'lte' + t_star1 + loggs2 + s_meta + '.txt', skip_header=skp_hdr)
                s_files['D'] = np.genfromtxt(directory + 'lte' + t_star2 + loggs2 + s_meta + '.txt', skip_header=skp_hdr)
            elif interp_Ts and interp_meta and not interp_loggs:
                if float(t_star1) < 26 or float(t_star2) < 26:
                    s_meta1 = '-' + str(0.0)
                    s_meta2 = '-' + str(0.0)
                s_files['A'] = np.genfromtxt(directory + 'lte' + t_star1 + loggs + s_meta1 + '.txt', skip_header=skp_hdr)
                s_files['B'] = np.genfromtxt(directory + 'lte' + t_star2 + loggs + s_meta1 + '.txt', skip_header=skp_hdr)
                s_files['C'] = np.genfromtxt(directory + 'lte' + t_star1 + loggs + s_meta2 + '.txt', skip_header=skp_hdr)
                s_files['D'] = np.genfromtxt(directory + 'lte' + t_star2 + loggs + s_meta2 + '.txt', skip_header=skp_hdr)
            elif interp_loggs and interp_meta and not interp_Ts:
                s_files['A'] = np.genfromtxt(directory + 'lte' + t_star + loggs1 + s_meta1 + '.txt', skip_header=skp_hdr)
                s_files['B'] = np.genfromtxt(directory + 'lte' + t_star + loggs2 + s_meta1 + '.txt', skip_header=skp_hdr)
                s_files['C'] = np.genfromtxt(directory + 'lte' + t_star + loggs1 + s_meta2 + '.txt', skip_header=skp_hdr)
                s_files['D'] = np.genfromtxt(directory + 'lte' + t_star + loggs2 + s_meta2 + '.txt', skip_header=skp_hdr)
            elif interp_Ts and not interp_loggs and not interp_meta:
                s_files['A'] = np.genfromtxt(directory + 'lte' + t_star1 + loggs + s_meta + '.txt', skip_header=skp_hdr)
                s_files['B'] = np.genfromtxt(directory + 'lte' + t_star2 + loggs + s_meta + '.txt', skip_header=skp_hdr)
            elif interp_loggs and not interp_Ts and not interp_meta:
                s_files['A'] = np.genfromtxt(directory + 'lte' + t_star + loggs1 + s_meta + '.txt', skip_header=skp_hdr)
                s_files['B'] = np.genfromtxt(directory + 'lte' + t_star + loggs2 + s_meta + '.txt', skip_header=skp_hdr)
            elif interp_meta and not interp_Ts and not interp_loggs:
                s_files['A'] = np.genfromtxt(directory + 'lte' + t_star + loggs + s_meta1 + '.txt', skip_header=skp_hdr)
                s_files['B'] = np.genfromtxt(directory + 'lte' + t_star + loggs + s_meta2 + '.txt', skip_header=skp_hdr)
            else:
                s_files['A'] = np.genfromtxt(directory + 'lte' + t_star + loggs + s_meta + '.txt', skip_header=skp_hdr)

            # plt.plot(s_files['D'][:, 0] * (10. ** -4.), s_files['D'][:, 1])
            # plt.plot(s_files['B'][:, 0] * (10. ** -4.), s_files['B'][:, 1])
            # plt.plot(s_files['C'][:, 0] * (10. ** -4.), s_files['C'][:, 1])
            # plt.plot(s_files['A'][:, 0] * (10. ** -4.), s_files['A'][:, 1])

            for s_keys in s_files.keys():
                if not param['light_star_mods']:
                    s_files[s_keys] = resize(s_files[s_keys], final=False)
                else:
                    s_files[s_keys] = s_files[s_keys][:,1]

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

            return resize(sp, final=True)

        if self.param['stellar_activity_parameters'] == int(3):
            st_het = take_star_spectrum(self.param, self.param['Ts_het'], meta=self.param['meta'])
            st_phot = take_star_spectrum(self.param, self.param['Ts_phot'], meta=self.param['meta'])
            st_het_effect = self.param['het_frac'] * (1.0 - (st_het / st_phot))
            return 1.0 / (1.0 - st_het_effect)

        elif self.param['stellar_activity_parameters'] == int(5):
            st_spot = take_star_spectrum(self.param, self.param['Ts_spot'], meta=self.param['meta'])
            st_fac = take_star_spectrum(self.param, self.param['Ts_fac'], meta=self.param['meta'])
            st_phot = take_star_spectrum(self.param, self.param['Ts_phot'], meta=self.param['meta'])
            st_spot_effect = self.param['spot_frac'] * (1.0 - (st_spot / st_phot))
            st_fac_effect = self.param['fac_frac'] * (1.0 - (st_fac / st_phot))
            return 1.0 / (1.0 - st_spot_effect - st_fac_effect)


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

        if param['fit_Mp']:
            param['Mp'] = evaluation['Mp']

        if param['fit_T']:
            param['Tp'] = evaluation['Tp']

        if param['fit_wtr_cld']:
            param['Pw_top'] = evaluation['pH2O']
            param['cldw_depth'] = evaluation['dH2O']
            param['CR_H2O'] = evaluation['crH2O']
        if param['fit_amn_cld']:
            param['Pa_top'] = evaluation['pNH3']
            param['clda_depth'] = evaluation['dNH3']
            param['CR_NH3'] = evaluation['crNH3']
        if param['fit_gen_cld']:
            param['P_top'] = evaluation['ptop']

        if param['incl_haze']:
            param['diam_haze'] = evaluation['dhaze']
            param['vmr_haze'] = evaluation['vmrhaze']

        if not param['bare_rock']:
            clr = {}
            for mol in param['fit_molecules']:
                clr[mol] = evaluation[mol]
            param = clr_to_vmr(param, clr)

        if param['incl_star_activity']:
            if param['stellar_activity_parameters'] == int(3):
                param['het_frac'] = evaluation['het_frac']
                param['Ts_het'] = evaluation['Ts_het']
                param['Ts_phot'] = evaluation['Ts_phot']
            if param['stellar_activity_parameters'] == int(5):
                param['spot_frac'] = evaluation['spot_frac']
                param['fac_frac'] = evaluation['fac_frac']
                param['Ts_spot'] = evaluation['Ts_spot']
                param['Ts_fac'] = evaluation['Ts_fac']
                param['Ts_phot'] = evaluation['Ts_phot']

    if not param['bare_rock']:
        param = cloud_pos(param)
        param = calc_mean_mol_mass(param)
        mod = FORWARD_MODEL(param)
        wl, trans = mod.atmospheric_structure()

        model = spectres(param['spectrum']['wl'][param['sorted_data_idx']], wl, trans, fill=False)
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
        return param['spectrum']['wl'][param['sorted_data_idx']], model

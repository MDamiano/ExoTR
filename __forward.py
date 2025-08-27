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
        self.param['Ray_scat_mol'] = ['H2O', 'CH4', 'NH3', 'SO2', 'CO', 'CO2', 'N2O', 'N2', 'He', 'H2']
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
                elif gas == 'CH4':
                    refidx = np.full(len(self.param['opacw'][0] * 1e6), 1.000444)
                elif gas == 'NH3':
                    refidx = np.full(len(self.param['opacw'][0] * 1e6), 1.000376)
                elif gas == 'CO':
                    refidx = np.full(len(self.param['opacw'][0] * 1e6), 1.000338)
                elif gas == 'CO2':
                    l = self.param['opacw'][0] * 1e6
                    l[np.where(l > 1.8)[0]] = 1.8
                    l[np.where(l < 0.48)[0]] = 0.48
                    w = 1.0 / (l ** 2.)
                    refidx = 1 + 1e-5 * (0.154489 / (0.0584738 - w) + 8309.1927 / (210.92417 - w) + 287.64190 / (60.122959 - w))
                elif gas == 'N2O':
                    refidx = np.full(len(self.param['opacw'][0] * 1e6), 1.000516)
                elif gas == 'SO2':
                    refidx = np.full(len(self.param['opacw'][0] * 1e6), 1.000646)
                elif gas == 'He':
                    refidx = 1 + 0.01470091 / (423.98 - (self.param['opacw'][0] * 1e6) ** (-2.0))

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

        n = len(P[:, 0])

        # Atmospheric Composition
        vmr = {}
        for mol in self.param['fit_molecules'] + [self.param['gas_fill']]:
            vmr[mol] = (self.param['vmr_' + mol][::-1]).reshape(n, 1)
        mu = np.array([self.param['mean_mol_weight'][::-1]]).T

        Z, g = np.zeros_like(P), np.zeros_like(P)
        for i in range(1, n):
            g[i - 1, 0] = 6.674E-11 * M0 / ((R0 + Z[i - 1, 0]) ** 2.0)
            H = 1.38E-23 * (T[i - 1, 0] + T[i, 0]) / (mu[i - 1, 0] + mu[i, 0]) / g[i - 1, 0] / 1.66E-27
            Z[i, 0] = Z[i - 1, 0] + np.log(P[i - 1, 0] / P[i, 0]) * H

        g[-1] = g[-2] - (g[-3] - g[-2])

        g = np.array([g[:-1, 0]]).T
        mul = np.array([(mu[0:n - 1, 0] + mu[1:n, 0]) / 2.0]).T
        if self.param['fit_T']:
            PL = np.array([np.sqrt(P[0:n - 1, 0] * P[1:n, 0])]).T
            TL = np.array([(T[0:n - 1, 0] + T[1:n, 0]) / 2.0]).T
        else:
            PL = self.param['forward']['PL'] + 0.0
            TL = self.param['forward']['TL'] + 0.0
        NL = PL / const.k_B.value / TL

        I1 = np.ones(len(PL)).reshape(len(PL), 1)
        I2 = np.ones(len(self.param['opacw'][0]))

        N = {}
        for mol in self.param['fit_molecules'] + [self.param['gas_fill']]:
            N[mol] = NL * np.sqrt(vmr[mol][0:n - 1] * vmr[mol][1:n])

        if self.param['fit_T']:
            # Calculate molecular opacities
            self.param['forward'] = {'S': {}}
            for mol in self.param['fit_molecules'] + [self.param['gas_fill']]:
                self.param['forward']['S'][mol] = (interpn((self.param['opacp'][0], self.param['opact'][0], self.param['opacw'][0]), self.param['opac' + mol.lower()], np.array([PL * I2, TL * I2, I1 * self.param['opacw'][0]]).T)).T
        else:
            pass

        op = np.zeros_like(I1 * I2)

        # CIA opacity
        (self.param['opaclc']).T[0][0] = 1e-7
        SCIA = (interpn(((self.param['opaclc']).T[0], self.param['opactc'][0]), self.param['opaccia'], np.array([I1 * self.param['opacw'][0], TL * I2]).T) * 1e-10).T
        # TODO: this line needs to be changed now that we can include h2 as a regular non-fill gas
        if self.param['gas_fill'] == 'H2' and self.param['CIA_contribution']:
            op += (SCIA * (N['H2'] * I2) * (N['H2'] * I2) * 0.8)

        # Rayleigh opacity
        if self.param['Rayleigh_contribution']:
            SR = np.zeros_like(I1 * I2)
            tot_vmr_cr = np.zeros_like(I1 * I2)
            for mol in self.param['fit_molecules'] + [self.param['gas_fill']]:
                if mol in self.param['Ray_scat_mol']:
                    SR += (self.param['vmr_' + mol][:-1].reshape(np.shape(I1)) * Ray_cross(mol))
                    tot_vmr_cr += (np.sqrt(vmr[mol][0:n - 1] * vmr[mol][1:n]) * I2)
                else:
                    pass

            if np.sum(tot_vmr_cr) != 0.0:
                SR /= tot_vmr_cr

            # SR = I1 * 8.49e-45 * (self.param['opacw'][0] * 100.0) ** (-4.0) * 1e-4
            op += (SR * (NL * I2))

        # Opacity of Atmosphere(m^-1)
        for mol in self.param['fit_molecules'] + [self.param['gas_fill']]:
            if self.param[mol + '_contribution']:
                op += (self.param['forward']['S'][mol] * (N[mol] * I2))
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
                    _, _, RA[i], VP[i] = self.particlesizef(g[i], TL[i], PL[i], mul[i], 18.0, 10.0, 0.01 * N['H2O'][i] / NL[i] * PL[i])
                    DA[i] = max(min(2.0 * RA[i], 100), 0.1)

                # Opacity
                if self.param['fit_wtr_cld']:
                    tck = interp2d(self.param['opacwa'][0], self.param['opacda'][0], self.param['opacaerh2o'])
                elif self.param['fit_amn_cld']:
                    tck = interp2d(self.param['opacwa'][0], self.param['opacda'][0], self.param['opacaernh3'])
                SAER = tck(self.param['opacw'][0], DA) * 1e-4

                # Mass
                VI = np.array([VP * 1e-3]).T
                NA = (vmr['H2O'][0:n - 1] - vmr['H2O'][1:n]) * 0.018 * PL / const.R.value / TL / VI

                op += (SAER * (NA * I2))
        else:
            pass

        # Tholin Haze
        if self.param['incl_haze'] and self.param['haze_contribution']:
            if self.param['fit_tholin']:
                tck_haze = interp2d(self.param['opacwa'][0], self.param['opacda'][0], self.param['opacaertholin'])
                d_haze = np.full(np.shape(NL), self.param['diam_tholin'])
                f_haze = self.param['vmr_tholin']
            elif self.param['fit_soot']:
                tck_haze = interp2d(self.param['opacwa'][0], self.param['opacda'][0], self.param['opacaersoot'])
                d_haze = np.full(np.shape(NL), self.param['diam_soot'])
                f_haze = self.param['vmr_soot']

            S_haze = tck_haze(self.param['opacw'][0], np.ndarray.flatten(d_haze)) * 1e-4
            # parameterize by fraction f_haze and mass (mp*mmm/volume/density) // density assumed: 800 kg/m^3
            N_haze = f_haze * NL * const.u.value * mul / (4 / 3 * math.pi * (d_haze / 2 / 1e6) ** 3 * 800)

            op += (S_haze * (N_haze * I2))

        # Total Opacity of Transit Path
        ta = np.zeros((n, len(self.param['opacw'][0])))
        for i in range(0, n - 1):
            ta[i, :] = 2.0 * op[i, :] * np.sqrt(2.0 * (Z[i + 1] - Z[i]) * (R0 + Z[i + 1]))
            for ip in range(i + 1, n - 1):
                dx = (Z[ip + 1] - Z[ip]) * (R0 + Z[ip]) / np.sqrt((Z[ip] - Z[i]) * (2 * R0 + Z[i] + Z[ip]))
                ta[i, :] = ta[i, :] + (2.0 * op[ip, :] * dx)

        # Apparent size of Planet
        sp = (math.pi * (R0 ** 2.0)) + np.zeros_like(self.param['opacw'])
        for i in range(1, n - 1):
            # sp += 2 * math.pi * (Z[i + 1] - Z[i]) * (R0 + Z[i]) * (1 - (np.exp(-ta[i]) + np.exp(-ta[i + 1])) / 2)
            sp += (2.0 * math.pi * (Z[i + 1] - Z[i]) * (R0 + Z[i]) * (1.0 - (np.exp(-ta[i, :]) + np.exp(-ta[i + 1, :])) / 2.0))

        # Transit Depth
        de = sp / math.pi / (RS ** 2.0)

        return self.param['opacw'][0] * 1000000., de[0]

    def particlesizef(self, g, T, P, M, MM, KE, deltaP):
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

    def stellar_activity(self):
        def take_star_spectrum(Ts, meta=None):
            # PHOENIX STELLAR SPECTRA FROM : http://svo2.cab.inta-csic.es/theory/newov2/index.php

            def resize(fl, final=False):
                wl_temp = reso_range(0.85 * self.param['min_wl'], 1.1 * self.param['max_wl'], res=20000, bins=False)
                if not final:
                    idx_start = find_nearest(fl[:, 0] * (10. ** -4.), 0.8 * self.param['min_wl'])
                    idx_stop = find_nearest(fl[:, 0] * (10. ** -4.), 1.15 * self.param['max_wl'])
                    return spectres(wl_temp, fl[idx_start:idx_stop, 0] * (10. ** -4.), fl[idx_start:idx_stop, 1], fill=False)
                else:
                    if self.param['light_star_mods']:
                        wl_temp = reso_range(0.5, 5.5, res=20000, bins=False)
                    return spectres(self.param['spectrum']['wl'][self.param['sorted_data_idx']], wl_temp, fl, fill=False)

            if not self.param['light_star_mods']:
                if self.param['stellar_spec_dir'] is not None:
                    directory = self.param['stellar_spec_dir'] + 'PHOENIX_models/'
                else:
                    directory = self.param['pkg_dir'] + 'PHOENIX_models/'
                skp_hdr = 6
            else:
                if self.param['stellar_spec_dir'] is not None:
                    directory = self.param['stellar_spec_dir'] + 'PHOENIX_models_light/'
                else:
                    directory = self.param['pkg_dir'] + 'PHOENIX_models_light/'
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
                self.param['Loggs'] += 0.0
            except KeyError:
                self.param['Loggs'] = np.log10(const.G.value * (self.param['Ms'] * const.M_sun.value) / ((self.param['Rs'] * const.R_sun.value) ** 2.)) + 2

            loggs = round(self.param['Loggs'] * 2.) / 2.
            interp_loggs = True
            if loggs - self.param['Loggs'] < 0:
                loggs1 = loggs + 0.0
                loggs2 = loggs + 0.5
            elif loggs - self.param['Loggs'] > 0:
                loggs1 = loggs - 0.5
                loggs2 = loggs + 0.0
            else:
                loggs = '-' + str(loggs)
                interp_loggs = False

            if interp_loggs:
                loggs2_ratio = 1.0 - ((loggs2 - self.param['Loggs']) / (loggs2 - loggs1))
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
                if not self.param['light_star_mods']:
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
            st_het = take_star_spectrum(self.param['Ts_het'], meta=self.param['meta'])
            st_phot = take_star_spectrum(self.param['Ts_phot'], meta=self.param['meta'])
            st_het_effect = self.param['het_frac'] * (1.0 - (st_het / st_phot))
            return 1.0 / (1.0 - st_het_effect)

        elif self.param['stellar_activity_parameters'] == int(5):
            st_spot = take_star_spectrum(self.param['Ts_spot'], meta=self.param['meta'])
            st_fac = take_star_spectrum(self.param['Ts_fac'], meta=self.param['meta'])
            st_phot = take_star_spectrum(self.param['Ts_phot'], meta=self.param['meta'])
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

        if param['fit_tholin']:
            param['diam_tholin'] = evaluation['dtholin']
            param['vmr_tholin'] = evaluation['vmrtholin']

        if param['fit_soot']:
            param['diam_soot'] = evaluation['dsoot']
            param['vmr_soot'] = evaluation['vmrsoot']

        if not param['bare_rock']:
            clogr = {}
            for mol in param['fit_molecules']:
                clogr[mol] = evaluation[mol]
            if param['gas_par_space'] == 'clr':
                param = clr_to_vmr(param, clogr)
            elif param['gas_par_space'] == 'vmr':
                for mol in param['fit_molecules']:
                    param['vmr_' + mol] = evaluation[mol]
                if param['gas_fill'] is not None:
                    param['vmr_' + param['gas_fill']] = evaluation[param['gas_fill']]

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
        if param['incl_clouds'] and param['cloud_type'] == 'water':
            param = cloud_pos(param)
        param = calc_mean_mol_mass(param)
        mod = FORWARD_MODEL(param)
        wl, trans = mod.atmospheric_structure()
        model = spectres(param['spectrum']['wl'][param['sorted_data_idx']], wl, trans, fill=False, verbose=False)
    else:
        model = np.ones(len(param['spectrum']['wl'])) * (((param['Rp'] * const.R_earth.value) / (param['Rs'] * const.R_sun.value)) ** 2.)

    if param['incl_star_activity'] and param['star_act_contribution']:
        if param['bare_rock']:
            mod = FORWARD_MODEL(param)
        star_act = mod.stellar_activity()
        model = model * star_act

    if platform.system() == 'Darwin' and np.isnan(np.sum(model)):
        print('NaN detected')

    if retrieval_mode:
        return model
    else:
        return param['spectrum']['wl'][param['sorted_data_idx']], model

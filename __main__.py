from .__basics import *
from .__utils import *
from .__forward import *

path = os.path.abspath(__file__)
pkg_dir = os.path.dirname(path) + '/'
param = take_standard_parameters(pkg_dir)
param['pkg_dir'] = pkg_dir


class RETRIEVAL:
    """
        This class executes the retrieval analysis using the selected optimization algorithm.

        Attributes
        ----------
        param: dict
            Dictionary containing all the necessary parameters for the retrieval.

        Methods
        -------
        __init__(param=param):
            Initializes the RETRIEVAL object with the provided parameters and sets the retrieval mode to True.

        run_retrieval(parfile):
            Executes the retrieval process using the optimizer defined in the parameters, currently only supports 'multinest'.
    """
    def __init__(self, param=param):
        self.param = param
        self.param['ret_mode'] = True

    def run_retrieval(self, parfile):
        """
            Executes the retrieval process using the optimizer defined in the parameters.

            Parameters
            ----------
            parfile: str
                The name of the file containing the necessary parameters for the retrieval.

            Returns
            -------
            None
        """
        self.param = read_parfile(self.param, parfile)
        if self.param['optimizer'] == 'multinest':
            from ExoTR.__multinest import MULTINEST
            bayes = MULTINEST(self.param)
            bayes.run_retrieval()
        else:
            pass


class CREATE_SPECTRUM:
    """
        This class is responsible for creating a synthetic spectrum of an exoplanet's atmosphere based on a set of
        parameters.

        Parameters
        ----------
        param: dict, optional
            A dictionary containing parameters for the simulation. If not provided, the 'param' object at the module
            level is used. Defaults to 'param'.

        plot_all: bool, optional
            A flag to decide whether to plot all intermediate results or not. Defaults to False.

        verbose: bool, optional
            A flag to decide whether to print verbose output to the console. Defaults to True.

        Methods
        -------
        __init__(param=param):
            Initializes the RETRIEVAL object with the provided parameters and sets the retrieval mode to True.
        run_forward(parfile):
            Executes the forward model calculation using the parameters specified in the file.
    """
    def __init__(self, param=param, plot_all=False, verbose=True):
        self.param = param
        self.param['ret_mode'] = False
        self.verbose = verbose
        self.plot_all = plot_all

    def run_forward(self, parfile):
        """
            Executes the forward model calculation using the parameters specified in the file.

            Parameters
            ----------
            parfile: str
                The name of the file containing the necessary parameters for the forward calculation.

            Returns
            -------
            None
        """
        self.param = read_parfile(self.param, parfile)
        self.param = par_and_calc(self.param)
        self.param = load_input_spectrum(self.param)
        self.param = pre_load_variables(param)

        if param['incl_star_activity'] and param['st_frac'] is None:
            raise KeyError('If you want to include the star activity in the calculation, please specify the fraction of the star covered by the activity through the parameter "st_frac"')
        if param['incl_star_activity'] and param['Ts_het'] is None:
            raise KeyError('If you want to include the star activity in the calculation, please specify the temperature of star hets through the parameter "Ts_het"')
        if param['incl_star_activity'] and param['Ts_phot'] is None:
            raise KeyError('If you want to include the star activity in the calculation, please specify the temperature of star through the parameter "Ts"')

        if self.verbose:
            print('Planet : ' + self.param['name_p'])
            print('Calculating the planetary transmission spectrum')
            print('Parameters:')

        for mol in self.param['supported_molecules']:
            try:
                self.param['vmr_' + mol] += 0.0
            except KeyError:
                self.param['vmr_' + mol] = 0.0

        if self.param['fit_wtr_cld']:
            self.param['Pw_top'] = 10. ** self.param['pH2O']
            self.param['cldw_depth'] = 10. ** self.param['dH2O']
            self.param['CR_H2O'] = 10. ** self.param['crH2O']
        if self.param['fit_amn_cld']:
            self.param['Pa_top'] = 10. ** self.param['pNH3']
            self.param['clda_depth'] = 10. ** self.param['dNH3']
            self.param['CR_NH3'] = 10. ** self.param['crNH3']
        if self.param['fit_g']:
            self.param['gp'] = (10. ** (self.param['g'] - 2.0))

        if self.param['fit_wtr_cld']:
            print('Log(H2O_Ptop) \t = \t' + str(self.param['pH2O']))
            print('Log(H2O_D) \t = \t' + str(self.param['dH2O']))
            print('Log(H2O_CR) \t = \t' + str(self.param['crH2O']))
        elif self.param['fit_amn_cld']:
            print('Log(NH3_Ptop) \t = \t' + str(self.param['pH2O']))
            print('Log(NH3_D) \t = \t' + str(self.param['dH2O']))
            print('Log(HN3_CR) \t = \t' + str(self.param['crH2O']))
        elif self.param['fit_gen_cld']:
            print('Log(Ptop) \t = \t' + str(self.param['P_top']))

        if self.param['incl_haze']:
            print('diam_haze \t = \t' + str(self.param['diam_haze']))
            print('vmr_haze \t = \t' + str(self.param['vmr_haze']))

        try:
            print('g \t\t = \t' + str(self.param['gp']))
        except KeyError:
            self.param['gp'] = (const.G.value * const.M_earth.value * param['Mp']) / ((const.R_earth.value * param['Rp']) ** 2.)  # g is in m/s2
            print('g \t\t = \t' + str(self.param['gp']))
        print('T \t\t = \t' + str(self.param['Tp']))

        self.param['fit_molecules'] = []
        for mol in self.param['supported_molecules']:
            try:
                if self.param['vmr_' + mol] != 0.0:
                    self.param['fit_molecules'].append(mol)
            except KeyError:
                pass

        for mol in self.param['fit_molecules']:
            if mol == 'N2' and self.param['gas_fill'] != 'N2' and self.param['vmr_N2'] != 0:
                print('VMR N2 \t\t = \t' + str(self.param['vmr_N2']))
            elif mol == 'N2' and self.param['gas_fill'] == 'N2':
                pass
            elif mol == 'O2' or mol == 'O3' or mol == 'CO':
                if self.param['vmr_' + mol] != 0.0:
                    print('VMR ' + mol + ' \t\t = \t' + str(self.param['vmr_' + mol]))
                else:
                    pass
            else:
                if self.param['vmr_' + mol] != 0.0:
                    print('VMR ' + mol + ' \t = \t' + str(self.param['vmr_' + mol]))
                else:
                    pass

        time1 = time.time()

        wl, model = forward(self.param, retrieval_mode=self.param['ret_mode'])

        time2 = time.time()

        if self.verbose:
            elapsed((time2 - time1) * (10 ** 9))

        data = np.array([wl, model]).T

        if self.param['add_noise']:
            data = add_noise(self.param, data)

        if self.param['spectrum']['bins'] and self.param['return_bins']:
            data = np.concatenate((np.array([self.param['spectrum']['wl_high']]).T, data), axis=1)
            data = np.concatenate((np.array([self.param['spectrum']['wl_low']]).T, data), axis=1)

        try:
            if not os.path.exists(self.param['out_dir']):
                os.mkdir(self.param['out_dir'])
            np.savetxt(self.param['out_dir'] + str(self.param['file_output_name']) + '.dat', data)
            print('The spectrum file has been saved in ' + self.param['out_dir'] + str(self.param['file_output_name']) + '.dat')
        except IOError or KeyError:
            if not os.path.exists(self.param['pkg_dir'] + 'Output/'):
                os.mkdir(self.param['pkg_dir'] + 'Output/')

            np.savetxt(self.param['pkg_dir'] + 'Output/spectrum.dat', data)
            print('The spectrum file has been saved in ' + self.param['pkg_dir'] + 'Output/spectrum.dat')

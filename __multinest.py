import numpy as np

from .__basics import *
from .__utils import *
from .__forward import *

# trying to initiate MPI parallelization
try:
    from mpi4py import MPI

    MPIrank = MPI.COMM_WORLD.Get_rank()
    MPIsize = MPI.COMM_WORLD.Get_size()
    MPIimport = True
except ImportError:
    MPIimport = False

if MPIimport:
    if MPIrank == 1:
        print('MPI enabled. Running on ' + str(MPIsize) + ' cores')
else:
    print('MPI disabled')

# checking for multinest library
try:
    import pymultinest

    multinest_import = True
except:
    multinest_import = False

if multinest_import:
    if MPIrank == 1:
        from pymultinest.run import lib_mpi

        print('MultiNest library: "' + str(lib_mpi) + '" correctly loaded.')
else:
    print('SOME ERRORS OCCURRED - MultiNest library is not loaded.')
    raise ImportError


class MULTINEST:
    """
        This class implements methods to perform a Bayesian retrieval analysis using the Multinest algorithm.

        Attributes
        ----------
        param: dict
            Dictionary containing all the necessary parameters for the retrieval.

        Methods
        -------
        __init__(param):
            Initializes the MULTINEST object with the provided parameters and sets up the necessary variables for the retrieval.

        run_retrieval():
            Executes the retrieval process using the nested sampling algorithm implemented in PyMultiNest.
    """
    def __init__(self, param):
        self.param = param
        self.param = par_and_calc(self.param)
        self.param = load_input_spectrum(self.param)
        self.param = pre_load_variables(param)
        self.param = ranges(self.param)

    def run_retrieval(self):
        """
            Executes the retrieval process using a nested sampling algorithm implemented in PyMultiNest.

            This method involves the following steps:
            1) Synchronizes all MPI processes.
            2) Broadcasts parameters to all MPI processes.
            3) Loads the prior for cloud to volume mixing ratio (clr) conversion.
            4) Runs the nested sampling algorithm using the loglike and prior methods for likelihood and prior evaluation.
            5) If rank is 0 (main process), it saves the parameter names, plots the spectra and posterior distribution functions.

            Notes
            -----
            This method uses MPI for parallelization. If MPI is not imported, the retrieval is run on a single processor.
            The retrieval is run using PyMultiNest, which performs Bayesian inference using the nested sampling algorithm.
        """
        if MPIimport and MPIrank == 0:
            print('Using ExoTR for transmission spectroscopy retrieval')

        if MPIimport:
            MPI.COMM_WORLD.Barrier()  # wait for everybody to synchronize here

        self.param = MPI.COMM_WORLD.bcast(self.param, root=0)

        if MPIimport:
            MPI.COMM_WORLD.Barrier()  # wait for everybody to synchronize here

        parameters, n_params = retrieval_par_and_npar(self.param)

        if len(self.param['fit_molecules']) > 0 and self.param['modified_clr_prior']:
            ppf = define_modified_clr_prior(len(self.param['fit_molecules']))

        def internal_model(cube, retrieval_mode=True):
            """
                Generates the forward model given a set of parameters.

                Parameters
                ----------
                cube : array_like
                    The input parameter array. The length and order of values depend on the parameters to be fitted.
                retrieval_mode : bool, optional
                    If True, the function will execute in retrieval mode. The default is True.

                Returns
                -------
                model : array_like
                    The forward model spectrum.

                Notes
                -----
                This function is primarily used in the context of retrieval, where it is used to generate the model spectrum
                for a given set of parameters in the MCMC or nested sampling process.
            """
            evaluation = {}
            par = 0
            if self.param['fit_offset']:
                for i in range(0, self.param['n_offsets']):
                    self.param['offset' + str(i + 1)] = cube[par] + 0.0  # offset of datasets
                    par += 1

            if self.param['fit_Rp']:
                evaluation['Rp'] = cube[par] + 0.0  # Rp, Planetary radius
                par += 1

            if self.param['fit_Mp']:
                evaluation['Mp'] = cube[par] + 0.0  # Mp, Planetary mass
                par += 1

            if self.param['fit_T']:
                evaluation['Tp'] = cube[par] + 0.0  # Planetary temperature or delta temperature
                par += 1

            if self.param['fit_wtr_cld']:
                evaluation['pH2O'], evaluation['dH2O'], evaluation['crH2O'] = (10. ** cube[par]), (10. ** cube[par + 1]), (10. ** cube[par + 2])  # pH2O, dH2O, crH2O
                par += 3
            if self.param['fit_amn_cld']:
                evaluation['pNH3'], evaluation['dNH3'], evaluation['crNH3'] = (10. ** cube[par]), (10. ** cube[par + 1]), (10. ** cube[par + 2])  # pNH3, dNH3, crNH3
                par += 3
            if self.param['fit_gen_cld']:
                evaluation['ptop'] = (10. ** cube[par])  # p_top
                par += 1

            if self.param['incl_haze']:
                evaluation['dhaze'], evaluation['vmrhaze'] = (10. ** cube[par]), (10. ** cube[par + 1])  # diameter and haze VMR
                par += 2

            for mol in self.param['fit_molecules']:
                evaluation[mol] = cube[par] + 0.0  # Molecules
                par += 1

            if self.param['incl_star_activity']:
                if self.param['stellar_activity_parameters'] == int(3):
                    evaluation['het_frac'] = cube[par] + 0.0  # heterogeneity coverage fraction
                    evaluation['Ts_het'] = cube[par + 1] + 0.0  # heterogeneity temperature
                    evaluation['Ts_phot'] = cube[par + 2] + 0.0  # star photosphere temperature
                    par += 3

                elif self.param['stellar_activity_parameters'] == int(5):
                    evaluation['spot_frac'] = cube[par] + 0.0  # spot coverage fraction
                    evaluation['fac_frac'] = cube[par + 1] + 0.0  # faculae coverage fraction
                    evaluation['Ts_spot'] = cube[par + 2] + 0.0  # spot temperature
                    evaluation['Ts_fac'] = cube[par + 3] + 0.0  # faculae temperature
                    evaluation['Ts_phot'] = cube[par + 4] + 0.0  # star photosphere temperature
                    par += 5

            if platform.system() == 'Darwin':
                return forward(self.param, evaluation=evaluation, retrieval_mode=retrieval_mode)
            else:
                try:
                    return forward(self.param, evaluation=evaluation, retrieval_mode=retrieval_mode)
                except:
                    if MPIimport:
                        MPI.Finalize()
                        sys.exit()
                    else:
                        print('Some errors occurred in during the calculation of the forward model.')
                        sys.exit()

        def prior(cube, ndim, nparams):
            """
                Applies prior transformations to the parameter array for a given number of gases and parameters.

                Parameters
                ----------
                cube : array_like
                    The input parameter array. The length and order of values depend on the parameters to be fitted.
                ndim : int
                    The number of dimensions (parameters).
                nparams : int
                    The number of parameters.
                n_gas : int
                    The number of gases in the model, default is from parameter dictionary.

                Returns
                -------
                None.

                Notes
                -----
                This function directly modifies the input `cube` based on the priors defined in the parameter dictionary.
                The priors are assumed to be uniform for most parameters, except the clr of the gases, which uses a modified prior.
            """
            par = 0
            if self.param['fit_offset']:
                for i in range(0, self.param['n_offsets']):
                    cube[par] = cube[par] * (self.param['off' + str(i + 1) + '_range'][1] - self.param['off' + str(i + 1) + '_range'][0]) + self.param['off' + str(i + 1) + '_range'][0]  # uniform prior -> offsets between differen datasets
                    par += 1

            if self.param['fit_Rp']:
                cube[par] = cube[par] * (self.param['rp_range'][1] - self.param['rp_range'][0]) + self.param['rp_range'][0]  # uniform prior -> Rp, Planetary Radius
                par += 1

            if self.param['fit_Mp']:
                if self.param['Mp_prior'] == 'uniform':
                    cube[par] = cube[par] * (self.param['mp_range'][1] - self.param['mp_range'][0]) + self.param['mp_range'][0]  # uniform prior -> Mp, Planetary Radius
                elif self.param['Mp_prior'] == 'gaussian':
                    Mp_range = np.linspace(self.param['mp_range'][0], self.param['mp_range'][1], num=10000, endpoint=True)                      # gaussian prior -> Mp, Planetary Radius
                    Mp_cdf = sp.stats.norm.cdf(Mp_range, self.param['Mp_orig'], self.param['Mp_err'])
                    Mp_cdf = np.array([0.0] + list(Mp_cdf) + [1.0])
                    Mp_range = np.array([Mp_range[0]] + list(Mp_range) + [Mp_range[-1]])
                    Mp_pri = interp1d(Mp_cdf, Mp_range)
                    cube[par] = Mp_pri(cube[par])
                par += 1

            if self.param['fit_T']:
                cube[par] = cube[par] * (self.param['tp_range'][1] - self.param['tp_range'][0]) + self.param['tp_range'][0]  # uniform prior -> Tp, Planetary temperature
                par += 1

            if self.param['fit_wtr_cld']:
                cube[par] = cube[par] * (self.param['ptopw_range'][1] - self.param['ptopw_range'][0]) + self.param['ptopw_range'][0]  # uniform prior -> P H2O cloud top [Pa]
                cube[par + 1] = cube[par + 1] * (self.param['dcldw_range'][1] - self.param['dcldw_range'][0]) + self.param['dcldw_range'][0]  # uniform prior -> D H2O cloud [Pa]
                cube[par + 2] = cube[par + 2] * (self.param['crh2o_range'][1] - self.param['crh2o_range'][0]) + self.param['crh2o_range'][0]  # uniform prior -> CR H2O
                par += 3
            if self.param['fit_amn_cld']:
                cube[par] = cube[par] * (self.param['ptopa_range'][1] - self.param['ptopa_range'][0]) + self.param['ptopa_range'][0]  # uniform prior -> P NH3 cloud top [Pa]
                cube[par + 1] = cube[par + 1] * (self.param['dclda_range'][1] - self.param['dclda_range'][0]) + self.param['dclda_range'][0]  # uniform prior -> D NH3 cloud [Pa]
                cube[par + 2] = cube[par + 2] * (self.param['crnh3_range'][1] - self.param['crnh3_range'][0]) + self.param['crnh3_range'][0]  # uniform prior -> CR NH3
                par += 3
            if self.param['fit_gen_cld']:
                cube[par] = cube[par] * (self.param['ptop_range'][1] - self.param['ptop_range'][0]) + self.param['ptop_range'][0]  # uniform prior -> P cloud top [Pa]
                par += 1

            if self.param['incl_haze']:
                cube[par] = cube[par] * (self.param['dhaze_range'][1] - self.param['dhaze_range'][0]) + self.param['dhaze_range'][0]
                cube[par + 1] = cube[par + 1] * (self.param['vmrhaze_range'][1] - self.param['vmrhaze_range'][0]) + self.param['vmrhaze_range'][0]
                par += 2

            for _ in self.param['fit_molecules']:
                if self.param['modified_clr_prior']:
                    cube[par] = ppf(cube[par])  # modified prior for clr
                else:
                    cube[par] = cube[par] * (self.param['gas_clr_range'][1] - self.param['gas_clr_range'][0]) + self.param['gas_clr_range'][0]  # standard clr range
                par += 1

            if self.param['incl_star_activity']:
                if self.param['stellar_activity_parameters'] == int(3):
                    cube[par] = cube[par] * (self.param['delta_range'][1] - self.param['delta_range'][0]) + self.param['delta_range'][0]  # uniform prior -> star activity coverage fraction
                    par += 1
                    cube[par], cube[par + 1] = Ts_prior(self.param, cube[par + 1], Ts_het_cube=cube[par])
                    par += 2

                elif self.param['stellar_activity_parameters'] == int(5):
                    cube[par] = cube[par] * (self.param['delta_range'][1] - self.param['delta_range'][0]) + self.param['delta_range'][0]  # uniform prior -> spots coverage fraction
                    par += 1
                    cube[par] = cube[par] * (self.param['delta_range'][1] - self.param['delta_range'][0]) + self.param['delta_range'][0]  # uniform prior -> faculae coverage fraction
                    par += 1
                    cube[par], cube[par + 1], cube[par + 2] = Ts_prior(self.param, cube[par + 2], Ts_spot_cube=cube[par], Ts_fac_cube=cube[par + 1])
                    par += 3

        def loglike(cube, ndim, nparams):
            """
                Computes the log-likelihood for given input parameters.

                Parameters
                ----------
                cube : array_like
                    The input parameter array. The length and order of values depend on the parameters to be fitted.
                ndim : int
                    The number of dimensions (parameters).
                nparams : int
                    The number of parameters.

                Returns
                -------
                float
                    The computed log-likelihood.

                Notes
                -----
                The likelihood is computed assuming Gaussian errors. The internal_model function is supposed to compute the model predictions based on the input parameters. The residuals (difference between the data and the model) are normalized by the errors and then used to compute the likelihood.
            """
            if self.param['fit_offset']:
                model = internal_model(cube)
                data_spec = self.param['spectrum1']['T_depth'] + 0.0
                for i in range(0, self.param['n_offsets']):
                    temp_s = self.param['spectrum' + str(i + 2)]['T_depth'] + self.param['offset' + str(i + 1)]
                    data_spec = np.append(data_spec, temp_s)

                chi = (data_spec[self.param['sorted_data_idx']] - model) / self.param['spectrum']['error_T'][self.param['sorted_data_idx']]
            else:
                chi = (self.param['spectrum']['T_depth'][self.param['sorted_data_idx']] - internal_model(cube)) / self.param['spectrum']['error_T'][self.param['sorted_data_idx']]

            loglikelihood = (-1.) * np.sum(np.log(self.param['spectrum']['error_T'][self.param['sorted_data_idx']] * np.sqrt(2.0 * math.pi))) - 0.5 * np.sum(chi * chi)

            return loglikelihood

        if MPIimport and MPIrank == 0:
            time1 = time.time()

        pymultinest.run(LogLikelihood=loglike,
                        Prior=prior,
                        n_dims=n_params,
                        multimodal=self.param['multimodal'],
                        max_modes=self.param['max_modes'],
                        outputfiles_basename=self.param['out_dir'] + self.param['name_p'] + '_',
                        importance_nested_sampling=False,
                        evidence_tolerance=self.param['ev_tolerance'],
                        n_live_points=self.param['nlive_p'],
                        resume=self.param['multinest_resume'],
                        verbose=self.param['multinest_verbose'],
                        init_MPI=False)

        if MPIimport and MPIrank == 0:  # Plot Nest_spectrum
            time2 = time.time()
            elapsed((time2 - time1) * (10 ** 9))

        prefix = self.param['out_dir'] + self.param['name_p'] + '_'
        if MPIimport and MPIrank == 0:
            json.dump(parameters, open(prefix + 'params.json', 'w'))  # save parameter names

        ### PRODUCE PLOTS FROM HERE --- POST-PROCESSING ###
        self.param['model_n_par'] = len(parameters)

        multinest_results = pymultinest.Analyzer(n_params=self.param['model_n_par'], outputfiles_basename=prefix, verbose=False)
        mc_samp = multinest_results.get_equal_weighted_posterior()[:, :-1]
        mds = len(multinest_results.get_mode_stats()['modes'])

        if self.param['calc_likelihood_data']:
            self.calc_spectra(mc_samp)

            if MPIimport:
                MPI.COMM_WORLD.Barrier()  # wait for everybody to synchronize here

            if MPIimport and MPIrank == 0:
                if platform.system() != 'Darwin':
                    time.sleep(600)
                rank_0 = np.loadtxt(self.param['out_dir'] + 'loglikelihood_per_datapoint/loglike_0.dat')
                rank_0_spec = np.loadtxt(self.param['out_dir'] + 'loglikelihood_per_datapoint/samples_0.dat')
                for i in range(1, MPIsize):
                    rank_n = np.loadtxt(self.param['out_dir'] + 'loglikelihood_per_datapoint/loglike_' + str(i) + '.dat')
                    rank_n_spec = np.loadtxt(self.param['out_dir'] + 'loglikelihood_per_datapoint/samples_' + str(i) + '.dat')
                    rank_0 = np.concatenate((rank_0, rank_n), axis=0)
                    rank_0_spec = np.concatenate((rank_0_spec, rank_n_spec[:, 1:]), axis=1)
                np.savetxt(self.param['out_dir'] + 'loglike_per_datapoint.dat', rank_0)
                np.savetxt(self.param['out_dir'] + 'random_samples.dat', rank_0_spec)
                os.system('rm -rf ' + self.param['out_dir'] + 'loglikelihood_per_datapoint/')

                self.param['spec_sample'] = rank_0_spec + 0.0
                del rank_0_spec, rank_0

        if MPIimport and MPIrank == 0:
            if self.param['plot_models']:
                self.param['chi_square_stat'] = {}
                if mds < 2:
                    # s = multinest_results.get_best_fit()
                    # cube = s['parameters']
                    # cube = np.array([cube, ]).T

                    # s = multinest_results.get_stats()
                    # cube = []
                    # for p, m in zip(parameters, s['marginals']):
                    #     cube.append(m['median'])
                    # cube = np.array([cube, ]).T

                    s = multinest_results.get_stats()
                    cube = np.ones((len(s['modes'][0]['maximum a posterior']), mds))
                    cube[:, 0] = list(s['modes'][0]['maximum a posterior'])

                    self.plot_nest_spec(cube[:, 0])
                    if not self.param['bare_rock']:
                        self.plot_P_profiles()
                        self.plot_contribution()
                else:
                    s = multinest_results.get_mode_stats()
                    cube = np.ones((len(s['modes'][0]['maximum a posterior']), mds))
                    for i in range(0, mds):
                        cube[:, i] = list(s['modes'][i]['maximum a posterior'])

                        self.plot_nest_spec(cube[:, i], solutions=i + 1)
                        if not self.param['bare_rock']:
                            self.plot_P_profiles(solutions=i + 1)
                            self.plot_contribution(solutions=i + 1)

                if self.param['spectrum']['bins']:
                    data_spec = np.array([self.param['spectrum']['wl_low'][self.param['sorted_data_idx']], self.param['spectrum']['wl_high'][self.param['sorted_data_idx']], self.param['spectrum']['wl'][self.param['sorted_data_idx']], self.param['spectrum']['T_depth'][self.param['sorted_data_idx']], self.param['spectrum']['error_T'][self.param['sorted_data_idx']]]).T
                else:
                    data_spec = np.array([self.param['spectrum']['wl'][self.param['sorted_data_idx']], self.param['spectrum']['T_depth'][self.param['sorted_data_idx']], self.param['spectrum']['error_T'][self.param['sorted_data_idx']]]).T
                np.savetxt(self.param['out_dir'] + 'data_spectrum.dat', data_spec)

            if self.param['plot_posterior']:
                self.plot_posteriors(prefix, multinest_results, parameters, mds)

        if MPIimport:
            MPI.Finalize()

    def cube_to_param(self, cube):
        """
            Converts the input cube (parameter array) to the actual parameter values used in the model.

            Parameters
            ----------
            cube : array_like
                The input parameter array. The length and order of values depend on the parameters to be fitted.

            Returns
            -------
            None. This function modifies self.param in-place.

            Notes
            -----
            self.param is a dictionary containing the model parameters.
        """

        par = 0
        if self.param['fit_offset']:
            for i in range(0, self.param['n_offsets']):
                self.param['offset' + str(i + 1)] = cube[par] + 0.0  # offset of datasets
                par += 1
            self.param['spectrum']['T_depth'] = self.param['spectrum1']['T_depth'] + 0.0
            for i in range(0, self.param['n_offsets']):
                temp_s = self.param['spectrum' + str(i + 2)]['T_depth'] + self.param['offset' + str(i + 1)]
                self.param['spectrum']['T_depth'] = np.append(self.param['spectrum']['T_depth'], temp_s)

        if self.param['fit_Rp']:
            self.param['Rp'] = cube[par] + 0.0
            par += 1

        if self.param['fit_Mp']:
            self.param['Mp'] = cube[par] + 0.0
            par += 1

        if self.param['fit_T']:
            self.param['Tp'] = cube[par] + 0.0  # Planetary temperature
            par += 1

        if self.param['fit_wtr_cld']:
            self.param['Pw_top'] = (10. ** cube[par])
            self.param['cldw_depth'] = (10. ** cube[par + 1])
            self.param['CR_H2O'] = (10. ** cube[par + 2])
            par += 3
        if self.param['fit_amn_cld']:
            self.param['Pa_top'] = (10. ** cube[par])
            self.param['clda_depth'] = (10. ** cube[par + 1])
            self.param['CR_NH3'] = (10. ** cube[par + 2])
            par += 3
        if self.param['fit_gen_cld']:
            self.param['P_top'] = (10. ** cube[par])  # p_top
            par += 1

        if self.param['incl_haze']:
            self.param['diam_haze'] = (10. ** cube[par])
            self.param['vmr_haze'] = (10. ** cube[par + 1])
            par += 2

        if not self.param['bare_rock']:
            clr = {}
            for i in self.param['fit_molecules']:
                clr[i] = cube[par] + 0.0
                par += 1
            self.param = clr_to_vmr(self.param, clr)
            self.param = cloud_pos(self.param)
            self.param = calc_mean_mol_mass(self.param)

        if self.param['incl_star_activity']:
            if self.param['stellar_activity_parameters'] == int(3):
                self.param['het_frac'] = cube[par] + 0.0
                self.param['Ts_het'] = cube[par + 1] + 0.0
                self.param['Ts_phot'] = cube[par + 2] + 0.0
                par += 3

            elif self.param['stellar_activity_parameters'] == int(5):
                self.param['spot_frac'] = cube[par] + 0.0
                self.param['fac_frac'] = cube[par + 1] + 0.0
                self.param['Ts_spot'] = cube[par + 2] + 0.0
                self.param['Ts_fac'] = cube[par + 3] + 0.0
                self.param['Ts_phot'] = cube[par + 4] + 0.0
                par += 5

    def plot_nest_spec(self, cube, solutions=None):
        """
            Plots the observed and modeled spectral data.

            Parameters
            ----------
            cube : array_like
                The input solution in the form of a cube (or parameter array).

            solutions : int, optional
                The specific solution number to be considered, if there are multiple solutions.

            Returns
            -------
            None. This function outputs a plot saved in the output directory.

            Notes
            -----
            The plot includes observed spectral data, the modeled data at the solution points, and the modeled data at higher resolution.
        """

        def _chi_square_stat(data, gen_mod, npar=self.param['model_n_par']):
            chi = (data[:, 1] - gen_mod) / data[:, 2]
            dict_stat = {'dof': len(data[:, 1]) - npar, 'chi2': np.sum(chi ** 2.)}
            dict_stat['chi2_red'] = dict_stat['chi2'] / dict_stat['dof']

            if solutions is None:
                self.param['chi_square_stat'] = dict_stat
            else:
                self.param['chi_square_stat']['solution_' + str(solutions)] = dict_stat

        self.cube_to_param(cube)
        fig = plt.figure(figsize=(8, 5))

        plt.errorbar(self.param['spectrum']['wl'], self.param['spectrum']['T_depth']*1e6, yerr=self.param['spectrum']['error_T']*1e6,
                     linestyle='', linewidth=0.5, color='black', marker='o', markerfacecolor='red', markersize=4, capsize=1.75, label='Data')

        _, model = forward(self.param, retrieval_mode=False)
        # plt.plot(self.param['spectrum']['wl'][self.param['sorted_data_idx']], model*1e6, linestyle='', color='black', marker='d', markerfacecolor='blue', markersize=4)

        _chi_square_stat(np.array([self.param['spectrum']['wl'][self.param['sorted_data_idx']], self.param['spectrum']['T_depth'][self.param['sorted_data_idx']], self.param['spectrum']['error_T'][self.param['sorted_data_idx']]]).T, model)

        new_wl = np.loadtxt(self.param['pkg_dir'] + 'Data/wl_bins/bins_02_20_R500.dat')
        new_wl_central = np.mean(new_wl, axis=1)
        start = find_nearest(new_wl_central, min(self.param['spectrum']['wl']) - 0.05)
        stop = find_nearest(new_wl_central, max(self.param['spectrum']['wl']) + 0.05)
        if self.param['spectrum']['bins']:
            temp = np.array([self.param['spectrum']['wl_low'], self.param['spectrum']['wl_high'], self.param['spectrum']['wl']]).T
        else:
            temp = self.param['spectrum']['wl'] + 0.0
        self.param['spectrum']['wl'] = new_wl_central[start:stop]
        self.param['spectrum']['wl_low'] = new_wl[start:stop, 0]
        self.param['spectrum']['wl_high'] = new_wl[start:stop, 1]
        temp_sorted_data_idx = copy.deepcopy(self.param['sorted_data_idx'])
        self.param['sorted_data_idx'] = np.argsort(self.param['spectrum']['wl'])

        wl, model = forward(self.param, retrieval_mode=False)
        plt.plot(wl, model*1e6, color='#404784', label='MAP solution R=500')

        best_fit = np.array([wl, model]).T

        if os.path.isfile(self.param['out_dir'] + 'random_samples.dat'):
            fl = np.loadtxt(self.param['out_dir'] + 'random_samples.dat')
            plt.fill_between(fl[:, 0], (best_fit[:, 1] + (np.quantile(fl[:, 1:], [0.00135, 0.99865], axis=1)[1] - np.quantile(fl[:, 1:], 0.5, axis=1))) * 1e6, (best_fit[:, 1] + (np.quantile(fl[:, 1:], [0.00135, 0.99865], axis=1)[0] - np.quantile(fl[:, 1:], 0.5, axis=1))) * 1e6, ec=('#404784', 0.0), fc=('#404784', 0.25))
            plt.fill_between(fl[:, 0], (best_fit[:, 1] + (np.quantile(fl[:, 1:], [0.0225, 0.9775], axis=1)[1] - np.quantile(fl[:, 1:], 0.5, axis=1))) * 1e6, (best_fit[:, 1] + (np.quantile(fl[:, 1:], [0.0225, 0.9775], axis=1)[0] - np.quantile(fl[:, 1:], 0.5, axis=1))) * 1e6, ec=('#404784', 0.0), fc=('#404784', 0.5))
            plt.fill_between(fl[:, 0], (best_fit[:, 1] + (np.quantile(fl[:, 1:], [0.16, 0.84], axis=1)[1] - np.quantile(fl[:, 1:], 0.5, axis=1))) * 1e6, (best_fit[:, 1] + (np.quantile(fl[:, 1:], [0.16, 0.84], axis=1)[0] - np.quantile(fl[:, 1:], 0.5, axis=1))) * 1e6, ec=('#404784', 0.0), fc=('#404784', 0.75))

            best_fit = np.concatenate((best_fit, np.array([best_fit[:, 1] + (np.quantile(fl[:, 1:], [0.16, 0.84], axis=1)[1] - np.quantile(fl[:, 1:], 0.5, axis=1))]).T), axis=1)
            best_fit = np.concatenate((best_fit, np.array([best_fit[:, 1] + (np.quantile(fl[:, 1:], [0.16, 0.84], axis=1)[0] - np.quantile(fl[:, 1:], 0.5, axis=1))]).T), axis=1)
            best_fit = np.concatenate((best_fit, np.array([best_fit[:, 1] + (np.quantile(fl[:, 1:], [0.0225, 0.9775], axis=1)[1] - np.quantile(fl[:, 1:], 0.5, axis=1))]).T), axis=1)
            best_fit = np.concatenate((best_fit, np.array([best_fit[:, 1] + (np.quantile(fl[:, 1:], [0.0225, 0.9775], axis=1)[0] - np.quantile(fl[:, 1:], 0.5, axis=1))]).T), axis=1)
            best_fit = np.concatenate((best_fit, np.array([best_fit[:, 1] + (np.quantile(fl[:, 1:], [0.00135, 0.99865], axis=1)[1] - np.quantile(fl[:, 1:], 0.5, axis=1))]).T), axis=1)
            best_fit = np.concatenate((best_fit, np.array([best_fit[:, 1] + (np.quantile(fl[:, 1:], [0.00135, 0.99865], axis=1)[0] - np.quantile(fl[:, 1:], 0.5, axis=1))]).T), axis=1)

            del fl

        if solutions is None:
            np.savetxt(self.param['out_dir'] + 'Best_fit.dat', best_fit)
        else:
            np.savetxt(self.param['out_dir'] + 'Best_fit_(solution ' + str(solutions) + ').dat', best_fit)

        if self.param['spectrum']['bins']:
            self.param['spectrum']['wl'] = temp[:, 2]
            self.param['spectrum']['wl_low'] = temp[:, 0]
            self.param['spectrum']['wl_high'] = temp[:, 1]
        else:
            self.param['spectrum']['wl'] = temp + 0.0

        self.param['sorted_data_idx'] = copy.deepcopy(temp_sorted_data_idx)

        plt.legend()
        plt.xlabel('Wavelength [$\mu$m]')
        plt.ylabel('Transit depth (R$_p$/R$_{\star}$)$^2$ [ppm]')
        fig.tight_layout()

        if solutions is None:
            plt.savefig(self.param['out_dir'] + 'Nest_spectrum.pdf')
        else:
            plt.savefig(self.param['out_dir'] + 'Nest_spectrum (solution ' + str(solutions) + ').pdf')
        plt.close()

    def plot_P_profiles(self, solutions=None):
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
        ### CHEMISTRY ###
        fig, ax = plt.subplots()

        for mol in self.param['fit_molecules'] + [self.param['gas_fill']]:
            print(str(mol) + ' -> Top: ' + str(self.param['vmr_' + mol][0]) + ', Bottom: ' + str(self.param['vmr_' + mol][-1]))
            ax.loglog(self.param['vmr_' + mol], self.param['P'], label=mol)

        ax.set_xlim((1e-18, 1.5))
        ax.set_xlabel('Molecular VMR')
        ax.set_ylabel('Pressure [Pa]')

        def pa_to_bar(y):
            return y / (10. ** 5.)

        def bar_to_pa(y):
            return y * (10. ** 5.)

        if self.param['fit_wtr_cld']:
            wtr2 = self.param['vmr_H2O'] - np.roll(self.param['vmr_H2O'], 1)
            wtr2[0] = 0.0
            plt.hlines(self.param['P'][min(np.where(wtr2 != 0.0)[0]) - 1], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black', label='H$_2$O cloud')
            plt.hlines(self.param['P'][max(np.where(wtr2 != 0.0)[0])], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black')
        if self.param['fit_amn_cld']:
            wtr2 = self.param['vmr_NH3'] - np.roll(self.param['vmr_NH3'], 1)
            wtr2[0] = 0.0
            plt.hlines(self.param['P'][min(np.where(wtr2 != 0.0)[0]) - 1], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black', label='NH$_3$ cloud')
            plt.hlines(self.param['P'][max(np.where(wtr2 != 0.0)[0])], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black')
        if self.param['fit_gen_cld']:
            indx = find_nearest(self.param['P'], self.param['P_top'])
            plt.hlines(self.param['P'][indx], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black', label='Optically thick cloud')

        ax.yaxis.set_ticks(10. ** np.arange(np.log10(self.param['P'][0]), 9, 1))
        secax_y = ax.secondary_yaxis('right', functions=(pa_to_bar, bar_to_pa))
        ax.set_ylim((self.param['P'][0], 1e6))
        plt.gca().invert_yaxis()
        secax_y.set_ylabel('Pressure [bar]')

        ax.legend(loc='lower left')
        if solutions is None:
            plt.savefig(self.param['out_dir'] + 'Chemistry.pdf')
        else:
            plt.savefig(self.param['out_dir'] + 'Chemistry (solution ' + str(solutions) + ').pdf')
        plt.close()

        ### MEAN MOLECULAR MASS ###
        fig, ax = plt.subplots()
        ax.semilogy(self.param['mean_mol_weight'], self.param['P'])
        ax.set_xlabel('Mean molecular weight')
        ax.set_ylabel('Pressure [Pa]')

        ax.set_xlim((ax.get_xlim()[0], ax.get_xlim()[1]))

        if self.param['fit_wtr_cld']:
            wtr2 = self.param['vmr_H2O'] - np.roll(self.param['vmr_H2O'], 1)
            wtr2[0] = 0.0
            plt.hlines(self.param['P'][min(np.where(wtr2 != 0.0)[0]) - 1], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black', label='H$_2$O cloud')
            plt.hlines(self.param['P'][max(np.where(wtr2 != 0.0)[0])], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black')
        if self.param['fit_amn_cld']:
            wtr2 = self.param['vmr_NH3'] - np.roll(self.param['vmr_NH3'], 1)
            wtr2[0] = 0.0
            plt.hlines(self.param['P'][min(np.where(wtr2 != 0.0)[0]) - 1], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black', label='NH$_3$ cloud')
            plt.hlines(self.param['P'][max(np.where(wtr2 != 0.0)[0])], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black')
        if self.param['fit_gen_cld']:
            indx = find_nearest(self.param['P'], self.param['P_top'])
            plt.hlines(self.param['P'][indx], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black', label='Optically thick cloud')

        ax.yaxis.set_ticks(10. ** np.arange(np.log10(self.param['P'][0]), 9, 1))
        secax_y = ax.secondary_yaxis('right', functions=(pa_to_bar, bar_to_pa))
        ax.set_ylim((self.param['P'][0], 1e6))
        plt.gca().invert_yaxis()
        secax_y.set_ylabel('Pressure [bar]')
        if solutions is None:
            plt.savefig(self.param['out_dir'] + 'MMM.pdf')
        else:
            plt.savefig(self.param['out_dir'] + 'MMM (solution ' + str(solutions) + ').pdf')
        plt.close()

        ### PLANETARY TEMPERATURE ###
        fig, ax = plt.subplots()
        if self.param['TP_profile'] is None:
            ax.semilogy(np.ones(len(self.param['P'])) * self.param['Tp'], self.param['P'])
        else:
            if self.param['fit_T']:
                T = np.ones(len(self.param['P']))
                for i in range(0, len(T)):
                    T[i] = min(max(100.0, (self.param['TP_profile'](self.param['P']) + self.param['Tp'])[i]), 2000.0)
                ax.semilogy(T, self.param['P'])
            else:
                ax.semilogy(self.param['TP_profile'](self.param['P']), self.param['P'])

        ax.set_xlabel('Planetary Temperature [K]')
        ax.set_ylabel('Pressure [Pa]')

        if self.param['TP_profile'] is None:
            ax.set_xlim((ax.get_xlim()[0], ax.get_xlim()[1]))
        else:
            if self.param['fit_T']:
                ax.set_xlim((ax.get_xlim()[0], self.param['TP_profile'](1e6) + self.param['Tp']))
            else:
                ax.set_xlim((ax.get_xlim()[0], self.param['TP_profile'](1e6)))

        if self.param['fit_wtr_cld']:
            wtr2 = self.param['vmr_H2O'] - np.roll(self.param['vmr_H2O'], 1)
            wtr2[0] = 0.0
            plt.hlines(self.param['P'][min(np.where(wtr2 != 0.0)[0]) - 1], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black', label='H$_2$O cloud')
            plt.hlines(self.param['P'][max(np.where(wtr2 != 0.0)[0])], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black')
        if self.param['fit_amn_cld']:
            wtr2 = self.param['vmr_NH3'] - np.roll(self.param['vmr_NH3'], 1)
            wtr2[0] = 0.0
            plt.hlines(self.param['P'][min(np.where(wtr2 != 0.0)[0]) - 1], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black', label='NH$_3$ cloud')
            plt.hlines(self.param['P'][max(np.where(wtr2 != 0.0)[0])], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black')
        if self.param['fit_gen_cld']:
            indx = find_nearest(self.param['P'], self.param['P_top'])
            plt.hlines(self.param['P'][indx], ax.get_xlim()[0], ax.get_xlim()[1], linestyle='--', color='black', label='Optically thick cloud')

        ax.yaxis.set_ticks(10. ** np.arange(np.log10(self.param['P'][0]), 9, 1))
        secax_y = ax.secondary_yaxis('right', functions=(pa_to_bar, bar_to_pa))
        ax.set_ylim((self.param['P'][0], 1e6))
        plt.gca().invert_yaxis()
        secax_y.set_ylabel('Pressure [bar]')
        if not self.param['fit_T']:
            plt.savefig(self.param['out_dir'] + 'Tp.pdf')
        else:
            if solutions is None:
                plt.savefig(self.param['out_dir'] + 'Tp.pdf')
            else:
                plt.savefig(self.param['out_dir'] + 'Tp (solution ' + str(solutions) + ').pdf')
        plt.close()

    def plot_contribution(self, solutions=None):
        """
            Plots the various contributions to the model result (e.g., different gas species, Rayleigh scattering, cloud effects, etc.)

            Parameters
            ----------
            cube : array_like
                The input solution in the form of a cube (or parameter array).

            solutions : int, optional
                The specific solution number to be considered, if there are multiple solutions.

            Returns
            -------
            None. This function outputs a plot saved in the output directory.

            Notes
            -----
            The plot includes individual lines for each contribution, as well as the total model result and the observed data.
        """
        if not os.path.exists(self.param['out_dir'] + 'Contr_spec/'):
            os.mkdir(self.param['out_dir'] + 'Contr_spec/')

        fig = plt.figure(figsize=(12, 5))

        new_wl = np.loadtxt(self.param['pkg_dir'] + 'Data/wl_bins/bins_02_20_R500.dat')
        new_wl_central = np.mean(new_wl, axis=1)
        is_bins = self.param['spectrum']['bins']
        self.param['spectrum']['bins'] = False
        start = find_nearest(new_wl_central, min(self.param['spectrum']['wl']) - 0.05)
        stop = find_nearest(new_wl_central, max(self.param['spectrum']['wl']) + 0.05)

        temp_wl = self.param['spectrum']['wl'] + 0.0
        self.param['spectrum']['wl'] = new_wl_central[start:stop]

        temp_sorted_data_idx = copy.deepcopy(self.param['sorted_data_idx'])
        self.param['sorted_data_idx'] = np.argsort(self.param['spectrum']['wl'])

        for mol in self.param['fit_molecules'] + [self.param['gas_fill']]:
            self.param[mol + '_contribution'] = False
        self.param['cld_contribution'] = False
        self.param['Rayleigh_contribution'] = False
        self.param['haze_contribution'] = False
        self.param['star_act_contribution'] = False

        if self.param['gas_fill'] == 'H2':
            print('Plotting the H2-H2 CIA contribution')
            self.param['CIA_contribution'] = True
            wl, model = forward(self.param, retrieval_mode=False)
            comp = np.array([wl, model]).T
            np.savetxt(self.param['out_dir'] + 'Contr_spec/contr_CIA.dat', comp)
            self.param['CIA_contribution'] = False

            plt.plot(wl, model * 1e6, linestyle='--', label='H$_2$-H$_2$ CIA')

        print('Plotting the Rayleigh scattering contribution')
        self.param['Rayleigh_contribution'] = True
        wl, model = forward(self.param, retrieval_mode=False)
        comp = np.array([wl, model]).T
        np.savetxt(self.param['out_dir'] + 'Contr_spec/contr_Rayleigh.dat', comp)
        self.param['Rayleigh_contribution'] = False

        plt.plot(wl, model * 1e6, linestyle='--', label='Rayleigh scattering')

        if self.param['incl_clouds']:
            print('Plotting the contribution of cloud')
            self.param['cld_contribution'] = True
            wl, model = forward(self.param, retrieval_mode=False)
            comp = np.array([wl, model]).T
            np.savetxt(self.param['out_dir'] + 'Contr_spec/contr_cloud.dat', comp)
            self.param['cld_contribution'] = False

            plt.plot(wl, model * 1e6, linestyle='--', label='Cloud')

        if self.param['incl_haze']:
            print('Plotting the contribution of haze')
            self.param['haze_contribution'] = True
            wl, model = forward(self.param, retrieval_mode=False)
            comp = np.array([wl, model]).T
            np.savetxt(self.param['out_dir'] + 'Contr_spec/contr_haze.dat', comp)
            self.param['haze_contribution'] = False

            plt.plot(wl, model * 1e6, linestyle='--', label='Haze')

        if self.param['incl_star_activity']:
            print('Plotting the contribution of the star activity')
            self.param['star_act_contribution'] = True
            wl, model = forward(self.param, retrieval_mode=False)
            comp = np.array([wl, model]).T
            np.savetxt(self.param['out_dir'] + 'Contr_spec/contr_star.dat', comp)
            self.param['star_act_contribution'] = False

            plt.plot(wl, model * 1e6, linestyle='--', label='Star heterogeneity')

        if not self.param['bare_rock']:
            print('Plotting the contribution of atmospheric gases')
        for mol in self.param['fit_molecules'] + [self.param['gas_fill']]:
            self.param[mol + '_contribution'] = True
            wl, model = forward(self.param, retrieval_mode=False)
            comp = np.array([wl, model]).T
            np.savetxt(self.param['out_dir'] + 'Contr_spec/contr_' + mol + '.dat', comp)
            self.param[mol + '_contribution'] = False

            plt.plot(wl, model * 1e6, linewidth=0.5, label=mol)

        for mol in self.param['fit_molecules'] + [self.param['gas_fill']]:
            self.param[mol + '_contribution'] = True
        self.param['cld_contribution'] = True
        self.param['Rayleigh_contribution'] = True
        self.param['CIA_contribution'] = True
        self.param['haze_contribution'] = True
        self.param['star_act_contribution'] = True

        wl, model = forward(self.param, retrieval_mode=False)

        plt.plot(wl, model * 1e6, color='black', label='MAP solution R=500')

        self.param['spectrum']['wl'] = temp_wl + 0.0
        self.param['sorted_data_idx'] = copy.deepcopy(temp_sorted_data_idx)

        plt.errorbar(self.param['spectrum']['wl'], self.param['spectrum']['T_depth'] * 1e6, yerr=self.param['spectrum']['error_T'] * 1e6,
                     linestyle='', linewidth=0.5, color='black', marker='o', markerfacecolor='red', markersize=4, capsize=1.75, label='Data')

        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.xlabel('Wavelength [$\mu$m]')
        plt.ylabel('Transit depth (R$_p$/R$_{\star}$)$^2$')
        fig.tight_layout()

        if solutions is None:
            plt.savefig(self.param['out_dir'] + 'Contribution.pdf')
        else:
            plt.savefig(self.param['out_dir'] + 'Contribution (solution ' + str(solutions) + ').pdf')
        plt.close()

        if is_bins:
            self.param['spectrum']['bins'] = True

    def calc_spectra(self, mc_samples):
        new_wl = np.loadtxt(self.param['pkg_dir'] + 'Data/wl_bins/bins_02_20_R500.dat')
        new_wl_central = np.mean(new_wl, axis=1)
        start = find_nearest(new_wl_central, min(self.param['spectrum']['wl']) - 0.05)
        stop = find_nearest(new_wl_central, max(self.param['spectrum']['wl']) + 0.05)
        wl_len = len(self.param['spectrum']['wl'])
        if self.param['spectrum']['bins']:
            temp = np.array([self.param['spectrum']['wl_low'], self.param['spectrum']['wl_high'], self.param['spectrum']['wl']]).T
        else:
            temp = self.param['spectrum']['wl'] + 0.0
        self.param['spectrum']['wl'] = new_wl_central[start:stop] + 0.0
        self.param['spectrum']['wl_low'] = new_wl[start:stop, 0] + 0.0
        self.param['spectrum']['wl_high'] = new_wl[start:stop, 1] + 0.0
        temp_sorted_data_idx = copy.deepcopy(self.param['sorted_data_idx'])
        self.param['sorted_data_idx'] = np.argsort(self.param['spectrum']['wl'])

        if mc_samples.shape[0] < self.param['n_likelihood_data']:
            self.param['n_likelihood_data'] = mc_samples.shape[0] - MPIsize
        else:
            pass

        samples = np.zeros((len(self.param['spectrum']['wl']), int(self.param['n_likelihood_data'] / MPIsize) + 1))
        samples[:, 0] = self.param['spectrum']['wl']
        loglike_data = np.zeros((int(self.param['n_likelihood_data'] / MPIsize), wl_len))

        if MPIrank == 0:
            print('\nCalculating the likelihood per data point')
            try:
                os.mkdir(self.param['out_dir'] + 'loglikelihood_per_datapoint/')
            except OSError:
                pass

        idx = np.random.choice(mc_samples.shape[0], int(self.param['n_likelihood_data']), replace=False)

        for i in range(int(self.param['n_likelihood_data'] / MPIsize)):
            out_cube = mc_samples[idx[i], :]
            self.cube_to_param(out_cube)

            _, samples[:, i + 1] = forward(self.param, retrieval_mode=False)

            if self.param['spectrum']['bins']:
                model = spectres(temp[temp_sorted_data_idx, 2], self.param['spectrum']['wl'], samples[:, i + 1], fill=False)
            else:
                model = spectres(temp[temp_sorted_data_idx], self.param['spectrum']['wl'], samples[:, i + 1], fill=False)

            # Calculate likelihood per single datapoint
            chi = (self.param['spectrum']['T_depth'][temp_sorted_data_idx] - model) / self.param['spectrum']['error_T'][temp_sorted_data_idx]
            loglike_data[i, :] = ((-1.) * np.log(self.param['spectrum']['error_T'][temp_sorted_data_idx] * np.sqrt(2.0 * math.pi))) - (0.5 * chi * chi)

        np.savetxt(self.param['out_dir'] + 'loglikelihood_per_datapoint/loglike_' + str(MPIrank) + '.dat', loglike_data)
        np.savetxt(self.param['out_dir'] + 'loglikelihood_per_datapoint/samples_' + str(MPIrank) + '.dat', samples)

        if self.param['spectrum']['bins']:
            self.param['spectrum']['wl'] = temp[:, 2] + 0.0
            self.param['spectrum']['wl_low'] = temp[:, 0] + 0.0
            self.param['spectrum']['wl_high'] = temp[:, 1] + 0.0
        else:
            self.param['spectrum']['wl'] = temp + 0.0

        self.param['sorted_data_idx'] = copy.deepcopy(temp_sorted_data_idx)

    def plot_posteriors(self, prefix, multinest_results, parameters, mds):
        """
            Generates and saves posterior distribution plots based on the given parameters and results of multinest analysis.

            Args:
                self: The instance of the class where this function is defined.
                prefix (str): The prefix string used in naming files.
                multinest_results: The results obtained from a multinest analysis.
                parameters: The parameters used in the analysis.
                mds (int): The number of modes to use for a multimodal solution.

            Returns:
                None. The function saves generated plots as PDFs and PNGs in the current directory.

            Notes:
                - The function also creates, moves, and deletes files in the current directory.
                - The function expects certain files to exist in the current directory with specific naming patterns based on the provided prefix.
                - If mds is less than 2, the function processes the results as a unimodal solution, otherwise it processes them as a multimodal solution.
        """
        from numpy import log
        from six.moves import range
        import logging
        import types
        from matplotlib.ticker import MaxNLocator, NullLocator
        from matplotlib.colors import LinearSegmentedColormap, colorConverter
        from matplotlib.ticker import ScalarFormatter
        from scipy.ndimage import gaussian_filter as norm_kde
        from scipy.stats import gaussian_kde

        try:
            str_type = types.StringTypes
            float_type = types.FloatType
            int_type = types.IntType
        except:
            str_type = str
            float_type = float
            int_type = int

        SQRTEPS = math.sqrt(float(np.finfo(np.float64).eps))

        def _posteriors_clr_to_vmr(prefix, modes=None):
            """
                Transforms the solutions from centered log-ratio (clr) to volume mixing ratio (vmr) and saves the results.

                Args:
                    param (dict): A dictionary that includes parameters of the system.
                    prefix (str): The prefix of the file names to be handled.
                    modes (int, optional): The mode index. If provided, the function will handle the solution file for the specific mode.
                                           Default is None.

                The function does not return anything.

                The function performs several steps:
                1. Copies the solution file to a new file with '_original' appended to the file name.
                2. Loads the solution data.
                3. Checks if clouds are included and adjusts data accordingly.
                4. Converts centered log-ratio to volume mixing ratio.
                5. Calculates the mean molar mass.
                6. Writes the transformed data back to the solution file.
                7. If the 'params.json' file exists, it leaves it as is. Otherwise, it renames 'params.json' to 'params_original.json'.
                8. Writes the parameter names to a new 'params.json' file.

                Please note that this function assumes the existence of specific files and writes to these files.
                Ensure that these files exist and that the file paths are accessible.
            """
            if modes is not None:
                os.system('cp ' + prefix + 'solution' + str(modes) + '.txt ' + prefix + 'solution' + str(modes) + '_original.txt')
                a = np.loadtxt(prefix + 'solution' + str(modes) + '.txt')
            else:
                os.system('cp ' + prefix + '.txt ' + prefix + 'original.txt')
                a = np.loadtxt(prefix + '.txt')

            if self.param['bare_rock']:
                b = np.ones((len(a[:, 0]), len(a[0, :])))
            elif len(self.param['fit_molecules']) < 1:
                b = np.ones((len(a[:, 0]), len(a[0, :]) + 1))
            else:
                b = np.ones((len(a[:, 0]), len(a[0, :]) + 2))

            if self.param['incl_clouds'] and self.param['fit_gen_cld']:
                z = 4
            elif self.param['incl_clouds'] and not self.param['fit_gen_cld']:
                z = 6
            else:
                z = 3
            if self.param['incl_haze']:
                z += 2
            if self.param['fit_offset']:
                z += self.param['n_offsets']
            if self.param['fit_Mp']:
                z += 1
            if self.param['fit_T']:
                z += 1

            b[:, 0:z] = a[:, 0:z] + 0.0

            if self.param['fit_offset']:
                for i in range(0, self.param['n_offsets']):
                    b[:, 2 + i] *= 1e6             # offsets in ppm

            if not self.param['bare_rock']:
                volume_mixing_ratio = {}
                if len(self.param['fit_molecules']) < 1:
                    i = -1
                    volume_mixing_ratio[self.param['gas_fill']] = 1.0
                    mmm = volume_mixing_ratio[self.param['gas_fill']] * self.param['mm'][self.param['gas_fill']]
                    b[:, z + i + 1] = np.array(mmm) + 0.0
                else:
                    centered_log_ratio = {}
                    mol_indx = z + 0
                    for mol in self.param['fit_molecules']:
                        centered_log_ratio[mol] = np.array(a[:, mol_indx])
                        mol_indx += 1
                    centered_log_ratio[self.param['gas_fill']] = np.zeros(len(a[:, 0]))

                    sumb = np.zeros(len(a[:, 0]))
                    for mol in self.param['fit_molecules']:
                        centered_log_ratio[self.param['gas_fill']] -= centered_log_ratio[mol]
                        sumb += np.exp(centered_log_ratio[mol])
                    sumb += np.exp(centered_log_ratio[self.param['gas_fill']])

                    for mol in centered_log_ratio.keys():
                        volume_mixing_ratio[mol] = np.exp(centered_log_ratio[mol]) / sumb
                    del centered_log_ratio, sumb

                    mmm = np.zeros(len(a[:, 0]))
                    for mol in volume_mixing_ratio.keys():
                        mmm += volume_mixing_ratio[mol] * self.param['mm'][mol]

                    for i, mol in enumerate(self.param['fit_molecules']):
                        b[:, z + i] = np.log10(volume_mixing_ratio[mol])

                    b[:, z + i + 1] = np.log10(volume_mixing_ratio[self.param['gas_fill']])
                    b[:, z + i + 2] = np.array(mmm) + 0.0

            if self.param['incl_star_activity']:
                if self.param['stellar_activity_parameters'] == int(3):
                    num = -3
                elif self.param['stellar_activity_parameters'] == int(5):
                    num = -5
                b[:, num:] = a[:, num:]

            if modes is not None:
                np.savetxt(prefix + 'solution' + str(modes) + '.txt', b)
            else:
                np.savetxt(prefix + '.txt', b)

            return z

        def _corner_parameters():
            if os.path.isfile(prefix + 'params_original.json'):
                pass
            else:
                os.system('mv ' + prefix + 'params.json ' + prefix + 'params_original.json')

            par = []
            if self.param['fit_offset']:
                for i in range(0, self.param['n_offsets']):
                    par.append("offset$_" + str(i + 1) + "$ [ppm]")
            par.append("R$_p$ [R$_{\oplus}$]")
            if self.param['fit_Mp']:
                par.append("M$_p$ [M$_{\oplus}$]")
            if self.param['fit_T'] and self.param['TP_profile'] is None:
                par.append("T$_p$ [K]")
            elif self.param['fit_T'] and self.param['TP_profile'] is not None:
                par.append("$\Delta$T$_p$ [K]")
            if self.param['fit_wtr_cld']:
                par.append("Log(P$_{top, H_2O}$)")
                par.append("Log(D$_{H_2O}$)")
                par.append("Log(CR$_{H_2O}$)")
            if self.param['fit_amn_cld']:
                par.append("Log(P$_{top, NH_3}$)")
                par.append("Log(D$_{NH_3}$)")
                par.append("Log(CR$_{NH_3}$)")
            if self.param['fit_gen_cld']:
                par.append("Log(P$_{top}$)")
            if self.param['incl_haze']:
                par.append("Log(d$_{haze}$)")
                par.append("Log(haze)")
            if not self.param['bare_rock']:
                for mol in self.param['fit_molecules']:
                    par.append(self.param['formatted_labels'][mol])
                if len(self.param['fit_molecules']) >= 1:
                    par.append(self.param['formatted_labels'][self.param['gas_fill']] + " (derived)")
                par.append("$\mu$ (derived)")
            if self.param['incl_star_activity']:
                if self.param['stellar_activity_parameters'] == int(3):
                    par.append("$\delta$" + "$_{het}$")
                    par.append("$T_{het}$")
                    par.append("$T_{phot}$")
                elif self.param['stellar_activity_parameters'] == int(5):
                    par.append("$\delta$" + "$_{spot}$")
                    par.append("$\delta$" + "$_{fac}$")
                    par.append("$T_{spot}$")
                    par.append("$T_{fac}$")
                    par.append("$T_{phot}$")
            json.dump(par, open(prefix + 'params.json', 'w'))

        def _quantile(x, q, weights=None):
            """
            Compute (weighted) quantiles from an input set of samples.
            Parameters
            ----------
            x : `~numpy.ndarray` with shape (nsamps,)
                Input samples.
            q : `~numpy.ndarray` with shape (nquantiles,)
               The list of quantiles to compute from `[0., 1.]`.
            weights : `~numpy.ndarray` with shape (nsamps,), optional
                The associated weight from each sample.
            Returns
            -------
            quantiles : `~numpy.ndarray` with shape (nquantiles,)
                The weighted sample quantiles computed at `q`.
            """

            # Initial check.
            x = np.atleast_1d(x)
            q = np.atleast_1d(q)

            # Quantile check.
            if np.any(q < 0.0) or np.any(q > 1.0):
                raise ValueError("Quantiles must be between 0. and 1.")

            if weights is None:
                # If no weights provided, this simply calls `np.percentile`.
                return np.percentile(x, list(100.0 * q))
            else:
                # If weights are provided, compute the weighted quantiles.
                weights = np.atleast_1d(weights)
                if len(x) != len(weights):
                    raise ValueError("Dimension mismatch: len(weights) != len(x).")
                idx = np.argsort(x)  # sort samples
                sw = weights[idx]  # sort weights
                cdf = np.cumsum(sw)[:-1]  # compute CDF
                cdf /= cdf[-1]  # normalize CDF
                cdf = np.append(0, cdf)  # ensure proper span
                quantiles = np.interp(q, cdf, x[idx]).tolist()
                return quantiles

        def _store_nest_solutions():
            """
                Stores the solutions of a nested sampling calculation and writes them to a pickle file.

                This function does not take any input arguments.

                Returns:
                NEST_out (dict): A dictionary containing the solutions of the nested sampling calculation.
                                  This includes nested sampling statistics, evidence, modes, weights, likelihoods, traces,
                                  and other fitting parameters.

                The function performs several steps:
                1. Loads data from a '.txt' file and fetches statistics from a multinest results object.
                2. Checks if the computation was multimodal, and if so, separates the different modes.
                3. Stores the modes, weights, and likelihoods in separate lists.
                4. Constructs a dictionary of results (mydict), which contains fit parameters and their associated statistics for each mode.
                5. If there is more than one solution, the function writes each solution to a separate '.txt' file.
                6. Finally, the function writes all solutions to a pickle file.

                Please note that the filename prefix is not provided as an argument to this function.
                This implies that the filename prefix is a global variable accessible to this function.
            """
            NEST_out = {'solutions': {}}
            data = np.loadtxt(prefix + '.txt')
            NEST_stats = multinest_results.get_stats()
            NEST_out['NEST_stats'] = NEST_stats
            NEST_out['global_logE'] = (NEST_out['NEST_stats']['global evidence'], NEST_out['NEST_stats']['global evidence error'])

            modes = []
            modes_weights = []
            modes_loglike = []
            chains = []
            chains_weights = []
            chains_loglike = []

            if self.param['multimodal'] and mds > 1:
                # separate modes. get individual samples for each mode
                # get parameter values and sample probability (=weight) for each mode
                with open(prefix + 'post_separate.dat') as f:
                    lines = f.readlines()
                    for idx, line in enumerate(lines):
                        if idx > 2:  # skip the first two lines
                            if lines[idx - 1] == '\n' and lines[idx - 2] == '\n':
                                modes.append(chains)
                                modes_weights.append(chains_weights)
                                modes_loglike.append(chains_loglike)
                                chains = []
                                chains_weights = []
                                chains_loglike = []
                        chain = [float(x) for x in line.split()[2:]]
                        if len(chain) > 0:
                            chains.append(chain)
                            chains_weights.append(float(line.split()[0]))
                            chains_loglike.append(float(line.split()[1]))
                    modes.append(chains)
                    modes_weights.append(chains_weights)
                    modes_loglike.append(chains_loglike)
                modes_array = []
                for mode in modes:
                    mode_array = np.zeros((len(mode), len(mode[0])))
                    for idx, line in enumerate(mode):
                        mode_array[idx, :] = line
                    modes_array.append(mode_array)
            else:
                # not running in multimode. Get chains directly from file prefix.txt
                modes_array = [data[:, 2:]]
                chains_weights = [data[:, 0]]
                modes_weights.append(chains_weights[0])
                chains_loglike = [data[:, 1]]
                modes_loglike.append(chains_loglike[0])
                modes = [0]

            for nmode in range(len(modes)):
                mydict = {'type': 'nest',
                          'local_logE': (NEST_out['NEST_stats']['modes'][nmode]['local log-evidence'], NEST_out['NEST_stats']['modes'][nmode]['local log-evidence error']),
                          'weights': np.asarray(modes_weights[nmode]),
                          'loglike': np.asarray(modes_loglike[nmode]),
                          'tracedata': modes_array[nmode],
                          'fit_params': {}}

                for idx, param_name in enumerate(parameters):
                    trace = modes_array[nmode][:, idx]
                    q_16, q_50, q_84 = _quantile(trace, [0.16, 0.5, 0.84], weights=np.asarray(modes_weights[nmode]))
                    mydict['fit_params'][param_name] = {
                        'value': q_50,
                        'sigma_m': q_50 - q_16,
                        'sigma_p': q_84 - q_50,
                        'nest_map': NEST_stats['modes'][nmode]['maximum a posterior'][idx],
                        'mean': NEST_stats['modes'][nmode]['mean'][idx],
                        'nest_sigma': NEST_stats['modes'][nmode]['sigma'][idx],
                        'trace': trace,
                    }

                NEST_out['solutions']['solution{}'.format(nmode)] = mydict

            if len(NEST_out['solutions']) > 1:
                for i in range(len(NEST_out['solutions'])):
                    fl = np.ones((len(NEST_out['solutions']['solution' + str(i)]['weights']), len(NEST_out['solutions']['solution' + str(i)]['tracedata'][0, :]) + 2))
                    fl[:, 0] = NEST_out['solutions']['solution' + str(i)]['weights']
                    fl[:, 1] = NEST_out['solutions']['solution' + str(i)]['loglike']
                    fl[:, 2:] = NEST_out['solutions']['solution' + str(i)]['tracedata']
                    np.savetxt(prefix + 'solution' + str(i) + '.txt', fl)

            return NEST_out

        def _plotting_bounds(results, gas_pos, modes=None):
            """
                Computes the boundaries for the plot based on the results from a dynamic nested sampling calculation.

                Parameters:
                results (dict): A dictionary containing the results of a dynamic nested sampling calculation.
                                This should include 'samples', 'logwt', and 'logz'.
                modes (list, optional): A list of modes to be considered for the computation of boundaries.
                                        Default is None, in which case all modes are considered.

                Returns:
                boundaries (list): A list containing the computed boundaries for the plot.

                The function performs several steps:
                1. Extract samples and weights from the results. If 'logwt' and 'logz' exist in results, use these to compute weights.
                2. Ensures that samples are 2D and transposed such that samples.shape[0] <= samples.shape[1].
                3. Checks the dimensionality and length of weights to ensure they are 1-D and have the same number of elements as samples.
                4. Defines an initial list of boundaries, then iterates over them, updating each boundary value
                   based on the quantile computed from samples and weights.
            """
            if modes is None:
                samples = results['samples']
                try:
                    weights = np.exp(results['logwt'] - results['logz'][-1])
                except:
                    weights = results['weights']

                # Deal with 1D results. A number of extra catches are also here
                # in case users are trying to plot other results besides the `Results`
                # instance generated by `dynesty`.
                samples = np.atleast_1d(samples)
                if len(samples.shape) == 1:
                    samples = np.atleast_2d(samples)
                else:
                    assert len(samples.shape) == 2, "Samples must be 1- or 2-D."
                    samples = samples.T
                assert samples.shape[0] <= samples.shape[1], "There are more " \
                                                             "dimensions than samples!"
                ndim, nsamps = samples.shape

                # Check weights.
                if weights.ndim != 1:
                    raise ValueError("Weights must be 1-D.")
                if nsamps != weights.shape[0]:
                    raise ValueError("The number of weights and samples disagree!")

                boundaries = [0.999999426697 for _ in range(len(results['samples'][0, :]))]
                boundaries = list(boundaries)
                for i, _ in enumerate(boundaries):
                    q = [0.5 - 0.5 * boundaries[i], 0.5 + 0.5 * boundaries[i]]
                    boundaries[i] = _quantile(samples[i], q, weights=weights)
                    if gas_pos - 2 <= i < gas_pos + len(self.param['fit_molecules']) and boundaries[i][0] < -12.0:
                        boundaries[i][0] = -12.0

                return boundaries
            else:
                boundaries = {}
                for sol in range(0, modes):
                    # Extract weighted samples.
                    samples = results[str(sol)]['samples']
                    try:
                        weights = np.exp(results[str(sol)]['logwt'] - results[str(sol)]['logz'][-1])
                    except:
                        weights = results[str(sol)]['weights']

                    # Deal with 1D results. A number of extra catches are also here
                    # in case users are trying to plot other results besides the `Results`
                    # instance generated by `dynesty`.
                    samples = np.atleast_1d(samples)
                    if len(samples.shape) == 1:
                        samples = np.atleast_2d(samples)
                    else:
                        assert len(samples.shape) == 2, "Samples must be 1- or 2-D."
                        samples = samples.T
                    assert samples.shape[0] <= samples.shape[1], "There are more " \
                                                                 "dimensions than samples!"
                    ndim, nsamps = samples.shape

                    # Check weights.
                    if weights.ndim != 1:
                        raise ValueError("Weights must be 1-D.")
                    if nsamps != weights.shape[0]:
                        raise ValueError("The number of weights and samples disagree!")

                    boundaries[str(sol)] = [0.999999426697 for _ in range(len(results[str(sol)]['samples'][0, :]))]
                    boundaries[str(sol)] = list(boundaries[str(sol)])
                    for i, _ in enumerate(boundaries[str(sol)]):
                        q = [0.5 - 0.5 * boundaries[str(sol)][i], 0.5 + 0.5 * boundaries[str(sol)][i]]
                        boundaries[str(sol)][i] = _quantile(samples[i], q, weights=weights)

                bound = []
                for i in range(ndim):
                    min_b, max_b = [], []
                    for j in boundaries.keys():
                        if gas_pos - 2 <= i < gas_pos + len(self.param['fit_molecules']) and boundaries[j][i][0] < -12.0:
                            min_b.append(-12.0)
                        else:
                            min_b.append(boundaries[j][i][0])
                        max_b.append(boundaries[j][i][1])
                    bound.append([min(min_b), max(max_b)])

                return list(bound)

        def _resample_equal(samples, weights, rstate=None):
            """
            Resample a new set of points from the weighted set of inputs
            such that they all have equal weight.
            Each input sample appears in the output array either
            `floor(weights[i] * nsamples)` or `ceil(weights[i] * nsamples)` times,
            with `floor` or `ceil` randomly selected (weighted by proximity).
            Parameters
            ----------
            samples : `~numpy.ndarray` with shape (nsamples,)
                Set of unequally weighted samples.
            weights : `~numpy.ndarray` with shape (nsamples,)
                Corresponding weight of each sample.
            rstate : `~numpy.random.RandomState`, optional
                `~numpy.random.RandomState` instance.
            Returns
            -------
            equal_weight_samples : `~numpy.ndarray` with shape (nsamples,)
                New set of samples with equal weights.
            Notes
            -----
            Implements the systematic resampling method described in `Hol, Schon, and
            Gustafsson (2006) <doi:10.1109/NSSPW.2006.4378824>`_.
            """

            if rstate is None:
                rstate = np.random

            if abs(np.sum(weights) - 1.) > SQRTEPS:  # same tol as in np.random.choice.
                raise ValueError("Weights do not sum to 1.")

            # Make N subdivisions and choose positions with a consistent random offset.
            nsamples = len(weights)
            positions = (rstate.random() + np.arange(nsamples)) / nsamples

            # Resample the data.
            idx = np.zeros(nsamples, dtype=int)
            cumulative_sum = np.cumsum(weights)
            i, j = 0, 0
            while i < nsamples:
                if positions[i] < cumulative_sum[j]:
                    idx[i] = j
                    i += 1
                else:
                    j += 1

            return samples[idx]

        def _hist2d(x, y, smooth=0.02, span=None, weights=None, levels=None,
                    ax=None, color='gray', plot_datapoints=False, plot_density=True,
                    plot_contours=True, no_fill_contours=False, fill_contours=True,
                    contour_kwargs=None, contourf_kwargs=None, data_kwargs=None,
                    **kwargs):
            """
            Internal function called by :meth:`cornerplot` used to generate
            a 2-D histogram/contour of samples.

            Parameters
            ----------
            x : interable with shape (nsamps,)
               Sample positions in the first dimension.

            y : iterable with shape (nsamps,)
               Sample positions in the second dimension.

            span : iterable with shape (ndim,), optional
                A list where each element is either a length-2 tuple containing
                lower and upper bounds or a float from `(0., 1.]` giving the
                fraction of (weighted) samples to include. If a fraction is provided,
                the bounds are chosen to be equal-tailed. An example would be::

                    span = [(0., 10.), 0.95, (5., 6.)]

                Default is `0.999999426697` (5-sigma credible interval).

            weights : iterable with shape (nsamps,)
                Weights associated with the samples. Default is `None` (no weights).

            levels : iterable, optional
                The contour levels to draw. Default are `[0.5, 1, 1.5, 2]`-sigma.

            ax : `~matplotlib.axes.Axes`, optional
                An `~matplotlib.axes.axes` instance on which to add the 2-D histogram.
                If not provided, a figure will be generated.

            color : str, optional
                The `~matplotlib`-style color used to draw lines and color cells
                and contours. Default is `'gray'`.

            plot_datapoints : bool, optional
                Whether to plot the individual data points. Default is `False`.

            plot_density : bool, optional
                Whether to draw the density colormap. Default is `True`.

            plot_contours : bool, optional
                Whether to draw the contours. Default is `True`.

            no_fill_contours : bool, optional
                Whether to add absolutely no filling to the contours. This differs
                from `fill_contours=False`, which still adds a white fill at the
                densest points. Default is `False`.

            fill_contours : bool, optional
                Whether to fill the contours. Default is `True`.

            contour_kwargs : dict
                Any additional keyword arguments to pass to the `contour` method.

            contourf_kwargs : dict
                Any additional keyword arguments to pass to the `contourf` method.

            data_kwargs : dict
                Any additional keyword arguments to pass to the `plot` method when
                adding the individual data points.

            """

            if ax is None:
                ax = plt.gca()

            # Determine plotting bounds.
            data = [x, y]
            if span is None:
                span = [0.999999426697 for i in range(2)]
            span = list(span)
            if len(span) != 2:
                raise ValueError("Dimension mismatch between samples and span.")
            for i, _ in enumerate(span):
                try:
                    xmin, xmax = span[i]
                except:
                    q = [0.5 - 0.5 * span[i], 0.5 + 0.5 * span[i]]
                    span[i] = _quantile(data[i], q, weights=weights)

            # The default "sigma" contour levels.
            if levels is None:
                levels = 1.0 - np.exp(-0.5 * np.arange(0.5, 2.1, 0.5) ** 2)

            # Color map for the density plot, over-plotted to indicate the
            # density of the points near the center.
            density_cmap = LinearSegmentedColormap.from_list(
                "density_cmap", [color, (1, 1, 1, 0)])

            # Color map used to hide the points at the high density areas.
            white_cmap = LinearSegmentedColormap.from_list(
                "white_cmap", [(1, 1, 1), (1, 1, 1)], N=2)

            # This "color map" is the list of colors for the contour levels if the
            # contours are filled.
            rgba_color = colorConverter.to_rgba(color)
            contour_cmap = [list(rgba_color) for l in levels] + [rgba_color]
            for i, l in enumerate(levels):
                contour_cmap[i][-1] *= float(i) / (len(levels) + 1)

            # Initialize smoothing.
            if (isinstance(smooth, int_type) or isinstance(smooth, float_type)):
                smooth = [smooth, smooth]
            bins = []
            svalues = []
            for s in smooth:
                if isinstance(s, int_type):
                    # If `s` is an integer, the weighted histogram has
                    # `s` bins within the provided bounds.
                    bins.append(s)
                    svalues.append(0.)
                else:
                    # If `s` is a float, oversample the data relative to the
                    # smoothing filter by a factor of 2, then use a Gaussian
                    # filter to smooth the results.
                    bins.append(int(round(2. / s)))
                    svalues.append(2.)

            # We'll make the 2D histogram to directly estimate the density.
            try:
                H, X, Y = np.histogram2d(x.flatten(), y.flatten(), bins=bins,
                                         range=list(map(np.sort, span)),
                                         weights=weights)
            except ValueError:
                raise ValueError("It looks like at least one of your sample columns "
                                 "have no dynamic range.")

            # Smooth the results.
            if not np.all(svalues == 0.):
                H = norm_kde(H, svalues)

            # Compute the density levels.
            Hflat = H.flatten()
            inds = np.argsort(Hflat)[::-1]
            Hflat = Hflat[inds]
            sm = np.cumsum(Hflat)
            sm /= sm[-1]
            V = np.empty(len(levels))
            for i, v0 in enumerate(levels):
                try:
                    V[i] = Hflat[sm <= v0][-1]
                except:
                    V[i] = Hflat[0]
            V.sort()
            m = (np.diff(V) == 0)
            if np.any(m) and plot_contours:
                logging.warning("Too few points to create valid contours.")
            while np.any(m):
                V[np.where(m)[0][0]] *= 1.0 - 1e-4
                m = (np.diff(V) == 0)
            V.sort()

            # Compute the bin centers.
            X1, Y1 = 0.5 * (X[1:] + X[:-1]), 0.5 * (Y[1:] + Y[:-1])

            # Extend the array for the sake of the contours at the plot edges.
            H2 = H.min() + np.zeros((H.shape[0] + 4, H.shape[1] + 4))
            H2[2:-2, 2:-2] = H
            H2[2:-2, 1] = H[:, 0]
            H2[2:-2, -2] = H[:, -1]
            H2[1, 2:-2] = H[0]
            H2[-2, 2:-2] = H[-1]
            H2[1, 1] = H[0, 0]
            H2[1, -2] = H[0, -1]
            H2[-2, 1] = H[-1, 0]
            H2[-2, -2] = H[-1, -1]
            X2 = np.concatenate([X1[0] + np.array([-2, -1]) * np.diff(X1[:2]), X1,
                                 X1[-1] + np.array([1, 2]) * np.diff(X1[-2:])])
            Y2 = np.concatenate([Y1[0] + np.array([-2, -1]) * np.diff(Y1[:2]), Y1,
                                 Y1[-1] + np.array([1, 2]) * np.diff(Y1[-2:])])

            # Plot the data points.
            if plot_datapoints:
                if data_kwargs is None:
                    data_kwargs = dict()
                data_kwargs["color"] = data_kwargs.get("color", color)
                data_kwargs["ms"] = data_kwargs.get("ms", 2.0)
                data_kwargs["mec"] = data_kwargs.get("mec", "none")
                data_kwargs["alpha"] = data_kwargs.get("alpha", 0.1)
                ax.plot(x, y, "o", zorder=-1, rasterized=True, **data_kwargs)

            # Plot the base fill to hide the densest data points.
            if (plot_contours or plot_density) and not no_fill_contours:
                ax.contourf(X2, Y2, H2.T, [V.min(), H.max()],
                            cmap=white_cmap, antialiased=False)

            if plot_contours and fill_contours:
                if contourf_kwargs is None:
                    contourf_kwargs = dict()
                contourf_kwargs["colors"] = contourf_kwargs.get("colors", contour_cmap)
                contourf_kwargs["antialiased"] = contourf_kwargs.get("antialiased",
                                                                     False)
                ax.contourf(X2, Y2, H2.T, np.concatenate([[0], V, [H.max() * (1 + 1e-4)]]),
                            **contourf_kwargs)

            # Plot the density map. This can't be plotted at the same time as the
            # contour fills.
            elif plot_density:
                ax.pcolor(X, Y, H.max() - H.T, cmap=density_cmap)

            # Plot the contour edge colors.
            if plot_contours:
                if contour_kwargs is None:
                    contour_kwargs = dict()
                contour_kwargs["colors"] = contour_kwargs.get("colors", color)
                ax.contour(X2, Y2, H2.T, V, **contour_kwargs)

            ax.set_xlim(span[0])
            ax.set_ylim(span[1])

        def traceplot(results, span=None, quantiles=[0.16, 0.5, 0.84], smooth=0.02,
                      post_color='blue', post_kwargs=None, kde=True, nkde=1000,
                      trace_cmap='plasma', trace_color=None, trace_kwargs=None,
                      connect=False, connect_highlight=10, connect_color='red',
                      connect_kwargs=None, max_n_ticks=5, use_math_text=False,
                      labels=None, label_kwargs=None,
                      show_titles=False, title_fmt=".2f", title_kwargs=None,
                      truths=None, truth_color='red', truth_kwargs=None,
                      verbose=False, fig=None):
            """
            Plot traces and marginalized posteriors for each parameter.

            Parameters
            ----------
            results : :class:`~dynesty.results.Results` instance
                A :class:`~dynesty.results.Results` instance from a nested
                sampling run. **Compatible with results derived from**
                `nestle <http://kylebarbary.com/nestle/>`_.

            span : iterable with shape (ndim,), optional
                A list where each element is either a length-2 tuple containing
                lower and upper bounds or a float from `(0., 1.]` giving the
                fraction of (weighted) samples to include. If a fraction is provided,
                the bounds are chosen to be equal-tailed. An example would be::

                    span = [(0., 10.), 0.95, (5., 6.)]

                Default is `0.999999426697` (5-sigma credible interval) for each
                parameter.

            quantiles : iterable, optional
                A list of fractional quantiles to overplot on the 1-D marginalized
                posteriors as vertical dashed lines. Default is `[0.16, 0.5, 0.84]`
                (the 68%/1-sigma credible interval). Use `[0.025, 0.5, 0.975]`
                for 95%/2-sigma credible interval.

            smooth : float or iterable with shape (ndim,), optional
                The standard deviation (either a single value or a different value for
                each subplot) for the Gaussian kernel used to smooth the 1-D
                marginalized posteriors, expressed as a fraction of the span.
                Default is `0.02` (2% smoothing). If an integer is provided instead,
                this will instead default to a simple (weighted) histogram with
                `bins=smooth`.

            post_color : str or iterable with shape (ndim,), optional
                A `~matplotlib`-style color (either a single color or a different
                value for each subplot) used when plotting the histograms.
                Default is `'blue'`.

            post_kwargs : dict, optional
                Extra keyword arguments that will be used for plotting the
                marginalized 1-D posteriors.

            kde : bool, optional
                Whether to use kernel density estimation to estimate and plot
                the PDF of the importance weights as a function of log-volume
                (as opposed to the importance weights themselves). Default is
                `True`.

            nkde : int, optional
                The number of grid points used when plotting the kernel density
                estimate. Default is `1000`.

            trace_cmap : str or iterable with shape (ndim,), optional
                A `~matplotlib`-style colormap (either a single colormap or a
                different colormap for each subplot) used when plotting the traces,
                where each point is colored according to its weight. Default is
                `'plasma'`.

            trace_color : str or iterable with shape (ndim,), optional
                A `~matplotlib`-style color (either a single color or a
                different color for each subplot) used when plotting the traces.
                This overrides the `trace_cmap` option by giving all points
                the same color. Default is `None` (not used).

            trace_kwargs : dict, optional
                Extra keyword arguments that will be used for plotting the traces.

            connect : bool, optional
                Whether to draw lines connecting the paths of unique particles.
                Default is `False`.

            connect_highlight : int or iterable, optional
                If `connect=True`, highlights the paths of a specific set of
                particles. If an integer is passed, :data:`connect_highlight`
                random particle paths will be highlighted. If an iterable is passed,
                then the particle paths corresponding to the provided indices
                will be highlighted.

            connect_color : str, optional
                The color of the highlighted particle paths. Default is `'red'`.

            connect_kwargs : dict, optional
                Extra keyword arguments used for plotting particle paths.

            max_n_ticks : int, optional
                Maximum number of ticks allowed. Default is `5`.

            use_math_text : bool, optional
                Whether the axis tick labels for very large/small exponents should be
                displayed as powers of 10 rather than using `e`. Default is `False`.

            labels : iterable with shape (ndim,), optional
                A list of names for each parameter. If not provided, the default name
                used when plotting will follow :math:`x_i` style.

            label_kwargs : dict, optional
                Extra keyword arguments that will be sent to the
                `~matplotlib.axes.Axes.set_xlabel` and
                `~matplotlib.axes.Axes.set_ylabel` methods.

            show_titles : bool, optional
                Whether to display a title above each 1-D marginalized posterior
                showing the 0.5 quantile along with the upper/lower bounds associated
                with the 0.025 and 0.975 (95%/2-sigma credible interval) quantiles.
                Default is `True`.

            title_fmt : str, optional
                The format string for the quantiles provided in the title. Default is
                `'.2f'`.

            title_kwargs : dict, optional
                Extra keyword arguments that will be sent to the
                `~matplotlib.axes.Axes.set_title` command.

            truths : iterable with shape (ndim,), optional
                A list of reference values that will be overplotted on the traces and
                marginalized 1-D posteriors as solid horizontal/vertical lines.
                Individual values can be exempt using `None`. Default is `None`.

            truth_color : str or iterable with shape (ndim,), optional
                A `~matplotlib`-style color (either a single color or a different
                value for each subplot) used when plotting `truths`.
                Default is `'red'`.

            truth_kwargs : dict, optional
                Extra keyword arguments that will be used for plotting the vertical
                and horizontal lines with `truths`.

            verbose : bool, optional
                Whether to print the values of the computed quantiles associated with
                each parameter. Default is `False`.

            fig : (`~matplotlib.figure.Figure`, `~matplotlib.axes.Axes`), optional
                If provided, overplot the traces and marginalized 1-D posteriors
                onto the provided figure. Otherwise, by default an
                internal figure is generated.

            Returns
            -------
            traceplot : (`~matplotlib.figure.Figure`, `~matplotlib.axes.Axes`)
                Output trace plot.

            """

            # Initialize values.
            if title_kwargs is None:
                title_kwargs = dict()
            if label_kwargs is None:
                label_kwargs = dict()
            if trace_kwargs is None:
                trace_kwargs = dict()
            if connect_kwargs is None:
                connect_kwargs = dict()
            if post_kwargs is None:
                post_kwargs = dict()
            if truth_kwargs is None:
                truth_kwargs = dict()

            # Set defaults.
            connect_kwargs['alpha'] = connect_kwargs.get('alpha', 0.7)
            post_kwargs['alpha'] = post_kwargs.get('alpha', 0.6)
            trace_kwargs['s'] = trace_kwargs.get('s', 3)
            truth_kwargs['linestyle'] = truth_kwargs.get('linestyle', 'solid')
            truth_kwargs['linewidth'] = truth_kwargs.get('linewidth', 2)

            # Extract weighted samples.
            samples = results['samples']
            logvol = results['logvol']
            try:
                weights = np.exp(results['logwt'] - results['logz'][-1])
            except:
                weights = results['weights']
            if kde:
                # Derive kernel density estimate.
                wt_kde = gaussian_kde(_resample_equal(-logvol, weights))  # KDE
                logvol_grid = np.linspace(logvol[0], logvol[-1], nkde)  # resample
                wt_grid = wt_kde.pdf(-logvol_grid)  # evaluate KDE PDF
                wts = np.interp(-logvol, -logvol_grid, wt_grid)  # interpolate
            else:
                wts = weights

            # Deal with 1D results. A number of extra catches are also here
            # in case users are trying to plot other results besides the `Results`
            # instance generated by `dynesty`.
            samples = np.atleast_1d(samples)
            if len(samples.shape) == 1:
                samples = np.atleast_2d(samples)
            else:
                assert len(samples.shape) == 2, "Samples must be 1- or 2-D."
                samples = samples.T
            assert samples.shape[0] <= samples.shape[1], "There are more " \
                                                         "dimensions than samples!"
            ndim, nsamps = samples.shape

            # Check weights.
            if weights.ndim != 1:
                raise ValueError("Weights must be 1-D.")
            if nsamps != weights.shape[0]:
                raise ValueError("The number of weights and samples disagree!")

            # Check ln(volume).
            if logvol.ndim != 1:
                raise ValueError("Ln(volume)'s must be 1-D.")
            if nsamps != logvol.shape[0]:
                raise ValueError("The number of ln(volume)'s and samples disagree!")

            # Check sample IDs.
            if connect:
                try:
                    samples_id = results['samples_id']
                    uid = np.unique(samples_id)
                except:
                    raise ValueError("Sample IDs are not defined!")
                try:
                    ids = connect_highlight[0]
                    ids = connect_highlight
                except:
                    ids = np.random.choice(uid, size=connect_highlight, replace=False)

            # Determine plotting bounds for marginalized 1-D posteriors.
            if span is None:
                span = [0.999999426697 for i in range(ndim)]
            span = list(span)
            if len(span) != ndim:
                raise ValueError("Dimension mismatch between samples and span.")
            for i, _ in enumerate(span):
                try:
                    xmin, xmax = span[i]
                except:
                    q = [0.5 - 0.5 * span[i], 0.5 + 0.5 * span[i]]
                    span[i] = _quantile(samples[i], q, weights=weights)

            # Setting up labels.
            if labels is None:
                labels = [r"$x_{" + str(i + 1) + "}$" for i in range(ndim)]

            # Setting up smoothing.
            if (isinstance(smooth, int_type) or isinstance(smooth, float_type)):
                smooth = [smooth for i in range(ndim)]

            # Setting up default plot layout.
            if fig is None:
                fig, axes = plt.subplots(ndim, 2, figsize=(12, 3 * ndim))
            else:
                fig, axes = fig
                try:
                    axes.reshape(ndim, 2)
                except:
                    raise ValueError("Provided axes do not match the required shape "
                                     "for plotting samples.")

            # Plotting.
            for i, x in enumerate(samples):

                # Plot trace.

                # Establish axes.
                if np.shape(samples)[0] == 1:
                    ax = axes[1]
                else:
                    ax = axes[i, 0]
                # Set color(s)/colormap(s).
                if trace_color is not None:
                    if isinstance(trace_color, str_type):
                        color = trace_color
                    else:
                        color = trace_color[i]
                else:
                    color = wts
                if isinstance(trace_cmap, str_type):
                    cmap = trace_cmap
                else:
                    cmap = trace_cmap[i]
                # Setup axes.
                ax.set_xlim([0., -min(logvol)])
                ax.set_ylim([min(x), max(x)])
                if max_n_ticks == 0:
                    ax.xaxis.set_major_locator(NullLocator())
                    ax.yaxis.set_major_locator(NullLocator())
                else:
                    ax.xaxis.set_major_locator(MaxNLocator(max_n_ticks))
                    ax.yaxis.set_major_locator(MaxNLocator(max_n_ticks))
                # Label axes.
                sf = ScalarFormatter(useMathText=use_math_text)
                ax.yaxis.set_major_formatter(sf)
                ax.set_xlabel(r"$-\ln X$", **label_kwargs)
                ax.set_ylabel(labels[i], **label_kwargs)
                # Generate scatter plot.
                ax.scatter(-logvol, x, c=color, cmap=cmap, **trace_kwargs)
                if connect:
                    # Add lines highlighting specific particle paths.
                    for j in ids:
                        sel = (samples_id == j)
                        ax.plot(-logvol[sel], x[sel], color=connect_color,
                                **connect_kwargs)
                # Add truth value(s).
                if truths is not None and truths[i] is not None:
                    try:
                        [ax.axhline(t, color=truth_color, **truth_kwargs)
                         for t in truths[i]]
                    except:
                        ax.axhline(truths[i], color=truth_color, **truth_kwargs)

                # Plot marginalized 1-D posterior.

                # Establish axes.
                if np.shape(samples)[0] == 1:
                    ax = axes[0]
                else:
                    ax = axes[i, 1]
                # Set color(s).
                if isinstance(post_color, str_type):
                    color = post_color
                else:
                    color = post_color[i]
                # Setup axes
                ax.set_xlim(span[i])
                if max_n_ticks == 0:
                    ax.xaxis.set_major_locator(NullLocator())
                    ax.yaxis.set_major_locator(NullLocator())
                else:
                    ax.xaxis.set_major_locator(MaxNLocator(max_n_ticks))
                    ax.yaxis.set_major_locator(NullLocator())
                # Label axes.
                sf = ScalarFormatter(useMathText=use_math_text)
                ax.xaxis.set_major_formatter(sf)
                ax.set_xlabel(labels[i], **label_kwargs)
                # Generate distribution.
                s = smooth[i]
                if isinstance(s, int_type):
                    # If `s` is an integer, plot a weighted histogram with
                    # `s` bins within the provided bounds.
                    n, b, _ = ax.hist(x, bins=s, weights=weights, color=color,
                                      range=np.sort(span[i]), **post_kwargs)
                    x0 = np.array(list(zip(b[:-1], b[1:]))).flatten()
                    y0 = np.array(list(zip(n, n))).flatten()
                else:
                    # If `s` is a float, oversample the data relative to the
                    # smoothing filter by a factor of 10, then use a Gaussian
                    # filter to smooth the results.
                    bins = int(round(10. / s))
                    n, b = np.histogram(x, bins=bins, weights=weights,
                                        range=np.sort(span[i]))
                    n = norm_kde(n, 10.)
                    x0 = 0.5 * (b[1:] + b[:-1])
                    y0 = n
                    ax.fill_between(x0, y0, color=color, **post_kwargs)
                ax.set_ylim([0., max(y0) * 1.05])
                # Plot quantiles.
                if quantiles is not None and len(quantiles) > 0:
                    qs = _quantile(x, quantiles, weights=weights)
                    for q in qs:
                        ax.axvline(q, lw=1, ls="dashed", color=color)
                    if verbose:
                        print("Quantiles:")
                        print(labels[i], [blob for blob in zip(quantiles, qs)])
                # Add truth value(s).
                if truths is not None and truths[i] is not None:
                    try:
                        [ax.axvline(t, color=truth_color, **truth_kwargs)
                         for t in truths[i]]
                    except:
                        ax.axvline(truths[i], color=truth_color, **truth_kwargs)
                # Set titles.
                if show_titles:
                    title = None
                    if title_fmt is not None:
                        ql, qm, qh = _quantile(x, [0.16, 0.5, 0.84], weights=weights)
                        q_minus, q_plus = qm - ql, qh - qm
                        fmt = "{{0:{0}}}".format(title_fmt).format
                        title = r"${{{0}}}_{{-{1}}}^{{+{2}}}$"
                        title = title.format(fmt(qm), fmt(q_minus), fmt(q_plus))
                        title = "{0} = {1}".format(labels[i], title)
                        ax.set_title(title, **title_kwargs)

            return fig, axes

        def cornerplot(results, span=None, quantiles=[0.16, 0.5, 0.84],
                       color='black', smooth=0.02, hist_kwargs=None,
                       hist2d_kwargs=None, labels=None, label_kwargs=None,
                       show_titles=False, title_fmt=".2f", title_kwargs=None,
                       truths=None, truth_color='red', truth_kwargs=None,
                       max_n_ticks=5, top_ticks=False, use_math_text=False,
                       verbose=False, fig=None):
            # type: (object, object, object, object, object, object, object, object, object, object, object, object, object, object, object, object, object, object, object, object) -> object
            """
            Generate a corner plot of the 1-D and 2-D marginalized posteriors.

            Parameters
            ----------
            results : :class:`~dynesty.results.Results` instance
                A :class:`~dynesty.results.Results` instance from a nested
                sampling run. **Compatible with results derived from**
                `nestle <http://kylebarbary.com/nestle/>`_.

            span : iterable with shape (ndim,), optional
                A list where each element is either a length-2 tuple containing
                lower and upper bounds or a float from `(0., 1.]` giving the
                fraction of (weighted) samples to include. If a fraction is provided,
                the bounds are chosen to be equal-tailed. An example would be::

                    span = [(0., 10.), 0.95, (5., 6.)]

                Default is `0.999999426697` (5-sigma credible interval).

            quantiles : iterable, optional
                A list of fractional quantiles to overplot on the 1-D marginalized
                posteriors as vertical dashed lines. Default is `[0.16, 0.5, 0.84]`
                (spanning the 68%/2-sigma credible interval). Use `[0.025, 0.5, 0.975]`
                for 95%/2-sigma credible interval.

            color : str or iterable with shape (ndim,), optional
                A `~matplotlib`-style color (either a single color or a different
                value for each subplot) used when plotting the histograms.
                Default is `'black'`.

            smooth : float or iterable with shape (ndim,), optional
                The standard deviation (either a single value or a different value for
                each subplot) for the Gaussian kernel used to smooth the 1-D and 2-D
                marginalized posteriors, expressed as a fraction of the span.
                Default is `0.02` (2% smoothing). If an integer is provided instead,
                this will instead default to a simple (weighted) histogram with
                `bins=smooth`.

            hist_kwargs : dict, optional
                Extra keyword arguments to send to the 1-D (smoothed) histograms.

            hist2d_kwargs : dict, optional
                Extra keyword arguments to send to the 2-D (smoothed) histograms.

            labels : iterable with shape (ndim,), optional
                A list of names for each parameter. If not provided, the default name
                used when plotting will follow :math:`x_i` style.

            label_kwargs : dict, optional
                Extra keyword arguments that will be sent to the
                `~matplotlib.axes.Axes.set_xlabel` and
                `~matplotlib.axes.Axes.set_ylabel` methods.

            show_titles : bool, optional
                Whether to display a title above each 1-D marginalized posterior
                showing the 0.5 quantile along with the upper/lower bounds associated
                with the 0.025 and 0.975 (95%/2-sigma credible interval) quantiles.
                Default is `True`.

            title_fmt : str, optional
                The format string for the quantiles provided in the title. Default is
                `'.2f'`.

            title_kwargs : dict, optional
                Extra keyword arguments that will be sent to the
                `~matplotlib.axes.Axes.set_title` command.

            truths : iterable with shape (ndim,), optional
                A list of reference values that will be overplotted on the traces and
                marginalized 1-D posteriors as solid horizontal/vertical lines.
                Individual values can be exempt using `None`. Default is `None`.

            truth_color : str or iterable with shape (ndim,), optional
                A `~matplotlib`-style color (either a single color or a different
                value for each subplot) used when plotting `truths`.
                Default is `'red'`.

            truth_kwargs : dict, optional
                Extra keyword arguments that will be used for plotting the vertical
                and horizontal lines with `truths`.

            max_n_ticks : int, optional
                Maximum number of ticks allowed. Default is `5`.

            top_ticks : bool, optional
                Whether to label the top (rather than bottom) ticks. Default is
                `False`.

            use_math_text : bool, optional
                Whether the axis tick labels for very large/small exponents should be
                displayed as powers of 10 rather than using `e`. Default is `False`.

            verbose : bool, optional
                Whether to print the values of the computed quantiles associated with
                each parameter. Default is `False`.

            fig : (`~matplotlib.figure.Figure`, `~matplotlib.axes.Axes`), optional
                If provided, overplot the traces and marginalized 1-D posteriors
                onto the provided figure. Otherwise, by default an
                internal figure is generated.

            Returns
            -------
            cornerplot : (`~matplotlib.figure.Figure`, `~matplotlib.axes.Axes`)
                Output corner plot.

            """

            # Initialize values.
            if quantiles is None:
                quantiles = []
            if truth_kwargs is None:
                truth_kwargs = dict()
            if label_kwargs is None:
                label_kwargs = dict()
            if title_kwargs is None:
                title_kwargs = dict()
            if hist_kwargs is None:
                hist_kwargs = dict()
            if hist2d_kwargs is None:
                hist2d_kwargs = dict()

            # Set defaults.
            hist_kwargs['alpha'] = hist_kwargs.get('alpha', 0.6)
            hist2d_kwargs['alpha'] = hist2d_kwargs.get('alpha', 0.6)
            truth_kwargs['linestyle'] = truth_kwargs.get('linestyle', 'solid')
            truth_kwargs['linewidth'] = truth_kwargs.get('linewidth', 2)
            truth_kwargs['alpha'] = truth_kwargs.get('alpha', 0.7)

            # Extract weighted samples.
            samples = results['samples']
            try:
                weights = np.exp(results['logwt'] - results['logz'][-1])
            except:
                weights = results['weights']

            # Deal with 1D results. A number of extra catches are also here
            # in case users are trying to plot other results besides the `Results`
            # instance generated by `dynesty`.
            samples = np.atleast_1d(samples)
            if len(samples.shape) == 1:
                samples = np.atleast_2d(samples)
            else:
                assert len(samples.shape) == 2, "Samples must be 1- or 2-D."
                samples = samples.T
            assert samples.shape[0] <= samples.shape[1], "There are more " \
                                                         "dimensions than samples!"
            ndim, nsamps = samples.shape

            # Check weights.
            if weights.ndim != 1:
                raise ValueError("Weights must be 1-D.")
            if nsamps != weights.shape[0]:
                raise ValueError("The number of weights and samples disagree!")

            # Set labels
            if labels is None:
                labels = [r"$x_{" + str(i + 1) + "}$" for i in range(ndim)]

            # Setting up smoothing.
            if (isinstance(smooth, int_type) or isinstance(smooth, float_type)):
                smooth = [smooth for i in range(ndim)]

            # Setup axis layout (from `corner.py`).
            factor = 2.0  # size of side of one panel
            lbdim = 0.5 * factor  # size of left/bottom margin
            trdim = 0.2 * factor  # size of top/right margin
            whspace = 0.05  # size of width/height margin
            plotdim = factor * ndim + factor * (ndim - 1.) * whspace  # plot size
            dim = lbdim + plotdim + trdim  # total size

            # Initialize figure.
            if fig is None:
                fig, axes = plt.subplots(ndim, ndim, figsize=(dim, dim))
            else:
                try:
                    fig, axes = fig
                    if np.shape(samples)[0] > 1:
                        axes = np.array(axes).reshape((ndim, ndim))
                except:
                    raise ValueError("Mismatch between axes and dimension.")

            # Format figure.
            lb = lbdim / dim
            tr = (lbdim + plotdim) / dim
            fig.subplots_adjust(left=lb, bottom=lb, right=tr, top=tr,
                                wspace=whspace, hspace=whspace)

            # Plotting.
            for i, x in enumerate(samples):
                if np.shape(samples)[0] == 1:
                    ax = axes
                else:
                    ax = axes[i, i]

                # Plot the 1-D marginalized posteriors.

                # Setup axes
                ax.set_xlim(span[i])
                if max_n_ticks == 0:
                    ax.xaxis.set_major_locator(NullLocator())
                    ax.yaxis.set_major_locator(NullLocator())
                else:
                    ax.xaxis.set_major_locator(MaxNLocator(max_n_ticks,
                                                           prune="lower"))
                    ax.yaxis.set_major_locator(NullLocator())
                # Label axes.
                sf = ScalarFormatter(useMathText=use_math_text)
                ax.xaxis.set_major_formatter(sf)
                if i < ndim - 1:
                    if top_ticks:
                        ax.xaxis.set_ticks_position("top")
                        [l.set_rotation(45) for l in ax.get_xticklabels()]
                    else:
                        ax.set_xticklabels([])
                else:
                    [l.set_rotation(45) for l in ax.get_xticklabels()]
                    ax.set_xlabel(labels[i], **label_kwargs)
                    ax.xaxis.set_label_coords(0.5, -0.3)
                # Generate distribution.
                sx = smooth[i]

                n, b, _ = ax.hist(x, bins=75, weights=weights, color=color, range=np.sort(span[i]), **hist_kwargs)

                # if isinstance(sx, int_type):
                #     # If `sx` is an integer, plot a weighted histogram with
                #     # `sx` bins within the provided bounds.
                #     n, b, _ = ax.hist(x, bins=sx, weights=weights, color=color,
                #                       range=np.sort(span[i]), **hist_kwargs)
                # else:
                #     # If `sx` is a float, oversample the data relative to the
                #     # smoothing filter by a factor of 10, then use a Gaussian
                #     # filter to smooth the results.
                #     bins = int(round(10. / sx))
                #     n, b = np.histogram(x, bins=bins, weights=weights,
                #                         range=np.sort(span[i]))
                #     n = norm_kde(n, 10.)
                #     b0 = 0.5 * (b[1:] + b[:-1])
                #     n, b, _ = ax.hist(b0, bins=b, weights=n,
                #                       range=np.sort(span[i]), color=color,
                #                       **hist_kwargs)

                if ax.get_ylim()[1] < max(n) * 1.05:
                    ax.set_ylim([0., max(n) * 1.05])
                else:
                    pass
                # Plot quantiles.
                if quantiles is not None and len(quantiles) > 0:
                    qs = _quantile(x, quantiles, weights=weights)
                    # for q in qs:
                    #     ax.axvline(q, lw=1, ls="dashed", color=color)
                    if verbose:
                        print("Quantiles:")
                        print(labels[i], [blob for blob in zip(quantiles, qs)])
                # Add truth value(s).
                if truths is not None and truths[i] is not None:
                    try:
                        [ax.axvline(t, color=truth_color, **truth_kwargs)
                         for t in truths[i]]
                    except:
                        ax.axvline(truths[i], color=truth_color, **truth_kwargs)
                # Set titles.
                if show_titles:
                    title = None
                    if title_fmt is not None:
                        ql, qm, qh = _quantile(x, [0.16, 0.5, 0.84], weights=weights)
                        q_minus, q_plus = qm - ql, qh - qm
                        fmt = "{{0:{0}}}".format(title_fmt).format
                        title = r"${{{0}}}_{{-{1}}}^{{+{2}}}$"
                        title = title.format(fmt(qm), fmt(q_minus), fmt(q_plus))
                        title = "{0} = {1}".format(labels[i], title)
                        ax.set_title(title, **title_kwargs)

                for j, y in enumerate(samples):
                    if np.shape(samples)[0] == 1:
                        ax = axes
                    else:
                        ax = axes[i, j]

                    # Plot the 2-D marginalized posteriors.

                    # Setup axes.
                    if j > i:
                        ax.set_frame_on(False)
                        ax.set_xticks([])
                        ax.set_yticks([])
                        continue
                    elif j == i:
                        continue

                    if max_n_ticks == 0:
                        ax.xaxis.set_major_locator(NullLocator())
                        ax.yaxis.set_major_locator(NullLocator())
                    else:
                        ax.xaxis.set_major_locator(MaxNLocator(max_n_ticks,
                                                               prune="lower"))
                        ax.yaxis.set_major_locator(MaxNLocator(max_n_ticks,
                                                               prune="lower"))
                    # Label axes.
                    sf = ScalarFormatter(useMathText=use_math_text)
                    ax.xaxis.set_major_formatter(sf)
                    ax.yaxis.set_major_formatter(sf)
                    if i < ndim - 1:
                        ax.set_xticklabels([])
                    else:
                        [l.set_rotation(45) for l in ax.get_xticklabels()]
                        ax.set_xlabel(labels[j], **label_kwargs)
                        ax.xaxis.set_label_coords(0.5, -0.3)
                    if j > 0:
                        ax.set_yticklabels([])
                    else:
                        [l.set_rotation(45) for l in ax.get_yticklabels()]
                        ax.set_ylabel(labels[i], **label_kwargs)
                        ax.yaxis.set_label_coords(-0.3, 0.5)
                    # Generate distribution.
                    sy = smooth[j]
                    check_ix = isinstance(sx, int_type)
                    check_iy = isinstance(sy, int_type)
                    if check_ix and check_iy:
                        fill_contours = False
                        plot_contours = False
                    else:
                        fill_contours = True
                        plot_contours = True
                    hist2d_kwargs['fill_contours'] = hist2d_kwargs.get('fill_contours',
                                                                       fill_contours)
                    hist2d_kwargs['plot_contours'] = hist2d_kwargs.get('plot_contours',
                                                                       plot_contours)
                    _hist2d(y, x, ax=ax, span=[span[j], span[i]],
                            weights=weights, color=color, smooth=[sy, sx],
                            **hist2d_kwargs)

                    # Add truth values
                    if truths is not None:
                        if truths[j] is not None:
                            try:
                                [ax.axvline(t, color=truth_color, **truth_kwargs)
                                 for t in truths[j]]
                            except:
                                ax.axvline(truths[j], color=truth_color,
                                           **truth_kwargs)
                        if truths[i] is not None:
                            try:
                                [ax.axhline(t, color=truth_color, **truth_kwargs)
                                 for t in truths[i]]
                            except:
                                ax.axhline(truths[i], color=truth_color,
                                           **truth_kwargs)

            return fig, axes

        print('Generating the Posterior Distribution Functions (PDFs) plot')
        _corner_parameters()
        if mds < 2:
            nest_out = _store_nest_solutions()
            gas_pos = _posteriors_clr_to_vmr(prefix, modes=None)

            data = np.loadtxt(prefix + '.txt')
            i = data[:, 1].argsort()[::-1]
            samples = data[i, 2:]
            weights = data[i, 0]
            loglike = data[i, 1]
            Z = nest_out['global_logE'][0]
            logvol = log(weights) + 0.5 * loglike + Z
            logvol = logvol - logvol.max()

            result = dict(samples=samples, weights=weights, logvol=logvol)

            parameters = json.load(open(prefix + 'params.json'))

            traceplot(result, labels=parameters, show_titles=True)
            plt.savefig(prefix + 'Nest_trace.pdf', bbox_inches='tight')
            plt.close()

            bound = _plotting_bounds(result, gas_pos)

            if self.param['truths'] is not None:
                tru = np.loadtxt(self.param['truths'])
                cornerplot(result, labels=parameters, show_titles=True, truths=list(tru), span=bound)
            else:
                cornerplot(result, labels=parameters, show_titles=True, span=bound)
            plt.savefig(prefix + 'Nest_posteriors.pdf', bbox_inches='tight')
            plt.close()

            os.system('mv ' + prefix + 'params.json ' + prefix + '_PostProcess.json')
            os.system('mv ' + prefix + 'params_original.json ' + prefix + 'params.json')
            os.system('mv ' + prefix + '.txt ' + prefix + '_PostProcess.txt')
            os.system('mv ' + prefix + 'original.txt ' + prefix + '.txt')

            with open(self.param['out_dir'] + 'stats_summary.txt', 'w') as fl:
                fl.write('############### SUMMARY STATISTICS ###############\n')
                fl.write('\n')
                if self.param['plot_models']:
                    fl.write('chi-square (d.o.f) = ' + str(round(self.param['chi_square_stat']['chi2'], 2)) + ' (' + str(self.param['chi_square_stat']['dof']) + ')\n')
                    fl.write('Reduced chi-square = ' + str(round(self.param['chi_square_stat']['chi2_red'], 2)) + '\n')
                fl.write('ln Z               = ' + str(round(Z, 2)) + ' +- ' + str(round(nest_out['global_logE'][1], 2)) + '\n')
                fl.write('\n')
                fl.write('##################################################\n')

        else:
            colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
            nest_out = _store_nest_solutions()
            result = {}

            for modes in range(0, mds):
                gas_pos = _posteriors_clr_to_vmr(prefix, modes=modes)

                data = np.loadtxt(prefix + 'solution' + str(modes) + '.txt')
                i = data[:, 1].argsort()[::-1]
                samples = data[i, 2:]
                weights = data[i, 0]
                loglike = data[i, 1]
                Z = nest_out['solutions']['solution' + str(modes)]['local_logE'][0]
                logvol = log(weights) + 0.5 * loglike + Z
                logvol = logvol - logvol.max()

                result[str(modes)] = dict(samples=samples, weights=weights, logvol=logvol)

                with open(self.param['out_dir'] + 'stats_summary.txt', 'a') as fl:
                    fl.write('*** SOLUTION ' + str(modes + 1) + ' ***\n')
                    fl.write('############### SUMMARY STATISTICS ###############\n')
                    fl.write('\n')
                    if self.param['plot_models']:
                        fl.write('chi-square (d.o.f) = ' + str(round(self.param['chi_square_stat']['solution_' + str(modes + 1)]['chi2'], 2)) + ' (' + str(self.param['chi_square_stat']['solution_' + str(modes + 1)]['dof']) + ')\n')
                        fl.write('Reduced chi-square = ' + str(round(self.param['chi_square_stat']['solution_' + str(modes + 1)]['chi2_red'], 2)) + '\n')
                    fl.write('ln Z               = ' + str(round(nest_out['solutions']['solution' + str(modes)]['local_logE'][0], 2)) + ' +- ' + str(round(nest_out['solutions']['solution' + str(modes)]['local_logE'][1], 2)) + '\n')
                    fl.write('\n')
                    fl.write('##################################################\n\n')

            parameters = json.load(open(prefix + 'params.json'))

            for modes in range(0, mds):
                traceplot(result[str(modes)], labels=parameters, show_titles=True)
                plt.savefig(prefix + 'Nest_trace (solution' + str(modes + 1) + ').pdf', bbox_inches='tight')
                plt.close()

            for modes in range(0, mds):
                bound = _plotting_bounds(result[str(modes)], modes=None)
                if self.param['truths'] is not None:
                    tru = np.loadtxt(self.param['truths'])
                    cornerplot(result[str(modes)], labels=parameters, show_titles=True, truths=list(tru), color=colors[modes], span=bound)
                else:
                    cornerplot(result[str(modes)], labels=parameters, show_titles=True, color=colors[modes], span=bound)

                plt.savefig(prefix + 'Nest_posteriors (solution' + str(modes + 1) + ').pdf', bbox_inches='tight')
                plt.close()

            try:
                if figu is not None:
                    pass
            except NameError:
                figu = None

            bound = _plotting_bounds(result, gas_pos, modes=mds)

            for modes in range(0, mds):
                if self.param['truths'] is not None:
                    tru = np.loadtxt(self.param['truths'])
                    figu = cornerplot(result[str(modes)], labels=parameters, show_titles=False, truths=list(tru), color=colors[modes], fig=figu, span=bound)
                else:
                    figu = cornerplot(result[str(modes)], labels=parameters, show_titles=False, color=colors[modes], fig=figu, span=bound)

                os.system('mv ' + prefix + 'solution' + str(modes) + '.txt ' + prefix + 'solution' + str(modes + 1) + '_PostProcess.txt')

            for modes in range(mds-1, -1, -1):
                os.system('mv ' + prefix + 'solution' + str(modes) + '_original.txt ' + prefix + 'solution' + str(modes + 1) + '.txt')

            plt.savefig(prefix + 'Nest_posteriors.pdf', bbox_inches='tight')
            plt.close()

            os.system('mv ' + prefix + 'params.json ' + prefix + '_PostProcess.json')
            os.system('mv ' + prefix + 'params_original.json ' + prefix + 'params.json')

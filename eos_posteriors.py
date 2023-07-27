from pathlib import Path
from tqdm.auto import tqdm
import os

import argparse

import numpy as np
from astropy.table import Table

import scipy.interpolate
import bilby


Neos = 5000
EOS_raw = np.arange(1, Neos+1)  # the gwem_resampling will add one to this


def lambdas2lambdaTs(lambda1, lambda2, q):
    eta = q / np.power(1. + q, 2.)
    eta2 = eta * eta
    eta3 = eta2 * eta
    root14eta = np.sqrt(1. - 4 * eta)

    lambdaT = (8. / 13.) * ((1. + 7 * eta - 31 * eta2) * (lambda1 + lambda2) + root14eta * (1. + 9 * eta - 11. * eta2) * (lambda1 - lambda2))
    dlambdaT = 0.5 * (root14eta * (1. - 13272. * eta / 1319. + 8944. * eta2 / 1319.) * (lambda1 + lambda2) + (1. - 15910. * eta / 1319. + 32850. * eta2 / 1319. + 3380. * eta3 / 1319.) * (lambda1 - lambda2))

    return lambdaT, dlambdaT


def get_parser():

    parser = argparse.ArgumentParser(
        description="Determine EOS_posteriors from Lamda_tilde."
    )
    parser.add_argument(
        "--outdir", 
        type=str, 
        default="outdir",
        help="Path to the output directory"
    )
    parser.add_argument(
        "--GWsamples",
        type=str,
        required=True,
        help="Data from the posterior of the GW injection",
    )
    parser.add_argument(
        "--EOSpath", 
        type=str, 
        required=True, 
        help="The EOS data"
    )
    return parser

def main(args=None):

    if args is None:
        parser = get_parser()
        args = parser.parse_args()
    
    outdir = args.outdir            
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    ### Read files with tabulated eos, and create interpolant for mass-lambda for each eos, saving it into a dictionary
    EOS_lambda_interp_dict = {}
    
    for i in EOS_raw:
        r, m, Lambda = np.loadtxt('{0}/{1}.dat'.format(args.EOSpath, i), usecols=[0, 1, 2], unpack=True)

        interp_lambda_of_m = scipy.interpolate.interp1d(m, Lambda, kind='linear', fill_value="extrapolate")
        EOS_lambda_interp_dict[i] = interp_lambda_of_m

    gwsamples = Path(args.GWsamples)
    print(gwsamples)
    
    simulation_id = args.GWsamples.split("/")[-1].split("_")[-2]
    run_num =       args.GWsamples.split("/")[-1].split("_")[0][3:] 
    print(simulation_id)
    print(run_num)
    
    
    ### Read gw samples
    result = bilby.result.read_in_result(gwsamples)
    pos = result.posterior
    print(pos.keys())

    ## posterior has 'chirp_mass', 'mass_ratio', 'chi_1', 'chi_2', 'luminosity_distance',
    ##           'dec', 'ra', 'cos_theta_jn', 'psi', 'phase', 'lambda_tilde'


    #####################################
    ### Method 2, compute for each sample
    #####################################
    lambdaT_obs = pos['lambda_tilde']

    mchirp = pos['chirp_mass']
    mratio = pos['mass_ratio']
    ldistance = pos['luminosity_distance']

    print("==============================================\n")
    print("Conversion luminosity distance to redshift")

    z =  bilby.gw.conversion.luminosity_distance_to_redshift(ldistance)
    zp1 = z + 1 
    print("==============================================\n")

    total_mass = bilby.gw.conversion.chirp_mass_and_mass_ratio_to_total_mass(mchirp, mratio)
    m1, m2 = bilby.gw.conversion.total_mass_and_mass_ratio_to_component_masses(mratio, total_mass)


    # convert detector masses to the sources masses 
    mass_1_source = m1/zp1
    mass_2_source = m2/zp1


    lambda_tilde_obs = []
    eos_id = []
    eos_samples_lik = []
    lambda_eos   = []
    eos_samples =[]
    mass_id = []
    with tqdm(total=len(mass_1_source)) as progress:

        for ii in range(len(mass_1_source)):

            lambdaT = []
            for eos in EOS_raw:
                try:
                    lambda1 = EOS_lambda_interp_dict[eos](mass_1_source[ii])
                    lambda2 = EOS_lambda_interp_dict[eos](mass_2_source[ii])

                    q = mass_1_source[ii]/mass_2_source[ii]
                    lambdaT_sample, dlamT_sample = lambdas2lambdaTs(lambda1, lambda2, q)

                    lambdaT.append(lambdaT_sample)
                    eos_id.append(eos)


                    del lambda1, lambda2, lambdaT_sample, dlamT_sample, q

                except ValueError as err:
                    print(err)

            if len(lambdaT)>0:

                dist = abs(np.array(lambdaT) - lambdaT_obs[ii])
                dmin_idx = np.where(dist==np.min(dist))[0][0]

                eos_samples.append(EOS_raw[dmin_idx])
                lambda_eos.append(lambdaT[dmin_idx])
                eos_samples_lik.append(np.exp(-(lambdaT[dmin_idx] - pos['lambda_tilde'][ii])))

                lambda_tilde_obs.append(pos['lambda_tilde'][ii])
                mass_id.append(ii)
            ###not sure likelihood is needed

                del dist, dmin_idx, lambdaT

            progress.update()


        tables_pos = Table({"mass_1_source":mass_1_source, "mass_2_source":mass_2_source, "chirp_mass":pos['chirp_mass'], "mass_ratio":pos['mass_ratio'], "luminosity_distance":pos['luminosity_distance'], "EOS":eos_samples, "redshift":z,  "chi_1":pos['chi_1'], "chi_2":pos['chi_2'],  "dec":pos['dec'], "ra":pos['ra'], "cos_theta_jn":pos['cos_theta_jn'], "psi":pos['psi'], "phase":pos['phase'], "lambda_tilde":pos['lambda_tilde']})
                        
        tables_pos.write(f"{outdir}/gw_simulation_{simulation_id}_posterior_samples_{run_num}.dat", format='ascii.tab', overwrite=True)


if __name__ == "__main__":
    main()

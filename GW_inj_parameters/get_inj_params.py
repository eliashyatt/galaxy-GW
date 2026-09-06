from bilby.core import prior
import pandas as pd
import scipy.stats as stats
import numpy as np
import argparse

def generate_inj_params(N):
    priors = prior.PriorDict().from_json('./BBH_15_parameters_O4_prior.json')

    injection_parameters = priors.sample(N)

    df = pd.DataFrame(injection_parameters)
    df['mass_2_source'] = df['mass_ratio'] * df['mass_1_source']

    GW_df_full = pd.read_csv("2700_GW_sources.csv")

    #gamma distribution parameters from GWTC-4
    gamma = stats.gamma
    GW_a, GW_loc, GW_scale = 2.4713875177782727, 0, 1193.2568157708797
    GW_probs = gamma.pdf(GW_df_full['luminosity_distance'], GW_a, GW_loc, GW_scale)
    GW_probs = GW_probs / np.sum(GW_probs)

    GW_df_indices = np.random.choice(len(GW_df_full), N, replace=False, p=GW_probs)
    GW_df = GW_df_full.iloc[GW_df_indices]

    df['luminosity_distance'] = np.array(GW_df['luminosity_distance'])
    df['dec'] = np.deg2rad(np.array(GW_df['Dec_deg']))
    df['ra'] = np.deg2rad(np.array(GW_df['RA_deg']))
    df['event'] = [str(f"{i:05d}") for i in range(1, N+1)]

    df.to_csv(f'inj_parameters_{N}.csv', index=False)
    #df.to_csv('/home/vboxuser/Desktop/Code_for_Thesis/galaxy-GW/GW_bayestar/v2/inj_parameters_all_15000.csv', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate injection parameters for GW events')
    parser.add_argument('--N', type=int, help='Number of injection parameters to generate. Suggested 3 times the number of GW events wanted', required=True)
    args = parser.parse_args()

    N = args.N

    print(f'Generating {N} injection parameters...')

    generate_inj_params(N)

    print('Done!')


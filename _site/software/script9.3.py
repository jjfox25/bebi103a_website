import itertools
import warnings

import numpy as np
import pandas as pd
import scipy
import scipy.stats as st
import tqdm

import iqplot

import bebi103

import bokeh.io
import bokeh.layouts
import colorcet
from bokeh.palettes import Category10_10 as palette
import bokeh.plotting

def main(data_path):
	df = pd.read_csv(data_path)
	bac1 = df.loc[df['bacterium'] == 1]
	bac2 = df.loc[df['bacterium'] == 2]
	bac_list = [bac1, bac2]
	bac_list_string = ['bac1','bac2']
	
	bac2['growth event'] = bac2['growth event'] + 20

	bac1_unique = list(bac1['growth event'].unique())
	bac2_unique = list(bac2['growth event'].unique())

	ge_unique = bac1_unique + bac2_unique

	df = pd.concat([bac1, bac2], axis = 0)

	master_list = []

	for i in ge_unique:
    		temp = df.loc[df['growth event']==i]
    		temp_list = []
    		for n in range(len(temp)):
        		temp_list.append(n+1)
   		master_list = master_list + temp_list

	df['growth_time'] = master_list

def theor_area_linear(a0, k, t):
    '''Linear model of bacterial growth.'''
    return a0*(1 + k*t)

def theor_area_exp(a0, k, t):
    '''Exponential model of bacterial growth.'''
    return a0*np.exp(k*t)

def gen_linear(params, t, size, rg):
    """Generate a new linear data set using the generative model."""
    return rg.normal(theor_area_linear(*params[:-1], t), params[-1])

def gen_exp(params, t, size, rg):
    """Generate a new exponential data set using the generative model."""
    return rg.normal(theor_area_exp(*params[:-1], t), params[-1])

def resid_linear(params, t, a):
    """Residual for linear model."""
    a0, k = params
    return a - theor_area_linear(*params, t)

def resid_exp(params, t, a):
    """Residual for exponential model."""
    a0, k = params
    return a - theor_area_exp(*params, t)

def mle(a, t):
    """Compute MLE of parameters for linear model."""
    # Perform least squares
    res = scipy.optimize.least_squares(
        resid_linear,
        np.array([1.4, 0.01]),
        args=(t, a),
        bounds=([0, 0], [np.inf, np.inf]),
    )

    # Compute residual sum of squares from optimal params
    rss_mle = np.sum(resid_linear(res.x, t, a) ** 2)

    # Compute MLE for sigma
    sigma_mle = np.sqrt(rss_mle / len(t))

    return np.concatenate((res.x, (sigma_mle,)))

def mle_exp(a, t):
    """Compute MLE of parameters for exponential model."""
    # Perform least squares
    res = scipy.optimize.least_squares(
        resid_exp,
        np.array([1.4, 0.01]),
        args=(t, a),
        bounds=([0, 0], [np.inf, np.inf]),
    )

    # Compute residual sum of squares from optimal params
    rss_mle = np.sum(resid_exp(res.x, t, a) ** 2)

    # Compute MLE for sigma
    sigma_mle = np.sqrt(rss_mle / len(t))

    return np.concatenate((res.x, (sigma_mle,)))

def log_likelihood_linear(params, t):
    """Log likelihood of linear model."""
    a0, k, sigma = params

    if a0 <= 0 or k <= 0:
        return -np.inf

    mu = theor_area_linear(a0, k, t)
    return np.sum(st.norm.logpdf(mu, sigma))

def log_likelihood_exp(params, t):
    """Log likelihood of linear model."""
    a0, k, sigma = params

    if k <= 0:
        return -np.inf

    mu = theor_area_exp(a0, k, t)
    return np.sum(st.norm.logpdf(mu, sigma))

df = pd.read_csv('https://s3.amazonaws.com/bebi103.caltech.edu/data/caulobacter_growth_events.csv')

bac1 = df.loc[df['bacterium'] == 1]
bac2 = df.loc[df['bacterium'] == 2]

bac_list = [bac1, bac2]
bac_list_string = ['bac1','bac2']

bac2['growth event'] = bac2['growth event'] + 20

bac1_unique = list(bac1['growth event'].unique())
bac2_unique = list(bac2['growth event'].unique())

ge_unique = bac1_unique + bac2_unique

df = pd.concat([bac1, bac2], axis = 0)

master_list = []

for i in ge_unique:
    temp = df.loc[df['growth event']==i]
    temp_list = []
    for n in range(len(temp)):
        temp_list.append(n+1)
    master_list = master_list + temp_list

df['growth_time'] = master_list

print('MLE linear')
results = np.empty((len(df["growth event"].unique()), 3))
for event, g in df.groupby("growth event"):
    results[event] = mle(
        g["area (µm²)"].values, g["growth_time"].values
    )
df_mle = pd.DataFrame(
    data=results,
    columns=["a0", "k", "σ (a.u.)"],
)
print('a0_avg_lin: ', df_mle['a0'].mean())
print('k_avg_lin: ', df_mle['k'].mean())

print('MLE exponential')
results_exp = np.empty((len(df["growth event"].unique()), 3))
for event, g in df.groupby("growth event"):
    results_exp[event] = mle_exp(
        g["area (µm²)"].values, g["growth_time"].values
    )
df_mle_exp = pd.DataFrame(
    data=results_exp,
    columns=["a0", "k", "σ (a.u.)"],
)
print('a0_avg_exp: ', df_mle_exp['a0'].mean())
print('k_avg_exp: ', df_mle_exp['k'].mean())

reps = {}
for event, g in tqdm.tqdm(df.groupby("growth event")):
    # Extract time points and area measurements
    t = g["growth_time"].values
    data = g["area (µm²)"].values

    # Generate bootstrap replicates
    reps[event] = bebi103.bootstrap.draw_bs_reps_mle(
        mle,
        gen_linear,
        data=data,
        mle_args=(t,),
        gen_args=(t,),
        size=1000,
        n_jobs=3,
    )

reps_exp = {}
for event, g in tqdm.tqdm(df.groupby("growth event")):
    # Extract time points and area measurements
    t = g["growth_time"].values
    data = g["area (µm²)"].values

    # Generate bootstrap replicates
    reps_exp[event] = bebi103.bootstrap.draw_bs_reps_mle(
        mle_exp,
        gen_exp,
        data=data,
        mle_args=(t,),
        gen_args=(t,),
        size=1000,
        n_jobs=3,
    )

AIC_lin = []
AIC_exp = []

for i in ge_unique:
    temp = df.loc[df['growth event']==i]
    
    # compute log likelihood
    temp_loglike_lin = log_likelihood_linear(df_mle.iloc[i,:], temp["growth_time"])
    temp_loglike_exp = log_likelihood_exp(df_mle_exp.iloc[i,:], temp["growth_time"])
    
    # compute AIC defined as -2 * (log_likelihood - number_of_parameters)
    temp_AIC_lin = -2 * (temp_loglike_lin - 3)
    temp_AIC_exp = -2 * (temp_loglike_exp - 3)
    
    AIC_lin.append(temp_AIC_lin)
    AIC_exp.append(temp_AIC_exp)

results = pd.DataFrame(np.array([ge_unique, 
                                 df_mle.iloc[:,0], 
                                 df_mle.iloc[:,1],
                                 df_mle.iloc[:,2],
                                 df_mle_exp.iloc[:,0], 
                                 df_mle_exp.iloc[:,1],
                                 df_mle_exp.iloc[:,2],
                                 AIC_lin, 
                                 AIC_exp
                                ]).T)
results.columns = ['growth event',
                   'a0_lin',
                   'k_lin',
                   'sigma_lin',
                   'a0_exp',
                   'k_exp',
                   'sigma_exp',
                   'AIC_lin',
                   'AIC_exp'
]

results['AIC_ratio_lin_to_exp'] = results['AIC_lin'] / results['AIC_exp']

results['w_linear'] = np.exp(-results['AIC_lin']/2) / (np.exp(-results['AIC_lin']/2) + np.exp(-results['AIC_exp']/2))
results['w_exp'] = np.exp(-results['AIC_exp']/2) / (np.exp(-results['AIC_exp']/2) + np.exp(-results['AIC_lin']/2))


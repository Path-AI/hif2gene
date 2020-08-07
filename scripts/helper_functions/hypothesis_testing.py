from builtins import list

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
from collections import Counter

from sklearn.decomposition import PCA
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score

from scipy.stats import mannwhitneyu
from scipy.stats import ks_2samp
from scipy.interpolate import UnivariateSpline
from scipy.stats import ranksums

from statsmodels.stats.proportion import proportion_confint
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import ConvergenceError
from lifelines.plotting import add_at_risk_counts
from lifelines.utils import k_fold_cross_validation
from lifelines.utils import concordance_index


def permutation_test(x,y,test=None):
    return 'Not implemented'


def fishers_exact_test(df,binary_1,binary_2):
    '''
    Compute Fisher's Exact test for two columns in a dataframe
    :param df: pandas dataframe
    :param binary_1: category column 1
    :param binary_2: category colum 2
    :return: odds_ratio and p value
    '''
    odds_ratio,p = fisher_exact(pd.crosstab(df[binary_1] > 0, df[binary_2] > 0).values)
    return odds_ratio,p


def compute_logrank_test(df1, df2, select_df1, select_df2,
                         surv_col='pfsinv', event_col='pfsinv_event_observed'):
    '''
    :param df1: pandas dataframe for group 1
    :param df2: pandas dataframe for group 1
    :param select_df1:
    :param select_df2:
    :param surv_col:
    :param event_col:
    :return:
    '''
    res = logrank_test(df1[surv_col].values[select_df1], df2[surv_col].values[select_df2],
                       df1[event_col].values[select_df1], df2[event_col].values[select_df2])
    p = np.round(res.p_value, 4)
    return p

def qvalue(pvals, alpha=0.05, verbose=True):
    """Function for estimating q-values from p-values using the Storey-
    Tibshirani q-value method (2003).
    :param pvals: numpy array of p-values
    :param alpha: desired FDR
    :return:
    significant - An array of flags indicating which p-values are significant.
    qvals - Q-values corresponding to the p-values.
    """

    """Count the p-values. Find indices for sorting the p-values into
    ascending order and for reversing the order back to original."""
    m, pvals = len(pvals), np.asarray(pvals)
    ind = np.argsort(pvals)
    rev_ind = np.argsort(ind)
    pvals = pvals[ind]

    # Estimate proportion of features that are truly null.
    kappa = np.arange(0, 0.96, 0.01)
    pik = [sum(pvals > k) / (m*(1-k)) for k in kappa]
    cs = UnivariateSpline(kappa, pik, k=3, s=None, ext=0)
    pi0 = float(cs(1.))
    if (verbose):
        print('The estimated proportion of truly null features is %.3f' % pi0)

    """The smoothing step can sometimes converge outside the interval [0, 1].
    This was noted in the published literature at least by Reiss and
    colleagues. There are at least two approaches one could use to
    attempt to fix the issue:
    (1) Set the estimate to 1 if it is outside the interval, which is the
        assumption in the classic FDR method.
    (2) Assume that if pi0 > 1, it was overestimated, and if pi0 < 0, it
        was underestimated. Set to 0 or 1 depending on which case occurs.

    Here we have chosen the first option, since it is the more conservative
    one of the two.
    """
    if (pi0 < 0 or pi0 > 1):
        pi0 = 1
        print('Smoothing estimator did not converge in [0, 1]')

    # Compute the q-values.
    qvals = np.zeros(np.shape(pvals))
    qvals[-1] = pi0*pvals[-1]
    for i in range(m-2, -1, -1):
        qvals[i] = min(pi0*m*pvals[i]/float(i+1), qvals[i+1])

    # Test which p-values are significant.
    significant = np.zeros(np.shape(pvals), dtype='bool')
    significant[ind] = qvals <= alpha

    """Order the q-values according to the original order of the p-values."""
    qvals = qvals[rev_ind]
    return significant, qvals


def bh_correction(pvals, alpha=0.05):
    """Function for correcting p-values using BH procedure.
        https://www.jstor.org/stable/2346101?seq=1
        :param pvals: numpy array of p-values
        :param alpha: desired FDR
        :return:
        significant - An array of flags indicating which p-values are significant.
        qvals - Q-values corresponding to the p-values.
        """
    results = multipletests(pvals,alpha,method='fdr_bh')
    significant = results[0]
    pvals_corrected = results[1]

    return significant, pvals_corrected

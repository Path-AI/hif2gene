import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import warnings


from sklearn import preprocessing
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

from statsmodels.distributions.empirical_distribution import ECDF #empirical cumulative distribution function
from scipy.special import chdtrc as chi2_cdf
from scipy.stats import pearsonr
from scipy.stats import spearmanr

import scipy.cluster.hierarchy as sch
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy import interpolate
from scipy.signal import argrelextrema



def dfun(u, v):
    C = spearmanr(u,v)
    return 1 - np.abs(C[0])

def cluster_features(scaled_df,t = None, return_linkage = False):
    '''
Input: scaled_df: normalized features dataframe.
       t: cut off for dendrogram.
Output: returns cluster assignments for features.
    '''
    d = sch.distance.pdist(scaled_df.values.transpose(),dfun)
    L = sch.linkage(d, method='ward')
    plt.figure(figsize=(10, 7))
    D = dendrogram(L,
                   orientation='top',
                   distance_sort='descending',
                   show_leaf_counts=True)

    cuts = np.linspace(0.02, 0.5, 1000) * np.max(np.stack(D['dcoord']))
    clusters_n = np.zeros_like(cuts)
    for idx, c in enumerate(cuts):
        ind = sch.fcluster(L, t=c, criterion='distance')
        clusters_n[idx] = np.max(ind)
        if np.max(ind) == 1:
            break

    clusters_n = clusters_n[:idx]
    cuts = cuts[:idx]
    fig, ax = plt.subplots(1)
    kn = KneeLocator(
        cuts,
        clusters_n,
        curve='convex',
        direction='decreasing',
        interp_method='polynomial',
    )
    if t is None:
        t = kn.knee
        ax.plot([kn.knee, kn.knee], [0, np.max(clusters_n)])
    else:
        ax.plot([t, t], [0, np.max(clusters_n)])

    ax.plot(cuts, clusters_n)
    print(t)
    clusters = sch.fcluster(L, t=t, criterion='distance') - 1
    if return_linkage:
        return clusters, L, D
    else:
        return clusters







def combine_p_values(clusters,p_vals,scaled_df,p_val_col = 'Pv.inter'):
    '''
Input: clusters: vector of cluster labels for each feature; 
       p_vals : dataframe of pvalues associating features with some outcome variable.
       scaled_df : dataframe of normalized features used for association
       p_val_col : column to use from p_vals dataframe

Output: Combined p-values for each cluster of features    
    '''
    ebm_p_values = []
    cluster_ids = list(set(clusters))
    for k in cluster_ids:
        # collect p-values for each cluster of features
        features = scaled_df.columns[clusters == k]
        pvals = p_vals.loc[features][p_val_col].values
        # check if there are more than one p-value to combine
        if len(pvals)>1:
        # get covariance matrix for columns
            DataMatrix = scaled_df.loc[:,features].transpose().values
            # compute combined p-value
            ebm_p_values.append(EmpiricalBrownsMethod(DataMatrix, np.array(pvals)))
        # if the cluster is a singleton use the p-value of the single feature
        elif len(pvals)>0:
            ebm_p_values.append(pvals[0])
        else:
        # if none of the features have a p-value return nan for this hypothesis
            print(features)
            ebm_p_values.append(np.nan)
    return np.array(ebm_p_values)


def scale_norm_transform_df(features):
    '''
Input: features: pandas dataframe of features
Output: return pandas dataframe of float type normalized features 
    '''
    # Only keeps float type features
    names = features.columns[features.dtypes==float]
    scaler = preprocessing.StandardScaler()
    scaled_df = scaler.fit_transform(features.loc[:,features.dtypes==float].values.astype(float))
    scaled_df = pd.DataFrame(scaled_df, columns=names)
    return scaled_df


def EmpiricalBrownsMethodArray(data_matrix, p_values_list, extra_info = False):
    '''
    EBM for combining multiple arrays of p-values assumes each array is drawn from the same datamatrix
    :param data_matrix: m x n numpy array each of m rows representing a variable and each of n columns representing a sample.
    :param p_values: list of p-value vectors to combine
    :param extra_info: set verbosity
    :return:
    '''
    # compute this once
    covar_matrix = CalculateCovariances(data_matrix)
    CombinedPValues = []
    for p_values in p_values_list:
        CombinedPValues.append(CombinePValues(covar_matrix,p_values,extra_info))
    return np.asarray(CombinedPValues)


'''
Emperical Brown's Method described here:
https://academic.oup.com/bioinformatics/article/32/17/i430/2450768
Code from William Poole @ https://github.com/IlyaLab/CombiningDependentPvaluesUsingEBM
'''
def EmpiricalBrownsMethod(data_matrix, p_values, extra_info = False):
    '''
    EBM for combining p-values
    :param data_matrix: m x n numpy array each of m rows representing a variable and each of n columns representing a sample.
    :param p_values: p-values to combine
    :param extra_info: set verbosity
    :return:
    '''
    covar_matrix = CalculateCovariances(data_matrix)
    return CombinePValues(covar_matrix, p_values, extra_info)


def TransformData(data_vector):
    '''
    :param data_vector
    :return centered and scaled
    '''
    m = np.mean(data_vector)
    sd = np.std(data_vector)
    s = [(d-m)/sd for d in data_vector]
    W = lambda x: -2*np.log(ECDF(s)(x))
    return np.array([W(x) for x in s])


def CalculateCovariances(data_matrix):
    '''
Input: An m x n data matrix with each of m rows representing a variable and each of n columns representing a sample. Should be of type numpy.array.
       Note: Method does not deal with missing values within the data.
Output: An m x m matrix of pairwise covariances between transformed raw data vectors
'''
    transformed_data_matrix = np.array([TransformData(f) for f in data_matrix])
    covar_matrix = np.cov(transformed_data_matrix)

    return covar_matrix
    
def CombinePValues(covar_matrix, p_values, extra_info = False):
    '''
Input: A m x m numpy array of covariances between transformed data vectors and a vector of m p-values to combine.
Output: A combined P-value. 
        If extra_info == True: also returns the p-value from Fisher's method, the scale factor c, and the new degrees of freedom from Brown's Method

    '''
    m = int(covar_matrix.shape[0])
    #print "m", m
    df_fisher = 2.0*m
    Expected = 2.0*m
    cov_sum = 0
    for i in range(m):
        for j in range(i+1, m):
            cov_sum += covar_matrix[i, j]
    
    #print "cov sum", cov_sum
    Var = 4.0*m+2*cov_sum
    c = Var/(2.0*Expected)
    df_brown = 2.0*Expected**2/Var
    if df_brown > df_fisher:
        df_brown = df_fisher
        c = 1.0

    x = 2.0*sum([-np.log(p) for p in p_values])
    #print "x", x
    p_brown = chi2_cdf(df_brown, 1.0*x/c)
    p_fisher = chi2_cdf(df_fisher, 1.0*x)
    
    if extra_info:
        return p_brown, p_fisher, c, df_brown
    else:
        return p_brown

def KostsMethod(data_matrix, p_values, extra_info = False):
    '''
Input: An m x n data matrix with each of m rows representing a variable and each of n columns representing a sample. Should be of type numpy.array
       A vector of m P-values to combine. May be a list or of type numpy.array.
Output: A combined P-value using Kost's Method.
        If extra_info == True: also returns the p-value from Fisher's method, the scale factor c, and the new degrees of freedom from Brown's Method
    '''
    covar_matrix = CalculateKostCovariance(data_matrix)
    return CombinePValues(covar_matrix, p_values, extra_info = extra_info)
    

def KostPolyFit(cor):
    '''
Input correlation between two n x n data vectors.
Output: Kost's approximation of the covariance between the -log cumulative distributions. This is calculated with a cubic polynomial fit.
    '''
    a1, a2, a3 = 3.263, .710, .027 #Kost cubic coeficients
    return a1*cor+a2*cor**2+a3*cor**3
    
def CalculateKostCovariance(data_matrix):
    '''
Input: An m x n data matrix with each of m rows representing a variable and each of n columns representing a sample. Should be of type numpy.array.
       Note: Method does not deal with missing values within the data.
Output: An m x m matrix of pairwise covariances between the data vectors calculated using Kost's polynomial fit and numpy's pearson correlation function.
    '''
    m = data_matrix.shape[0]
    covar_matrix = np.zeros((m, m))
    for i in range(m):
        for j in range(i+1, m):
            cor, p_val = pearsonr(data_matrix[i, :], data_matrix[j, :])
            covar = KostPolyFit(cor)
            covar_matrix[i, j] = covar
            covar_matrix[j, i] = covar
    return covar_matrix
    





class KneeLocator(object):
    def __init__(
        self,
        x,
        y,
        S=1.0,
        curve="concave",
        direction="increasing",
        interp_method="interp1d",
        online=False,
    ):
        """
        Once instantiated, this class attempts to find the point of maximum
        curvature on a line. The knee is accessible via the `.knee` attribute.
        :param x: x values.
        :type x: list or array.
        :param y: y values.
        :type y: list or array.
        :param S: Sensitivity, original paper suggests default of 1.0
        :type S: float
        :param curve: If 'concave', algorithm will detect knees. If 'convex', it
            will detect elbows.
        :type curve: string
        :param direction: one of {"increasing", "decreasing"}
        :type direction: string
        :param interp_method: one of {"interp1d", "polynomial"}
        :type interp_method: string
        :param online: Will correct old knee points if True, will return first knee if False
        :type online: bool
        """
        # Step 0: Raw Input
        self.x = np.array(x)
        self.y = np.array(y)
        self.curve = curve
        self.direction = direction
        self.N = len(self.x)
        self.S = S
        self.all_knees = set()
        self.all_norm_knees = set()
        self.online = online

        # Step 1: fit a smooth line
        if interp_method == "interp1d":
            uspline = interpolate.interp1d(self.x, self.y)
            self.Ds_y = uspline(self.x)
        elif interp_method == "polynomial":
            pn_model = PolynomialFeatures(7)
            xpn = pn_model.fit_transform(self.x.reshape(-1, 1))
            regr_model = LinearRegression()
            regr_model.fit(xpn, self.y)
            self.Ds_y = regr_model.predict(
                pn_model.fit_transform(self.x.reshape(-1, 1))
            )
        else:
            warnings.warn(
                "{} is an invalid interp_method parameter, use either 'interp1d' or 'polynomial'".format(
                    interp_method
                )
            )
            return

        # Step 2: normalize values
        self.x_normalized = self.__normalize(self.x)
        self.y_normalized = self.__normalize(self.Ds_y)

        # Step 3: Calculate the Difference curve
        self.x_normalized, self.y_normalized = self.transform_xy(
            self.x_normalized, self.y_normalized, self.direction, self.curve
        )
        # normalized difference curve
        self.y_difference = self.y_normalized - self.x_normalized
        self.x_difference = self.x_normalized.copy()

        # Step 4: Identify local maxima/minima
        # local maxima
        self.maxima_indices = argrelextrema(self.y_difference, np.greater_equal)[0]
        self.x_difference_maxima = self.x_difference[self.maxima_indices]
        self.y_difference_maxima = self.y_difference[self.maxima_indices]

        # local minima
        self.minima_indices = argrelextrema(self.y_difference, np.less_equal)[0]
        self.x_difference_minima = self.x_difference[self.minima_indices]
        self.y_difference_minima = self.y_difference[self.minima_indices]

        # Step 5: Calculate thresholds
        self.Tmx = self.y_difference_maxima - (
            self.S * np.abs(np.diff(self.x_normalized).mean())
        )

        # Step 6: find knee
        self.knee, self.norm_knee = self.find_knee()

    @staticmethod
    def __normalize(a):
        """normalize an array
        :param a: The array to normalize
        :type a: array
        """
        return (a - min(a)) / (max(a) - min(a))

    @staticmethod
    def transform_xy(x, y, direction, curve):
        """transform x and y to concave, increasing based on given direction and curve"""
        # convert elbows to knees
        if curve == "convex":
            x = x.max() - x
            y = y.max() - y
        # flip decreasing functions to increasing
        if direction == "decreasing":
            y = np.flip(y, axis=0)

        if curve == "convex":
            x = np.flip(x, axis=0)
            y = np.flip(y, axis=0)

        return x, y

    def find_knee(self,):
        """This function finds and sets the knee value and the normalized knee value. """
        if not self.maxima_indices.size:
            warnings.warn(
                "No local maxima found in the difference curve\n"
                "The line is probably not polynomial, try plotting\n"
                "the difference curve with plt.plot(knee.x_difference, knee.y_difference)\n"
                "Also check that you aren't mistakenly setting the curve argument",
                RuntimeWarning,
            )
            return None, None

        # placeholder for which threshold region i is located in.
        maxima_threshold_index = 0
        minima_threshold_index = 0
        # traverse the difference curve
        for i, x in enumerate(self.x_difference):
            # skip points on the curve before the the first local maxima
            if i < self.maxima_indices[0]:
                continue

            j = i + 1

            # reached the end of the curve
            if x == 1.0:
                break

            # if we're at a local max, increment the maxima threshold index and continue
            if (self.maxima_indices == i).any():
                threshold = self.Tmx[maxima_threshold_index]
                threshold_index = i
                maxima_threshold_index += 1
            # values in difference curve are at or after a local minimum
            if (self.minima_indices == i).any():
                threshold = 0.0
                minima_threshold_index += 1

            if self.y_difference[j] < threshold:
                if self.curve == "convex":
                    if self.direction == "decreasing":
                        knee = self.x[threshold_index]
                        norm_knee = self.x_normalized[threshold_index]
                    else:
                        knee = self.x[-(threshold_index + 1)]
                        norm_knee = self.x_normalized[-(threshold_index + 1)]

                elif self.curve == "concave":
                    if self.direction == "decreasing":
                        knee = self.x[-(threshold_index + 1)]
                        norm_knee = self.x_normalized[-(threshold_index + 1)]
                    else:
                        knee = self.x[threshold_index]
                        norm_knee = self.x_normalized[threshold_index]

                self.all_knees.add(knee)
                self.all_norm_knees.add(norm_knee)

                # if detecting in offline mode, return the first knee found
                if self.online is False:
                    return knee, norm_knee

        if self.all_knees == set():
            warnings.warn("No knee/elbow found")
            return None, None

        return knee, norm_knee

    def plot_knee_normalized(self,):
        """Plot the normalized curve, the difference curve (x_difference, y_normalized) and the knee, if it exists."""
        import matplotlib.pyplot as plt

        plt.figure(figsize=(6, 6))
        plt.title("Normalized Knee Point")
        plt.plot(self.x_normalized, self.y_normalized, "b", label="normalized curve")
        plt.plot(self.x_difference, self.y_difference, "r", label="difference curve")
        plt.xticks(
            np.arange(self.x_normalized.min(), self.x_normalized.max() + 0.1, 0.1)
        )
        plt.yticks(
            np.arange(self.y_difference.min(), self.y_normalized.max() + 0.1, 0.1)
        )

        plt.vlines(
            self.norm_knee,
            plt.ylim()[0],
            plt.ylim()[1],
            linestyles="--",
            label="knee/elbow",
        )
        plt.legend(loc="best")

    def plot_knee(self,):
        """Plot the curve and the knee, if it exists"""
        import matplotlib.pyplot as plt

        plt.figure(figsize=(6, 6))
        plt.title("Knee Point")
        plt.plot(self.x, self.y, "b", label="data")
        plt.vlines(
            self.knee, plt.ylim()[0], plt.ylim()[1], linestyles="--", label="knee/elbow"
        )
        plt.legend(loc="best")

    # Niceties for users working with elbows rather than knees
    @property
    def elbow(self):
        return self.knee

    @property
    def norm_elbow(self):
        return self.norm_knee

    @property
    def all_elbows(self):
        return self.all_knees

    @property
    def all_norm_elbows(self):
        return self.all_norm_knees
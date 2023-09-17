"""
Various statistical tests for hypothesis testing.
"""


import numpy as np
from scipy.stats import t
from scipy.stats import norm
from scipy.stats import f
from scipy.stats import chi2

# pylint: disable=C0103
# pylint: disable=R0913

def z_test(a, b, alpha=0.05, alter="two-tailed") -> tuple:
    """
    Perform a Z-test for two samples.

    Args:
        a (numpy.ndarray): The first sample.
        b (numpy.ndarray): The second sample.
        alpha (float, optional): The significance level. Default value 0.05.
        alter (str, optional): The alternative hypothesis. Default "two-tailed", 
           choose from [right, left, two-detailed].

    Returns:
        Tuple: (test statstic: float, p_value: float, alpha: float, alternative: str).
    """
    n1, n2, x_bar1, x_bar2, sigma1, sigma2 = _sample_statistics(a, b)

    test_statistic = (x_bar1 - x_bar2) / np.sqrt(sigma1**2/n1 + sigma2**2/n2)

    if alter == "two-tailed":
        p_value = 2 * (1 - norm.cdf(abs(test_statistic)))
    elif alter == 'left':
        p_value =  1 - norm.cdf(abs(test_statistic))

    elif alter == 'right':
        p_value =  norm.cdf(abs(test_statistic))
    else:
        raise ValueError("Invalid value")

    return (test_statistic, p_value, alpha, alter)

def t_independent(a, b, alpha=0.05, alter="two-tailed", eq_varience=True) -> tuple:
    """
    Perform an independent two-sample T-test.

    Args:
        a (numpy.ndarray): The first sample.
        b (numpy.ndarray): The second sample.
        alpha (float, optional): The significance level. Default value 0.05.
        alter (str, optional): The alternative hypothesis. Default "two-tailed", 
           choose from [right, left, two-detailed].

    Returns:
        Tuple: (test statstic: float, p_value: float, alpha: float, alternative: str).
    """

    if eq_varience:
        test_statistic, dof = _t_equal_var(a, b)
    else:
        test_statistic, dof = _t_unequal_var(a, b)

    if alter == "two-tailed":
        p_value = 2 * (1 - t.cdf(abs(test_statistic), df=dof))
    elif alter == 'left':
        p_value = 1 - t.cdf(abs(test_statistic), df=dof)
    elif alter == 'right':
        p_value = t.cdf(abs(test_statistic), df=dof)
    else:
        raise ValueError("Invalud value")


    return (test_statistic, p_value, alpha, alter)

def proportion_test(obs1, n1, obs2, n2, alpha=0.05, method="pooled", alter="two-tailed") -> tuple:
    """
    Perform a two-sample proportion test.

    Args:
        obs1 (int): The number of successes in the first sample.
        n1 (int): The size of the first sample.
        obs2 (int): The number of successes in the second sample.
        n2 (int): The size of the second sample.
        alpha (float, optional): The significance level. Default value = 0.05.
        method (str, optional): The method for calculating test statistic. Defaul value "pooled".
           choose from ['pooled', '']
        alter (str, optional): The alternative hypothesis. Default "two-tailed", 
           choose from [right, left, two-detailed].

    Returns:
        Tuple: (test statstic: float, p_value: float, alpha: float, alternative: str).
    """
    # DELETE - update - in each file can be 0 and 1, where one ->  and 0 -> fail, then find the \
    # length and only the prob of a
    # and b (numerator only 1 and denom all)
    # n1 = len(a)
    # n2 = len(b)
    # p1 = np.count_nonzero(a==1)
    # p2 = np.count_nonzero(a==1)

    p1 = obs1 / n1
    p2 = obs2 / n2

    if method == "pooled":
        test_statistic = _proportion_test_p_pooled(p1, n1, p2, n2)
    elif method == "unpooled":
        test_statistic = _proportion_test_p_unpooled(p1, n1, p2, n2)
    else:
        raise ValueError("Invalud value")

    if alter == "two-tailed":
        p_value = 2 * (1 - norm.cdf(abs(test_statistic)))

    elif alter == 'right':
        p_value =  norm.cdf(abs(test_statistic))

    elif alter == 'left':
        p_value =  1 - norm.cdf(abs(test_statistic))
    else:
        raise ValueError("Invalud value")

    return (test_statistic, p_value, alpha, alter)



def paired_test(a, b, alpha=0.05, alter='two-tailed') -> tuple:
    """
    Perform a paired T-test.

    Args:
        a (numpy.ndarray): The first paired sample.
        b (numpy.ndarray): The second paired sample.
        alpha (float, optional): The significance level. Default value 0.05.
        alter (str, optional): The alternative hypothesis. Default "two-tailed", 
           choose from [right, left, two-detailed].

    Returns:
        Tuple: (test statstic: float, p_value: float, alpha: float, alternative: str).
    """
    n1 = len(a)
    n2 = len(b)

    if n1 != n2:
        raise ValueError("The two input arrays must have the same length.")
    dof = n1 - 1
    diff = a - b
    d_bar = np.mean(diff)
    std = np.std(diff, ddof=1)
    test_statistic = d_bar / (std/np.sqrt(n1))

    if alter == 'two-tailed':
        p_value = 2 * (1 - t.cdf(abs(test_statistic), dof))
    elif alter == 'left':
        p_value = 1 - t.cdf(abs(test_statistic), dof)
    elif alter == 'right':
        p_value = t.cdf(abs(test_statistic), dof)
    else:
        raise ValueError("Invalud value")

    return (test_statistic, p_value, alpha, alter)

def f_test(a, b, alpha=0.05) -> tuple:
    """
    Perform F-test for comparing variances of two samples.

    Args:
        a (numpy.ndarray): The first sample.
        b (numpy.ndarray): The second sample.
        alpha (float, optional): The significance level. Default value 0.05.
    Returns:
        Tuple: (test statstic: float, p_value: float, alpha: float).
    """
    var1 = np.var(a, ddof=1) if np.var(a, ddof=1) < np.var(b, ddof=1) else np.var(b, ddof=1)
    var2 =  np.var(a, ddof=1) if np.var(a, ddof=1) > np.var(b, ddof=1) else np.var(b, ddof=1)
    dof1 = len(a) - 1
    dof2 = len(b) -1
    test_statistic = var1 / var2
    p_value = 1 - f.cdf(test_statistic, dof1, dof2)

    return (test_statistic, p_value, alpha)

def chi2_conting(observed, alpha=0.05) -> tuple:
    """
    Perform a chi-squared contingency test on 2x2 table.

    Args:
        observed (numpy.ndarray): The observed 2x2 table.

    Returns:
        Tuple: (test statstic: float, p_value: float, alpha: float).
    """
    # “fourfold table” or “2 x 2 contingency table”

    test_statistic = _chi2_calc(observed)
    p_value = 1 - chi2.cdf(test_statistic, df=1)

    return (test_statistic, p_value, alpha)

def _sample_statistics(a, b) -> tuple:
    """Calculate the length, mean and standard deviation of two samples."""

    n1 = len(a)
    n2 = len(b)
    x_bar1 = a.mean()
    x_bar2 = b.mean()
    std1 = a.std(ddof=1)
    std2 = b.std(ddof=1)

    return (n1, n2, x_bar1, x_bar2, std1, std2)

def _t_equal_var(a, b):
    """Calculate the test statistic of 2 equal variences t-distribution."""

    n1, n2, x_bar1, x_bar2, std1, std2 = _sample_statistics(a, b)
    pooled_std = np.sqrt(((n1-1)*(std1**2) + (n2-1)*(std2**2)) / (n1+n2-2))
    test_statistic = (x_bar1-x_bar2) / (pooled_std*np.sqrt((1/n1)+(1/n2)))
    dof = n1 + n2 - 2
    return test_statistic, dof

def _t_unequal_var(a, b) -> tuple:
    """Calculate the test statistic of 2 unequal variences t-distribution."""

    n1, n2, x_bar1, x_bar2, std1, std2 = _sample_statistics(a, b)
    numerator = (std1**2/n1 + std2**2/n2)**2
    denominator = ((std1**2/n1)**2/(n1-1)) + ((std2**2/n2)**2/(n2-1))

    dof = np.floor(numerator / denominator)
    test_statistic =  (x_bar1 - x_bar2) / np.sqrt(std1**2/n1 + std2**2/n2)

    return (test_statistic, dof)

def _proportion_test_p_pooled(p1, n1, p2, n2) -> float:
    """Calculate the test statistic of 2 pooled proportions sample."""

    p_pooled = (p1*n1 + p2*n2) / (n1+n2)
    test_statistic = (p1-p2) / np.sqrt(p_pooled*(1 - p_pooled) * ((1/n1)+(1/n2)))
    return test_statistic

def _proportion_test_p_unpooled(p1, n1, p2, n2) -> float:
    """Calculate the test statistic of 2 unpooled proportions sample."""

    test_statistic = (p1-p2) / np.sqrt(p1*(1-p1)/n1 + p2*(1-p2)/n2)
    return test_statistic

def _chi2_calc(observed) -> float:
    """Calculate the test statistic of the 2*2 chi2 contingency table."""

    det = observed[0, 0] * observed[1, 1] - observed[1, 0] * observed[0, 1]
    test_statistic = (
        det / observed[0].sum()
        * det / observed[1].sum()
        * observed.sum() / (observed[:, 0].sum() * observed[:, 1].sum())
    )
    return test_statistic

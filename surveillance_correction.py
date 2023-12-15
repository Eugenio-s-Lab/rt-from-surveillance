import numpy as np
import pandas as pd
import scipy.linalg
from epyestim import bagging_r
from typing import Optional,List,Iterable
from datetime import date
import warnings


def Reproduction_Operator(colocation_matrix: pd.DataFrame, departments_population:pd.Series) -> pd.DataFrame:
    """
    Compute the reproduction operator from colocation matrix and populations series
    
    """
    return colocation_matrix.mul(departments_population,axis=0)
 
    
def left_eigenvector(R: pd.DataFrame)->pd.Series:
    """
    Compute v_star, left eigenvector associated to the spectral radius of the reproduction operator

    """
    w,vl,vr=scipy.linalg.eig(R.to_numpy(),left=True)
    v_star=vl[:,w.argmax()]/np.sum(vl[:,w.argmax()])
    v_star=pd.Series(v_star.real)  
    v_star=v_star.set_axis(R.columns)
    return v_star



def infections_corrected(new_infections_departments: pd.DataFrame, v_star: pd.Series)->pd.Series:
    """
    Compute corrected infections reweighting incident infections in communities

    """
    N=v_star.shape[0]
    new_infections_corr = new_infections_departments.copy()
    new_infections_corr.iloc[:,1:]=new_infections_corr.iloc[:,1:].mul(v_star*N,axis=1)
    new_infections_corr['total'] = new_infections_corr.iloc[:,1:].sum(axis=1) 
    new_infections_corr['date']=pd.to_datetime(new_infections_corr['date'])
    new_infections_corr=new_infections_corr.set_index('date')['total']
    return new_infections_corr


def infections_sum(new_infections_departments: pd.DataFrame)->pd.Series:
    """
    Compute global incident infections from community-level ones
    
    """
    new_infections = new_infections_departments.copy()
    new_infections['total'] = new_infections.iloc[:,1:].sum(axis=1) 
    new_infections['date']=pd.to_datetime(new_infections['date'])
    new_infections = new_infections.set_index('date')['total']
    return new_infections


def corrected_rt_from_surveillance(
        colocation_matrix: pd.DataFrame,
        departments_population: pd.Series,
        new_infections_departments: pd.DataFrame,
        gt_distribution: np.ndarray,
        delay_distribution: np.ndarray,
        correction: Optional[bool] =True,
        a_prior: Optional[float]=1,
        b_prior: Optional[float]=5,
        smoothing_window: Optional[int] = None,
        r_window_size: Optional[int] = None,
        r_interval_dates: Optional[List[date]] = None,
        n_samples: int = 100,
        quantiles: Iterable[float] = (0.025, 0.5, 0.975),
        auto_cutoff: bool = True
                                    )->pd.DataFrame: 
    """
    Estimate the corrected reproduction ratio.

    Parameters
    ----------
    colocation_matrix : pd.DataFrame
        colocation matrix.
    departments_population : pd.Series
        populations of departments.
    new_infections_departments : pd.DataFrame
        time series of infection numbers.
    gt_distribution : np.ndarray
        the generation time distribution.
    delay_distribution : np.ndarray
        discretised delay distribution (delay_distrb[j] = probability
        that an infection gets detected with a delay of j days). Array must contain no zero values.
    correction : Optional[bool], optional
        choose whether to correct estimates or not. The default is True.
    a_prior : Optional[float], optional
        prior for the Gamma shape parameter for R. The default is 1.
    b_prior : Optional[float], optional
        prior for the Gamma scale parameter for R. The default is 5.
    smoothing_window : Optional[int], optional
        width of moving time frame (in days) for local regression. The default is None.
    r_window_size : Optional[int], optional
        Size of the rolling window. The default is None.
    r_interval_dates : Optional[List[date]], optional
        boundaries of the intervals for which R should be estimated. The default is None.
    n_samples : int, optional
        number of samples to compute quantiles. The default is 100.
    quantiles : Iterable[float], optional
        desired quantiles. The default is (0.025, 0.5, 0.975).
    auto_cutoff : bool, optional
        The default is True.

    Returns
    -------
    estimated_R : pd.DataFrame
        The estimated reproduction number.

    """
    
    """ Verify that inputs have correct type, dimension and indices. """
    assert type(colocation_matrix)==pd.core.frame.DataFrame, 'colocation_matrix must be a dataframe'
    assert type(departments_population)==pd.core.series.Series, 'departments_population must be a series'
    assert type(new_infections_departments)==pd.core.frame.DataFrame, 'new_infections_departments must be a dataframe'
    assert new_infections_departments.columns[0]=='date', 'first column of new_infections_departments must be \'date\''
    assert colocation_matrix.shape[0]==colocation_matrix.shape[1], 'colocation_matrix must be a square dataframe'
    assert departments_population.shape[0]==colocation_matrix.shape[0]==new_infections_departments.shape[1]-1, 'Shapes must match between colocation_matrix, departments_population and columns of new_infections_departments (except \'date\')'
    assert len(set(departments_population.index)-set(colocation_matrix.index))==0 and len(set(departments_population.index)-set(new_infections_departments.columns[1:]))==0 and len(set(new_infections_departments.columns[1:])-set(colocation_matrix.index))==0, 'Indices and columns of colocation_matrix, indices of departments_population and columns of new_infections_departments (except \'date\') must match'
    
    print('Estimating...',end='', flush=True)
    
    R = Reproduction_Operator(colocation_matrix,departments_population)
    v_star = left_eigenvector(R)
    
    if correction==True:
        
        new_infections=infections_corrected(new_infections_departments, v_star)
    
        estimated_R = bagging_r(confirmed_cases=new_infections,
                                 gt_distribution=gt_distribution,
                                 delay_distribution=delay_distribution,
                                 a_prior=a_prior,
                                 b_prior=b_prior,
                                 smoothing_window=smoothing_window,
                                 r_window_size=r_window_size,
                                 r_interval_dates=r_interval_dates,
                                 n_samples=n_samples,              
                                 quantiles=quantiles,
                                 auto_cutoff=auto_cutoff)
    
    else:
        
        new_infections=infections_sum(new_infections_departments)
    
        estimated_R = bagging_r(confirmed_cases=new_infections,
                                 gt_distribution=gt_distribution,
                                 delay_distribution=delay_distribution,
                                 a_prior=a_prior,
                                 b_prior=b_prior,
                                 smoothing_window=smoothing_window,
                                 r_window_size=r_window_size,
                                 r_interval_dates=r_interval_dates,
                                 n_samples=n_samples,              
                                 quantiles=quantiles,
                                 auto_cutoff=auto_cutoff)
    
    print("\r", end='', flush=True)
    
    return estimated_R
    
    
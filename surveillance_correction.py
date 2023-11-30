import numpy as np
import pandas as pd
import scipy.linalg
import epyestim.covid19 as covid19
from epyestim import bagging_r
from typing import Optional,List,Iterable
from datetime import date
import warnings

def Reproduction_Operator(colocation_matrix: pd.DataFrame,departments_population:pd.Series)->pd.DataFrame:
    return colocation_matrix.mul(departments_population,axis=0)

def left_eigenvector(R: pd.DataFrame)->pd.Series:
    w,vl,vr=scipy.linalg.eig(R.to_numpy(),left=True)
    v_star=vl[:,w.argmax()]/np.sum(vl[:,w.argmax()])
    v_star=pd.Series(v_star.real)  
    v_star.set_axis(R.columns,inplace=True)
    return v_star

def infections_corrected(new_infections_departments: pd.DataFrame, v_star: pd.Series)->pd.Series:
    N=v_star.shape[0]
    new_infections_corr = new_infections_departments.copy()
    new_infections_corr.iloc[:,1:]=new_infections_corr.iloc[:,1:].mul(v_star*N,axis=1)
    new_infections_corr['total'] = new_infections_corr.iloc[:,1:].sum(axis=1) 
    new_infections_corr['date']=pd.to_datetime(new_infections_corr['date'])
    new_infections_corr=new_infections_corr.set_index('date')['total']
    return new_infections_corr

def infections_sum(new_infections_departments: pd.DataFrame)->pd.Series:
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
    
    
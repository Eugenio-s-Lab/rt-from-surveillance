from surveillance_correction import corrected_rt_from_surveillance

import pandas as pd
import epyestim.covid19 as covid19
import matplotlib.pyplot as plt

# Colocation data
colocation_matrix=pd.read_csv('Matrix_Example',index_col=0)    
colocation_matrix.index=colocation_matrix.index.astype(str)

# Populations of departments
departments_population=pd.read_csv('Populations_Example',index_col=0)
departments_population=pd.Series(departments_population.iloc[:,0])
departments_population.index=departments_population.index.astype(str)

# New infections in departments
new_infections_departments = pd.read_csv('Cases_Example', index_col=0)

# Insert your generation interval distribution. As an example, here we set a generation interval distribution for COVID-19.
si_distrb = covid19.generate_standard_si_distribution()
# Insert your distribution of delay from infection to reporting. As an example, here we set a distribution of delay for COVID-19
delay_distrb = covid19.generate_standard_infection_to_reporting_distribution()

estimated_R=corrected_rt_from_surveillance(colocation_matrix,
                                     departments_population,
                                     new_infections_departments,
                                     si_distrb,
                                     delay_distrb,
                                     correction=False,
                                     smoothing_window=7,
                                     r_window_size=7)

estimated_R_corr=corrected_rt_from_surveillance(colocation_matrix,
                                     departments_population,
                                     new_infections_departments,
                                     si_distrb,
                                     delay_distrb,
                                     correction=True,
                                     smoothing_window=7,
                                     r_window_size=7)

# plot results
fig,ax=plt.subplots(figsize=(10,5))
ax.plot(estimated_R.index[7:],estimated_R['R_mean'][7:],color='orchid',label='EpiEstim')
ax.plot(estimated_R_corr.index[7:],estimated_R_corr['R_mean'][7:],color='tab:blue',label=r'EpiEstim$^{corr}$')
ax.tick_params(axis='x',labelsize=15,rotation=90)
ax.tick_params(axis='y',labelsize=15)
ax.legend(fontsize=15)
fig.tight_layout()
plt.show()


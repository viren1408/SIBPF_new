import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read the pulsar candidates  This is a dummy data . This step will also be generalized for multiple observation I am still figuring out how to directly script pybdsf in this 
candidates_day_1 = pd.read_csv('Pulsar_candidates_G033.0-5.0.csv')
candidates_day_2 = pd.read_csv('Pulsar_candidates_G033.0-5.0 _1.csv')
candidates_day_3 = pd.read_csv('Pulsar_candidates_G033.0-5.0 _2.csv')
candidates_day_4 = pd.read_csv('Pulsar_candidates_G033.0-5.0 _3.csv')

# Merge the DataFrames based on 'Xposn' and 'Yposn' columns 
merged_df = candidates_day_1.merge(candidates_day_2, on=[' Xposn', ' Yposn'], suffixes=('_day1', '_day2'))
merged_df = merged_df.merge(candidates_day_3, on=[' Xposn', ' Yposn'])
merged_df = merged_df.merge(candidates_day_4, on=[' Xposn', ' Yposn'], suffixes=('_day3', '_day4'))

# Initialize lists to store data
flux_values = []
error_values = []

# Loop through candidates
days = [1,2,3,4] #This is a dummy time axes and will be adjusted according to our observation , will of course be generalized once data is available 
for i in range(len(merged_df)):
    # Extract flux and error values for each day
    flux = [merged_df[ ' Total_flux_day1'][i],merged_df[ ' Total_flux_day2'][i],merged_df[ ' Total_flux_day3'][i],merged_df[ ' Total_flux_day4'][i]] 
    error = [merged_df[ ' E_Total_flux_day1'][i],merged_df[ ' E_Total_flux_day2'][i],merged_df[ ' E_Total_flux_day3'][i],merged_df[ ' E_Total_flux_day4'][i]] 

    # Calculate weighted mean (sbar_i) and weighted standard deviation
    sbar_i = np.sum(np.array(flux) / np.array(error) ** 2) / np.sum(1 / np.array(error) ** 2)
    weighted_std = np.sqrt(1 / np.sum(1 / np.array(error) ** 2))

    #chi-squared 
    chi_squared = np.sum(((np.array(flux) - sbar_i) / np.array(error)) ** 2)

    #degrees of freedom
    degrees_of_freedom = len(flux) - 1 

    #reduced chi-squared (chi_sq_r)
    chi_sq_r = chi_squared / degrees_of_freedom

    # Append chi-squared and reduced chi-squared to lists
    flux_values.append(sbar_i)
    error_values.append(weighted_std)

    print(f"Chi-squared for Candidate {i+1}: {chi_squared}")
    print(f"Reduced Chi-squared for Candidate {i+1}: {chi_sq_r}")
    plt.scatter(days,flux)
    plt.errorbar(days,flux,yerr = error,alpha=0.7,ecolor='black',capsize=2,ls = 'none')
    plt.plot(days,flux, linestyle='-', color='red')
    plt.title(f"$\\chi^2$ for candidate {i+1}: {chi_squared}")
    plt.grid()
    plt.xlabel('Days')
    plt.ylabel('Flux(Jy)')
    plt.show()
# Now you have lists containing the weighted mean (flux_values) and weighted standard deviation (error_values)
# for each candidate, as well as the chi-squared and reduced chi-squared values for each candidate.
# You can further analyze or plot these values as needed.

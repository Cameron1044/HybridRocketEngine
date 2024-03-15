import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import os
import matplotlib.pyplot as plt

# Define the fitting equation
def fit_equation(x, E1, E2, E3, E4, E5):    
    return E1 + E2*x + E3*x**2 + E4*x**3 + E5*x**4

initial_guesses = [6.7556e4, 5.4373e1, 0, 0, 0]

# Prepare a DataFrame to store Q values for each pressure level
e_values_df = pd.DataFrame(columns=['Pressure', 'E1', 'E2', 'E3', 'E4', 'E5'])

# Folder containing the data files
data_folder = 'src/Model/NistDataProc/data'
plt.figure(figsize=(10, 6))
# Iterate over the files in increments of 100 PSI from 2000 to 3500
for pressure in range(2000, 3600, 100):
    file_path = os.path.join(data_folder, f'fluid{pressure}.txt')
    dataRaw = pd.read_csv(file_path, sep='\t')
    
    # Convert density from lbm/ft^3 to kg/m^3
    dataRaw['Cp (J/mol*K)'] = dataRaw['Cp (J/mol*K)'] * 1000
    
    # Filter rows where the temperature values are whole numbers
    data = dataRaw[dataRaw['Temperature (K)'].apply(lambda x: x == round(x))]
    
    x_data = data['Temperature (K)']
    y_data = data['Cp (J/mol*K)']
    
    # Perform the curve fitting
    try:
        popt, pcov = curve_fit(fit_equation, x_data, y_data, p0=initial_guesses, maxfev=10000)
    except RuntimeError:
        popt, pcov = curve_fit(fit_equation, x_data, y_data, p0=initial_guesses, maxfev=50000)
    
    # Use the current fitted parameters as the initial guesses for the next pressure level
    initial_guesses = popt
    print(initial_guesses)
    
    # Append the fitted parameters to the DataFrame
    temp_df = pd.DataFrame({'Pressure': [pressure], 'E1': [popt[0]], 'E2': [popt[1]], 'E3': [popt[2]], 'E4': [popt[3]], 'E5': [popt[4]]})
    e_values_df = pd.concat([e_values_df, temp_df], ignore_index=True)
    print(e_values_df)
    plt.scatter(x_data, y_data, label='Original Data', color='blue')
    temperature_range = np.linspace(min(x_data), max(x_data), 500)
    plt.plot(temperature_range, fit_equation(temperature_range, *popt), label='Fitted Curve', color='red', linestyle='--')
    plt.title('Density vs. Temperature with Fitted Curve')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Density (kg/m³)')
    plt.legend()
    plt.grid(True)
plt.show()

# Save the DataFrame of Q values to a CSV file
e_values_csv_path = 'src/Model/NistDataProc/data/E_values.csv'
e_values_df.to_csv(e_values_csv_path, index=False)
print(f'Q values saved to {e_values_csv_path}')
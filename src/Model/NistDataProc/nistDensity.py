import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import os
import matplotlib.pyplot as plt

# Define the fitting equation
def fit_equation(x, Q1, Q2, Q3, Q4):
    return 44.013 / (Q2**(1 + (1 - x / Q3)**Q4) / Q1)

# Initial guesses for the parameters for the first file
initial_guesses = [1.999391250, 0.230325058, 335.475601882, 0.270838709]

# Prepare a DataFrame to store Q values for each pressure level
q_values_df = pd.DataFrame(columns=['Pressure', 'Q1', 'Q2', 'Q3', 'Q4'])

# Folder containing the data files
data_folder = 'src/Model/NistDataProc/data'
plt.figure(figsize=(10, 6))
# Iterate over the files in increments of 100 PSI from 2000 to 3500
for pressure in range(2000, 3600, 100):
    file_path = os.path.join(data_folder, f'fluid{pressure}.txt')
    dataRaw = pd.read_csv(file_path, sep='\t')
    
    # Convert density from lbm/ft^3 to kg/m^3
    dataRaw['Density (kg/m3)'] = dataRaw['Density (lbm/ft3)'] * 16.0185
    
    # Filter rows where the temperature values are whole numbers
    data = dataRaw[dataRaw['Temperature (K)'].apply(lambda x: x == round(x))]
    
    x_data = data['Temperature (K)']
    y_data = data['Density (kg/m3)']
    
    # Perform the curve fitting
    try:
        popt, pcov = curve_fit(fit_equation, x_data, y_data, p0=initial_guesses, maxfev=10000)
    except RuntimeError:
        popt, pcov = curve_fit(fit_equation, x_data, y_data, p0=initial_guesses, maxfev=50000)
    
    # Use the current fitted parameters as the initial guesses for the next pressure level
    initial_guesses = popt
    print(initial_guesses)
    
    # Append the fitted parameters to the DataFrame
    temp_df = pd.DataFrame({'Pressure': [pressure], 'Q1': [popt[0]], 'Q2': [popt[1]], 'Q3': [popt[2]], 'Q4': [popt[3]]})
    q_values_df = pd.concat([q_values_df, temp_df], ignore_index=True)
    print(q_values_df)
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
q_values_csv_path = 'src/Model/NistDataProc/data/Q_values.csv'
q_values_df.to_csv(q_values_csv_path, index=False)
print(f'Q values saved to {q_values_csv_path}')
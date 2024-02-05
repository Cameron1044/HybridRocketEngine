import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


df = pd.read_csv("src/Model/Test_data/01_30_2024_Collected_data/CSV_Files/Burn_2.csv")

# Time (s)
# Thrust (lbf)
# Combustion Chamber Temperature 1 (R)
# N20 Tank Temperature 1 (R)
# N20 Tank Pressure (psi)
# Combustion Chamber Pressure (psi)
# Combustion Chamber Temperature 2 (R)
# N20 Tank Temperature 2 (R)

#find index of max thrust
max_thrust = df['Thrust (lbf)'].max()
#cut off dataframe from index 0 to max_index
max_index = df['Thrust (lbf)'].idxmax()
# df = df.iloc[max_index+3:-1]

plt.figure()
plt.plot(df['Time (s)'], df['Thrust (lbf)'] + 171.165, linewidth=0.1)
plt.xlabel('Time (s)')
plt.ylabel('Thrust (lbf)')
plt.title('Thrust vs Time')

plt.figure()
plt.plot(df['Time (s)'], df['N20 Tank Pressure (psi)'], linewidth=0.1)
plt.xlabel('Time (s)')
plt.ylabel('Pressure (psi)')
plt.title('Tank Pressure vs Time')

plt.figure()
plt.plot(df['Time (s)'], df['Combustion Chamber Pressure (psi)'], linewidth=0.1)
plt.xlabel('Time (s)')
plt.ylabel('Pressure (psi)')
plt.title('Chamber Pressure vs Time')

plt.show()
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from nptdms import TdmsFile

# Read the TDMS file
tdms_file = TdmsFile.read('src/Model/Test_data/Burn_1.tdms')

# Access the specific group
group = tdms_file['Characteristic_Fire_01_26_2024']

# Extracting data for each channel in the group
data = {channel.name: channel.data for channel in group.channels()}

# Creating a Pandas DataFrame from the extracted data
df = pd.DataFrame(data)

# Plotting Force against its index
# Assuming that the slope of the pressure ducer is 5000 PSI / 4 V with a weird off set of 1250 PSI
plt.grid(which='both', linestyle='--', linewidth=0.5)
plt.plot(df.index, (df['Temperature_1']), label='Force',linewidth=0.5)
plt.axhline(y=0, linewidth=2, linestyle='--', color='k')
plt.xlabel('Index (Assumed Equal Time Interval)')
plt.ylabel('Tank Pressure')
plt.title('N20 PSI vs. Index')
plt.legend()
plt.show()
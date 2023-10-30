import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression

# Read the CSV file
df = pd.read_csv('src/GDL/results.csv')
def FuelProperties(OF, Pc):

    # Calculate the Euclidean distance between the given OF, Pc and all rows in the DataFrame
    distances = np.sqrt((df['OF'] - OF)**2 + (df['Pc'] - Pc)**2)
    
    # Find the index of the minimum distance
    idx = distances.idxmin()
    
    # Extract the corresponding T, k, and M values
    T_pred = df.loc[idx, 'T']
    k_pred = df.loc[idx, 'k']
    M_pred = df.loc[idx, 'M']
    
    return T_pred, k_pred, M_pred

print(FuelProperties(100, 0))

# Plotting
fig = plt.figure(figsize=(18, 6))

# Plot for T
ax1 = fig.add_subplot(131, projection='3d')
ax1.scatter(df['OF'], df['Pc'], df['T'], c='r', marker='o', label='T Data')
ax1.set_xlabel('OF')
ax1.set_ylabel('Pc')
ax1.set_zlabel('T')
ax1.set_title('T')

# Plot for k
ax2 = fig.add_subplot(132, projection='3d')
ax2.scatter(df['OF'], df['Pc'], df['k'], c='g', marker='^', label='k Data')
ax2.set_xlabel('OF')
ax2.set_ylabel('Pc')
ax2.set_zlabel('k')
ax2.set_title('k')

# Plot for M
ax3 = fig.add_subplot(133, projection='3d')
ax3.scatter(df['OF'], df['Pc'], df['M'], c='b', marker='s', label='M Data')
ax3.set_xlabel('OF')
ax3.set_ylabel('Pc')
ax3.set_zlabel('M')
ax3.set_title('M')

plt.tight_layout()
plt.show()



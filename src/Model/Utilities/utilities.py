import pandas as pd
import matplotlib.pyplot as plt

def units(value, conversion, n=1):
    conversionDict = {
        "m_to_in": 39.37,
        "in_to_m": 1/39.37,

        "m2_to_in2": 1550,
        "in2_to_m2": 1/1550,

        "m3_to_in3": 61020,
        "in3_to_m3": 1/61020,

        "psi_to_pa": 6895,
        "pa_to_psi": 1/6895,

        "lbf_to_N": 4.448,
        "N_to_lbf": 1/4.448,

        "m3_to_L": 1000,
        "L_to_m3": 1/1000,

        "ai_to_am": (0.0254**(1 + 2*(n)))*(0.453592**(-n))
    }
    return value * conversionDict[conversion]

def plot_as_individual(df):
    # Thrust Profile
    plt.figure(figsize=(5, 4))
    plt.grid(True)
    plt.plot(df['time'], df['thrust'], linewidth=2)
    plt.title('Thrust vs. Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Thrust (lbf)')
    plt.show()

    # Oxidizer Tank and Combustion Chamber Pressure
    plt.figure(figsize=(5, 4))
    plt.grid(True)
    plt.plot(df['time'], df['p_chmb'], linewidth=2, label='Chamber Pressure')
    plt.plot(df['time'], df['p_tank'], linewidth=2, label='Tank Pressure')
    plt.title('Oxidizer Tank and Combustion Chamber Pressure vs. Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure (psi)')
    plt.legend()
    plt.show()

    # Oxidizer Tank Equilibrium of Nitrous Oxide Mass
    plt.figure(figsize=(5, 4))
    plt.grid(True)
    plt.plot(df['time'], df['n_go'], 'b', linewidth=2, label='Mass of N2O gas')
    plt.plot(df['time'], df['n_lo'], 'g', linewidth=2, label='Mass of N2O liquid')
    plt.title('Mass of N2O vs. Time')
    plt.xlabel('Time [s]')
    plt.ylabel('Mass of N2O [Moles]')
    plt.legend()
    plt.show()

    # OF Profile
    plt.figure(figsize=(5, 4))
    plt.grid(True)
    plt.plot(df['time'], df['OF'], linewidth=2)
    plt.title('Oxidizer to Fuel Ratio vs. Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Oxidizer to Fuel Ratio (O/F)')
    plt.show()

    # Fuel Grain Radius Profile
    plt.figure(figsize=(5, 4))
    plt.grid(True)
    plt.plot(df['time'], df['r'], linewidth=2)
    plt.title('Fuel Grain Burn over Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Fuel Grain Burn Rate (r) [in]')
    plt.show()

    # Combustion Chamber Profile
    plt.figure(figsize=(5, 4))
    plt.grid(True)
    plt.plot(df['time'], df['V_chmb'], linewidth=2)
    plt.title('Combustion Chamber Volume vs. Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Combustion Chamber Volume [in^3]')
    plt.show()

def plot_as_subplots(df):
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))  # 2 rows, 3 columns
    fig.suptitle('Rocket Analysis Plots')  # Overall title

    # thrust Profile
    axs[0, 0].grid(True)
    axs[0, 0].plot(df['time'], df['thrust'], linewidth=2)
    axs[0, 0].set_title('thrust vs. time')
    axs[0, 0].set_xlabel('time (s)')
    axs[0, 0].set_ylabel('thrust (lbf)')

    # Oxidizer Tank and Combustion Chamber Pressure
    axs[0, 1].grid(True)
    axs[0, 1].plot(df['time'], df['p_chmb'], linewidth=2, label='Chamber Pressure')
    axs[0, 1].plot(df['time'], df['p_tank'], linewidth=2, label='Tank Pressure')
    axs[0, 1].set_title('Oxidizer Tank and Combustion Chamber Pressure vs. time')
    axs[0, 1].set_xlabel('time (s)')
    axs[0, 1].set_ylabel('Pressure (psi)')
    axs[0, 1].legend()

    # Oxidizer Tank Equilibrium of Nitrous Oxide Mass
    axs[0, 2].grid(True)
    axs[0, 2].plot(df['time'], df['n_go'], 'b', linewidth=2, label='Mass of N2O gas')
    axs[0, 2].plot(df['time'], df['n_lo'], 'g', linewidth=2, label='Mass of N2O liquid')
    axs[0, 2].set_title('Mass of N2O vs. time')
    axs[0, 2].set_xlabel('time [s]')
    axs[0, 2].set_ylabel('Mass of N2O [Moles]')
    axs[0, 2].legend()

    # OF Profile
    axs[1, 0].grid(True)
    axs[1, 0].plot(df['time'], df['OF'], linewidth=2)
    axs[1, 0].set_title('Oxidizer to Fuel Ratio vs. time')
    axs[1, 0].set_xlabel('time (s)')
    axs[1, 0].set_ylabel('Oxidizer to Fuel Ratio (O/F)')

    # Fuel Grain Radius Profile
    axs[1, 1].grid(True)
    axs[1, 1].plot(df['time'], df['r'], linewidth=2)
    axs[1, 1].set_title('Fuel Grain Burn over time')
    axs[1, 1].set_xlabel('time (s)')
    axs[1, 1].set_ylabel('Fuel Grain Burn Rate (r) [in]')

    # Combustion Chamber Profile
    axs[1, 2].grid(True)
    axs[1, 2].plot(df['time'], df['V_chmb'], linewidth=2)
    axs[1, 2].set_title('Combustion Chamber Volume vs. time')
    axs[1, 2].set_xlabel('time (s)')
    axs[1, 2].set_ylabel('Combustion Chamber Volume [in^3]')

    plt.tight_layout()
    plt.subplots_adjust(top=0.90)  # To make space for the overall title
    plt.show()
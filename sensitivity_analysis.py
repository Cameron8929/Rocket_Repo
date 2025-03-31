import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

def HEM(c_d, diameter, num_orf, temp, exit_pressures):
    """Calculate N2O mass flow rate using the Homogeneous Equilibrium Model (HEM).

    Args:
        c_d (float): Discharge coefficient.
        diameter (float): Orifice diameter (metres).
        num_orf (int): Number of orifices.
        temp (float): Temperature (Kelvin).
        exit_pressures (array): Exit pressure values (Pascals).

    Returns:
        array: Mass flow rates (kg/s) for each exit pressure.
    """
    P_i = CP.PropsSI('P', 'T', temp, 'Q', 0, 'N2O')
    h_i = CP.PropsSI("H", "T", temp, "Q", 0, "N2O")
    s_i = CP.PropsSI("S", "T", temp, "Q", 0, "N2O")
    mass_flow_rate = []
    max_val = -10000
    for i in exit_pressures:
        h_o = CP.PropsSI("H", "P", i, "S", s_i, "N2O")
        rho_o = CP.PropsSI("D", "P", i, "S", s_i, "N2O")
        mass_flow_curr = c_d * rho_o * np.sqrt(2 * (h_i - h_o)) * np.pi * ((diameter/2) ** 2)
        mass_flux = c_d * rho_o * np.sqrt(2 * (h_i - h_o))
        if mass_flow_curr > max_val:
            max_val = mass_flux
            mass_flow_rate.append(mass_flux)
        else:
            mass_flow_rate.append(mass_flux)
    return np.array(mass_flow_rate)

def SPI(c_d, diameter, num_orf, temp, exit_pressure):
    """Calculate N2O mass flow rate using the Single-Phase Incompressible (SPI) model.

    Args:
        c_d (float): Discharge coefficient.
        diameter (float): Orifice diameter (metres).
        num_orf (int): Number of orifices.
        temp (float): Temperature (Kelvin).
        exit_pressure (array): Exit pressure values (Pascals).

    Returns:
        array: Mass flow rates (kg/s) for each exit pressure.
    """
    P_i = CP.PropsSI('P', 'T', temp, 'Q', 0, 'N2O')
    rho_l = CP.PropsSI('D', 'T', temp, 'Q', 0, 'N2O')
    mass_flow_rate = []
    for i in exit_pressure:
        mass_flow_curr = c_d * np.sqrt(2 * rho_l * (P_i - i)) 
        mass_flow_rate.append(mass_flow_curr)

    return np.array(mass_flow_rate)

def DYER(c_d, diameter, num_orf, temp, exit_pressure):
    """Calculate N2O mass flow rate using the Dyer model (average of HEM and SPI).

    Args:
        c_d (float): Discharge coefficient.
        diameter (float): Orifice diameter (metres).
        num_orf (int): Number of orifices.
        temp (float): Temperature (Kelvin).
        exit_pressure (array): Exit pressure values (Pascals).

    Returns:
        list: Mass flow rates (kg/s) for each exit pressure.
    """
    area = num_orf * np.pi * ((diameter/2) ** 2)
    mdot_spi = SPI(c_d, diameter, num_orf, temp, exit_pressure)
    mdot_hem = HEM(c_d, diameter, num_orf, temp, exit_pressure)
    mass_flow_rate = []
    for i, j, k in zip(exit_pressure, mdot_spi, mdot_hem):
        # Because at saturated
        mdot_dyer = ((1/2) * j * area) + ((1/2) * k * area)
        mass_flow_rate.append(mdot_dyer)
    
    return mass_flow_rate

def DYER_validate(c_d, diameter, num_orf):
    """Validate the Dyer model with a range of exit pressures.

    Args:
        c_d (float): Discharge coefficient.
        diameter (float): Orifice diameter (metres).
        num_orf (int): Number of orifices.

    Returns:
        array: Predicted mass flow rates (kg/s) for validation.
    """
    predicted_mass_flow = DYER(c_d, diameter, num_orf, 280, np.linspace(100000, 3e6, 100))
    return predicted_mass_flow

def sensitivity_analysis():
    """Perform sensitivity analysis of the Dyer model by varying discharge coefficient.

    Plots:
        Graph of mass flow rate (kg/s) vs. pressure drop (MPa) for different C_d values at 280 K.
    """
    delta_P_MPa = np.linspace(0.1, 3.0, 100)
    
    # Calculate saturation pressure at 280K
    P_i = CP.PropsSI('P', 'T', 280, 'Q', 0, 'N2O')
    
    # Calculate exit pressures by subtracting ΔP from saturation pressure
    exit_pressures_Pa = P_i - (delta_P_MPa * 1e6)
    
    # Rest of the function remains the same...
    m_1 = DYER(0.567, 0.0015, 1, 280, exit_pressures_Pa)
    m_2 = DYER(0.630, 0.0015, 1, 280, exit_pressures_Pa)
    m_3 = DYER(0.693, 0.0015, 1, 280, exit_pressures_Pa)
    
    # Plot directly with delta_P_MPa
    plt.figure(figsize=(10, 6))
    plt.plot(delta_P_MPa, m_1, label=r'$C_d$ = 0.567')
    plt.plot(delta_P_MPa, m_2, label=r'$C_d$ = 0.630')
    plt.plot(delta_P_MPa, m_3, label=r'$C_d$ = 0.693')
    plt.xlabel(r'$\Delta P$ [MPa]')
    plt.ylabel(r'$\dot{m}$ [kg/s]')
    plt.title(r'Impact of Discharge Coefficient ($C_d$) on N$_2$O Mass Flow Rate')
    plt.legend()
    plt.grid()
    plt.savefig("plots/sense_analysis.pdf")
    plt.show()

def plot_dyer_standard():
    """Plot the Dyer model mass flow rate with standard parameters.

    Plots:
        Graph of mass flow rate (kg/s) vs. pressure drop (MPa) at 280 K with C_d=0.630.
    """
    delta_P_MPa = np.linspace(0.1, 3.0, 100)
    P_i = CP.PropsSI('P', 'T', 280, 'Q', 0, 'N2O')
    exit_pressures_Pa = P_i - (delta_P_MPa * 1e6)
    m_2 = DYER(0.630, 0.0015, 1, 280, exit_pressures_Pa)
    # Plot directly with delta_P_MPa
    plt.figure(figsize=(10, 6))
    plt.plot(delta_P_MPa, m_2, label="Dyer")
    plt.xlabel(r'$\Delta P$ [MPa]')
    plt.ylabel(r'$\dot{m}$ [kg/s]')
    plt.title(r'Impact of Discharge Coefficient ($C_d$) on N$_2$O Mass Flow Rate')
    plt.legend()
    plt.grid()
    plt.savefig("plots/sense_analysis.pdf")
    plt.show()

def plot_dyer_with_user_input(cd, dia, num_orf, temp):
    """Plot the Dyer model mass flow rate with user-defined parameters.

    Args:
        cd (float): Discharge coefficient.
        dia (float): Orifice diameter (metres).
        num_orf (int): Number of orifices.
        temp (float): Temperature (Kelvin).

    Plots:
        Graph of mass flow rate (kg/s) vs. pressure drop (MPa) for specified conditions.
    """
    delta_P_MPa = np.linspace(0.1, 3.0, 100)
    P_i = CP.PropsSI('P', 'T', temp, 'Q', 0, 'N2O')
    exit_pressures_Pa = P_i - (delta_P_MPa * 1e6)
    m_2 = DYER(cd, dia, num_orf, temp, exit_pressures_Pa)
    # Plot directly with delta_P_MPa
    plt.figure(figsize=(10, 6))
    plt.plot(delta_P_MPa, m_2, label=rf"$C_d={cd}$")
    plt.xlabel(r'$\Delta P$ [MPa]')
    plt.ylabel(r'$\dot{m}$ [kg/s]')
    plt.title(rf'Test Conditions: $C_d={cd}$, Diameter = {dia} m, Number of Orifices = {num_orf}, Temperature = {temp} K')
    plt.legend()
    plt.grid()
    plt.savefig("plots/user_input.pdf")
    plt.show()

def main_sensitive():
    """Run the N2O flow model with user-selected options.

    Options:
        1. Sensitivity analysis with default parameters.
        2. Custom parameters for the Dyer model.
    """
    print("\033[1m------------------N2O Flow Model - Options:------------------\033[0m")
    print("1. Run sensitivity analysis with default values (Cd=0.630, d=0.0015m, 1 orifice, T=280K)")
    print("2. Enter custom parameters for DYER model")
    
    choice = input("Enter your choice (1 or 2): ")
    
    if choice == "1":
        sensitivity_analysis()
    elif choice == "2":
        try:
            cd = float(input("Enter discharge coefficient (Cd): "))
            diameter = float(input("Enter orifice diameter (in meters): "))
            num_orifices = int(input("Enter number of orifices: "))
            temperature = float(input("Enter temperature (in Kelvin): "))
            
            plot_dyer_with_user_input(cd, diameter, num_orifices, temperature)
        except ValueError:
            print("Error: Please enter valid numerical values.")
            main_sensitive()
    else:
        print("Invalid choice. Please enter 1 or 2.")
        main_sensitive()
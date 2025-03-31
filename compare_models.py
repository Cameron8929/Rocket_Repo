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
        delta_p_values (array): Pressure drop values (Pascals).

    Returns:
        list: Mass flow rates (kg/s) for each pressure drop.
    """
    P_i = CP.PropsSI('P', 'T', temp, 'Q', 0, 'N2O')
    h_i = CP.PropsSI("H", "T", temp, "Q", 0, "N2O")
    s_i = CP.PropsSI("S", "T", temp, "Q", 0, "N2O")
    mass_flow_rate = []
    max_val = -10000
    reached_max = False
    
    for i in exit_pressures:
        if i <= 0:
            mass_flow_rate.append(max_val if max_val > 0 else 0)
            continue
            
        try:
            h_o = CP.PropsSI("H", "P", i, "S", s_i, "N2O")
            rho_o = CP.PropsSI("D", "P", i, "S", s_i, "N2O")
            
            mass_flux = c_d * rho_o * np.sqrt(2 * (h_i - h_o))
            
            if reached_max:
                mass_flow_rate.append(max_val)
            elif mass_flux > max_val:
                max_val = mass_flux
                mass_flow_rate.append(mass_flux)
            else:
                reached_max = True
                mass_flow_rate.append(max_val)
        except:
            mass_flow_rate.append(max_val if max_val > 0 else 0)
            
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

def plot_all_models_comparison():
    """Plot a comparison of Dyer, HEM, and SPI models for N2O mass flow rate.

    Plots:
        Graph of mass flow rate (kg/s) vs. pressure drop (bar) at 280 K.
    """
    # Set parameters
    cd = 0.63
    diameter = 0.0015
    num_orf = 1
    temp = 280
    
    # Convert bar to Pa for pressure range (0.1 to 30 bar)
    delta_P_bar = np.linspace(0.1, 30, 100)
    delta_P_Pa = delta_P_bar * 1e5
    
    # Get saturation pressure at specified temperature
    P_i = CP.PropsSI('P', 'T', temp, 'Q', 0, 'N2O')
    
    # Calculate exit pressures based on pressure drops
    exit_pressures_Pa = P_i - delta_P_Pa
    
    # Apply minimum pressure check
    exit_pressures_Pa = np.maximum(exit_pressures_Pa, 100)  # Avoid negative or very low pressures
    
    # Calculate flow rates for each model
    area = num_orf * np.pi * ((diameter/2) ** 2)
    mdot_dyer = DYER(cd, diameter, num_orf, temp, exit_pressures_Pa)
    mdot_hem_flux = HEM(cd, diameter, num_orf, temp, exit_pressures_Pa)
    mdot_spi_flux = SPI(cd, diameter, num_orf, temp, exit_pressures_Pa)
    
    # Convert flux to flow rate for HEM and SPI
    mdot_hem = mdot_hem_flux * area
    mdot_spi = mdot_spi_flux * area
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(delta_P_bar, mdot_dyer, label="DYER Model", linewidth=2)
    plt.plot(delta_P_bar, mdot_hem, label="HEM Model", linewidth=2)
    plt.plot(delta_P_bar, mdot_spi, label="SPI Model", linewidth=2)
    
    plt.xlabel(r'$\Delta P$ [bar]')
    plt.ylabel(r'$\dot{m}$ [kg/s]')
    plt.title(r'Comparison of Flow Models for N$_2$O at Saturation Conditions')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig("plots/model_comparison.pdf")
    plt.show()


import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

def HEM_with_delta_p_given(c_d, diameter, num_orf, temp, delta_p_values):
    """Calculate N2O mass flow rate using the Homogeneous Equilibrium Model (HEM).

    Args:
        c_d (float): Discharge coefficient.
        diameter (float): Orifice diameter (metres).
        num_orf (int): Number of orifices.
        temp (float): Temperature (Kelvin).
        delta_p_values (array): Pressure drop values (Pascals).

    Returns:
        array: Mass flow rates (kg/s) for each pressure drop.
    """
    P_i = CP.PropsSI('P', 'T', temp, 'Q', 0, 'N2O')
    h_i = CP.PropsSI("H", "T", temp, "Q", 0, "N2O")
    s_i = CP.PropsSI("S", "T", temp, "Q", 0, "N2O")
    
    mass_flow_rate = []
    max_val = -10000
    for delta_p in delta_p_values:
        P_o = P_i - delta_p
        h_o = CP.PropsSI("H", "P", P_o, "S", s_i, "N2O")
        rho_o = CP.PropsSI("D", "P", P_o, "S", s_i, "N2O")

        mass_flow_curr = c_d * rho_o * np.sqrt(2 * (h_i - h_o)) * np.pi * ((diameter/2) ** 2)
        mass_flux = c_d * rho_o * np.sqrt(2 * (h_i - h_o))

        if mass_flow_curr > max_val:
            max_val = mass_flux
            mass_flow_rate.append(mass_flux)
        else:
            mass_flow_rate.append(mass_flux)
        
    return np.array(mass_flow_rate)

def SPI_with_delta_p_given(c_d, diameter, num_orf, temp, exit_pressure):
    """Calculate N2O mass flow rate using the Single-Phase Incompressible (SPI) model.

    Args:
        c_d (float): Discharge coefficient.
        diameter (float): Orifice diameter (metres).
        num_orf (int): Number of orifices.
        temp (float): Temperature (Kelvin).
        exit_pressure (array): Pressure drop values (Pascals).

    Returns:
        list: Mass flow rates (kg/s) for each pressure drop.
    """
    rho_l = CP.PropsSI('D', 'T', temp, 'Q', 0, 'N2O')
    mass_flow_rate = []
    for i in exit_pressure:
        mass_flow_curr = c_d * np.sqrt(2 * rho_l * (i)) 
        mass_flow_rate.append(mass_flow_curr)

    return mass_flow_rate

def DYER_with_delta_p_given(c_d, diameter, num_orf, temp, exit_pressure):
    """Calculate N2O mass flow rate using the Dyer model (average of HEM and SPI).

    Args:
        c_d (float): Discharge coefficient.
        diameter (float): Orifice diameter (metres).
        num_orf (int): Number of orifices.
        temp (float): Temperature (Kelvin).
        exit_pressure (array): Pressure drop values (Pascals).

    Returns:
        list: Mass flow rates (kg/s) for each pressure drop.
    """
    area = num_orf * np.pi * ((diameter/2) ** 2)
    mdot_spi = SPI_with_delta_p_given(c_d, diameter, num_orf, temp, exit_pressure)
    mdot_hem = HEM_with_delta_p_given(c_d, diameter, num_orf, temp, exit_pressure)
    mass_flow_rate = []
    for i, j, k in zip(exit_pressure, mdot_spi, mdot_hem):
        # Because at saturated
        mdot_dyer = ((1/2) * j * area) + ((1/2) * k * area)
        mass_flow_rate.append(mdot_dyer)
    
    return mass_flow_rate

def calculate_mape_value(experimental, predicted):
    """Calculate the Mean Absolute Percentage Error (MAPE) between experimental and predicted values.

    Args:
        experimental (array): Experimental mass flow rates (kg/s).
        predicted (array): Predicted mass flow rates (kg/s).

    Returns:
        float: MAPE value as a percentage.
    """
    mape_total = 0
    for exp, pred in zip(experimental, predicted):
        curr = np.abs((exp-pred)/exp)
        mape_total += curr
    return 1/(len(experimental)) * 100 * mape_total

def DYER_validate(c_d, diameter, num_orf):
    """Validate the Dyer model by comparing predicted and experimental N2O mass flow rates.

    Args:
        c_d (float): Discharge coefficient.
        diameter (float): Orifice diameter (metres).
        num_orf (int): Number of orifices.

    Plots:
        Graph of mass flow rate (kg/s) vs. pressure drop (MPa) at 0°C and P_i = 3.12 MPa.
    """
    plt.figure(figsize=(10, 8))
    delta_p, mass_flow_rate_experimental = read_data("DYER_DATA/dyer.txt")
    predicted_mass_flow = DYER_with_delta_p_given(c_d, diameter, num_orf, 273, delta_p)
    print("MAPE for 3.12 MPa: " + f"{calculate_mape_value(mass_flow_rate_experimental, predicted_mass_flow):.2f}%")
    plt.plot(delta_p/1e6, mass_flow_rate_experimental, label = "Experimental P_i = 3.12 MPa")
    plt.plot(delta_p/1e6, predicted_mass_flow, label = "Predicted P_i = 3.12 MPa")
    plt.xlabel(r'$\Delta P$ [MPa]')
    plt.ylabel(r'$\dot{m}$ [kg/s]')
    plt.title(r'DYER Model: N₂O Mass Flow Rate vs. Pressure Drop (T = 0°C, $P_i$ = 3.12 MPa)')
    plt.legend()
    plt.grid()
    plt.savefig("plots/dyer_model.pdf")
    plt.show()

def read_data(file_name):
    """Read pressure drop and mass flow rate data from a text file.

    Args:
        file_name (str): Path to the text file containing data.

    Returns:
        tuple: Arrays of pressure drops (Pascals) and mass flow rates (kg/s).
    """
    lines = []
    with open(file_name, mode ='r') as file:    
        lines = file.readlines()
    
    pressure_drop = []
    mass_flow_rate = []
    i = 0
    while i < len(lines):
        l_current = lines[i].strip()
        l_current = l_current.split(", ")
        pressure_drop.append(float(l_current[0]))
        mass_flow_rate.append(float(l_current[1]))
        i += 1
    pressure_drop = np.array(pressure_drop)
    mass_flow_rate = np.array(mass_flow_rate)
    pressure_drop = pressure_drop * 1e6
    mass_flow_rate = mass_flow_rate 
    return pressure_drop, mass_flow_rate

def dyer_validate_main():
    """Run Dyer model validation with default parameters.

    Calls DYER_validate with c_d=0.75, diameter=0.0015 m, and num_orf=1.
    """
    DYER_validate(0.75, 0.0015, 1)
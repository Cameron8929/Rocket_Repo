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
        list: Mass flow rates (kg/s) for each pressure drop.
    """
    P_i = CP.PropsSI('P', 'T', temp, 'Q', 0, 'N2O')
    h_i = CP.PropsSI("H", "T", temp, "Q", 0, "N2O")
    s_i = CP.PropsSI("S", "T", temp, "Q", 0, "N2O")
    
    mass_flow_rate = []
    # Negative large value to be immediately changed to
    max_val = -10000
    for delta_p in delta_p_values:
        P_o = P_i - delta_p
        h_o = CP.PropsSI("H", "P", P_o, "S", s_i, "N2O")
        rho_o = CP.PropsSI("D", "P", P_o, "S", s_i, "N2O")
        area = num_orf * np.pi * ((diameter/2) ** 2)
        mass_flow_curr = c_d * area * rho_o * np.sqrt(2 * (h_i - h_o))
        if mass_flow_curr > max_val:
            max_val = mass_flow_curr
            mass_flow_rate.append(mass_flow_curr)
        else:
            mass_flow_rate.append(max_val)
        
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


def hem_validate(c_d, diameter, num_orf):
    """Validate the HEM model by comparing predicted and experimental N2O mass flow rates.
    Data used can be found in HEM_DATA/
    Args:
        c_d (float): Discharge coefficient.
        diameter (float): Orifice diameter (metres).
        num_orf (int): Number of orifices.

    Plots:
        Graph of mass flow rate (kg/s) vs. pressure drop (MPa) at 0°C and P_i = 3.12 MPa.
    """
    delta_p_hem, mass_flow_rate_hem = read_data("HEM_DATA/hem.txt")
    predicted_mass_flow_rate = HEM_with_delta_p_given(c_d, diameter, num_orf, 273, delta_p_hem)
    print("MAPE for P_i = 3.12 MPa: " + f"{calculate_mape_value(mass_flow_rate_hem, predicted_mass_flow_rate):.2f}%")
    
    plt.figure(figsize=(10, 6))
    
    if delta_p_hem[0] > 0:
        delta_p_with_zero = np.concatenate(([0], delta_p_hem))
        exp_with_zero = np.concatenate(([0], mass_flow_rate_hem))
        pred_with_zero = np.concatenate(([0], predicted_mass_flow_rate))
    else:
        delta_p_with_zero = delta_p_hem
        exp_with_zero = mass_flow_rate_hem
        pred_with_zero = predicted_mass_flow_rate
    
    plt.plot(delta_p_with_zero/1e6, exp_with_zero, color='blue', label="Experimental")
    plt.plot(delta_p_with_zero/1e6, pred_with_zero, color='red', label="HEM Model")
    
    plt.xlim(0, max(delta_p_hem/1e6))
    plt.ylim(0, max(max(mass_flow_rate_hem), max(predicted_mass_flow_rate)) * 1.1)
    
    plt.xlabel(r'$\Delta P$ [MPa]')
    plt.ylabel(r'$\dot{m}$ [kg/s]')
    plt.title(r'HEM Model: N₂O Mass Flow Rate vs. Pressure Drop (T = 0°C, $P_i$ = 3.12 MPa)')
    plt.legend()
    plt.grid(True)
    plt.savefig("plots/hem_model.pdf")
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

def calculate_mape_value(experimental, predicted):
    """Calculates the mape value

    Args:
        experimental (list): Experimental values
        predicted (list): Predicted values by the HEM model.

    Returns:
        float: Returns the mean absolute percentage error.
    """
    mape_total = 0
    for exp, pred in zip(experimental, predicted):
        curr = np.abs((exp-pred)/exp)
        mape_total += curr
    return 1/(len(experimental)) * 100 * mape_total

def HEM_validate_main():
    hem_validate(0.75, 0.0015, 1)
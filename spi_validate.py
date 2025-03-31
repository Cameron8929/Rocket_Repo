import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

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
        mass_flow_curr = c_d * np.sqrt(2 * rho_l * (i)) * num_orf * np.pi * ((diameter/2) ** 2)
        mass_flow_rate.append(mass_flow_curr)

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

def SPI_validate(c_d, diameter, num_orf):
    """Validate the SPI model by comparing predicted and experimental N2O mass flow rates.

    Args:
        c_d (float): Discharge coefficient.
        diameter (float): Orifice diameter (metres).
        num_orf (int): Number of orifices.

    Plots:
        Graph of mass flow rate (kg/s) vs. pressure drop (MPa) for five inlet pressures.
    """
    # 600 PSI
    delta_p_600_psi, mass_flow_rate_600_psi_experimental = read_data("SPI_DATA/600_psi.txt")
    predicted_600_mass_flow_rate = SPI_with_delta_p_given(c_d, diameter, num_orf, 284, delta_p_600_psi)
    print("MAPE for 284 K, P_i = 4.14 MPa: " + f"{calculate_mape_value(mass_flow_rate_600_psi_experimental, predicted_600_mass_flow_rate):.2f}%")
    plt.plot(delta_p_600_psi/1e6, mass_flow_rate_600_psi_experimental, label = "Experimental 1")
    plt.plot(delta_p_600_psi/1e6, predicted_600_mass_flow_rate, label = "Predicted 1")
    
    # 700 PSI
    delta_p_700_psi, mass_flow_rate_700_psi_experimental = read_data("SPI_DATA/700_psi.txt")
    predicted_700_mass_flow_rate = SPI_with_delta_p_given(c_d, diameter, num_orf, 291, delta_p_700_psi)
    print("MAPE for 291 K, P_i = 4.83 MPa: " + f"{calculate_mape_value(mass_flow_rate_700_psi_experimental, predicted_700_mass_flow_rate):.2f}%")
    plt.plot(delta_p_700_psi/1e6, mass_flow_rate_700_psi_experimental, label = "Experimental 2")
    plt.plot(delta_p_700_psi/1e6, predicted_700_mass_flow_rate, label = "Predicted 2")

    # 800 PSI
    delta_p_800_psi, mass_flow_rate_800_psi_experimental = read_data("SPI_DATA/800_psi.txt")
    predicted_800_mass_flow_rate = SPI_with_delta_p_given(c_d, diameter, num_orf, 297, delta_p_800_psi)
    print("MAPE for 297 K, P_i = 5.51 MPa: " + f"{calculate_mape_value(mass_flow_rate_800_psi_experimental, predicted_800_mass_flow_rate):.2f}%")
    plt.plot(delta_p_800_psi/1e6, mass_flow_rate_800_psi_experimental, label = "Experimental 3")
    plt.plot(delta_p_800_psi/1e6, predicted_800_mass_flow_rate, label = "Predicted 3")

    # 900 PSI
    delta_p_900_psi, mass_flow_rate_900_psi_experimental = read_data("SPI_DATA/900_psi.txt")
    predicted_900_mass_flow_rate = SPI_with_delta_p_given(c_d, diameter, num_orf, 302, delta_p_900_psi)
    print("MAPE for 302 K, P_i = 6.2 MPa: " + f"{calculate_mape_value(mass_flow_rate_900_psi_experimental, predicted_900_mass_flow_rate):.2f}%")
    plt.plot(delta_p_900_psi/1e6, mass_flow_rate_900_psi_experimental, label = "Experimental 4")
    plt.plot(delta_p_900_psi/1e6, predicted_900_mass_flow_rate, label = "Predicted 4")

    # 1000 PSI
    delta_p_1000_psi, mass_flow_rate_1000_psi_experimental = read_data("SPI_DATA/1000_psi.txt")
    predicted_1000_mass_flow_rate = SPI_with_delta_p_given(c_d, diameter, num_orf, 307, delta_p_1000_psi)
    plt.plot(delta_p_1000_psi/1e6, mass_flow_rate_1000_psi_experimental, label = "Experimental 5")
    plt.plot(delta_p_1000_psi/1e6, predicted_1000_mass_flow_rate, label = "Predicted 5")
    print("MAPE for 307 K, P_i = 6.89 MPa: " + f"{calculate_mape_value(mass_flow_rate_1000_psi_experimental, predicted_1000_mass_flow_rate):.2f}%")
    plt.xlabel(r'$\Delta P$ [MPa]')
    plt.ylabel(r'$\dot{m}$ [kg/s]')
    plt.title(r'SPI Model: N₂O Mass Flow Rate vs. Pressure Drop for 5 Inlet Pressures')
    plt.grid()
    plt.legend()
    plt.savefig("plots/spi_model.pdf")
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

def SPI_validate_main():
    SPI_validate(0.75, 0.0015, 1)
# N₂O Flow Models

## Overview
This repository contains code to model and analyse N2O flow through orifices using three different approaches:
- Single-Phase Incompressible (SPI) model
- Homogeneous Equilibrium Model (HEM)
- Dyer model (a hybrid approach combining SPI and HEM)

Each model is validated against experimental data and compared for accuracy.

## Dependancies 
The program uses the following libraries: 
### matplotlib.pyplot 
Used for creating the validation plots for SPI
### CoolProp
Used for calculation of thermodynamic properties of nitrous oxide (N2O).
### numpy
Mathematical operations and array handling.

## Data Used for Validation and Plots Produced
Data was gathered for validation using automeris from graphs for saturated N2O from [1].
- `SPI_DATA/`: Experimental data for SPI model validation from Waxman Figure 2.2.
- `HEM_DATA/`: Experimental data for HEM model validation from Waxman Figure 2.20.
- `DYER_DATA/`: Experimental data for Dyer model validation from Waxman Figure 2.44.
- `plots/`: Generated plots and visualisations

## Models
All models are validated against experimental data. The Mean Absolute Percentage Error (MAPE) is calculated to quantify model accuracy for each test condition. All models use the following input conditions for validation - 1.5 mm, 0.75, 1 and 273 K for the diameter, number of orifices, discharge coefficient and ambient temperature respectively. 
### SPI Model
Models liquid flow assuming incompressible behaviour for saturated N2O.
`spi_validate.py`uses the same input parameters [1]. 

### HEM Model
Models two-phase flow with homogeneous equilibrium assumptions. Accounts for phase change effects.`hem_validate.py`uses the same input parameters [1].

### Dyer Model
A hybrid approach that combines SPI and HEM predictions (50% weight to each) because of saturation pressure. `dyer_validate.py`uses the same input parameters [1]. 

## Sensitivity Analysis
The mass flow rate vs. pressure drop is plotted in `sensitivity_analysis.py`. The user can perform the sensitivity analysis by typing '1' into the terminal or produce a plot of the mass flow rate vs. pressure drop using their own discharge coefficient, number of orifices, diameter (m) and ambient temperature (K). 

## Comparison of SPI, HEM and Dyer Model
The SPI, HEM and Dyer models are plotted against each other in `compare_models.py` to demonstrate how the models perform for the same input conditions. The conditions are - 1.5 mm, 0.63, 1 and 280 K for the diameter, number of orifices, discharge coefficient and ambient temperature respectively. 

## References
[1] B. S. Waxman, "An investigation of injectors for use with high vapor pressure propellants with applications to hybrid rockets," Ph.D. dissertation, Dept. Aeronautics Astronautics, Stanford Univ., Stanford, CA, USA, 2014. [Online]. Available: https://stacks.stanford.edu/file/druid:ng346xh6244/BenjaminWaxmanFinal-augmented.pdf
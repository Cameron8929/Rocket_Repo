import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
from dyer_validate import *
from hem_validate import *
from spi_validate import *
from sensitivity_analysis import *
from compare_models import *

def main():
    print("\n\033[1m------------------Validation------------------\033[0m")
    print("\n\033[1mSPI Validation\033[0m")
    SPI_validate_main()
    print("\n\033[1mHEM Validation\033[0m")
    HEM_validate_main()
    print("\n\033[1mDYER Validation\033[0m")
    dyer_validate_main()
    print("\n\033[1mSPI Sensitivity Analysis\033[0m")
    main_sensitive()
    print("\n\033[1mModel Comparison Plots\033[0m")
    plot_all_models_comparison()
    
if __name__ == "__main__":
    main()
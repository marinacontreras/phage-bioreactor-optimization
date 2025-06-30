# phage-bioreactor-optimization

This repository contains the full implementation of a dynamic modeling and optimization framework for bacteriophage production in batch bioreactors. The model is implemented in Python using [CasADi](https://web.casadi.org/) and includes calibration, sensitivity analysis, and control optimization using metaheuristics.

## Repository Contents

### Main Notebook
- `Código_completo_G3.ipynb`: Main Jupyter notebook containing all steps of the modeling workflow:
  - Model definition
  - Parameter estimation
  - Sensitivity analysis
  - Optimization using Grey Wolf Optimizer (GWO)
  - Monte Carlo simulations

### Input Data
- `Datos_modelo_10min.xlsx`, `Datos_modelo_15min.xlsx`, `Datos_modelo_30min.xlsx`: Synthetic experimental datasets for model calibration, with different sampling frequencies (10, 15, and 30 minutes). Each file includes:
  - `Tiempo`: Time (h)
  - `cfu desde OD`: Estimated bacterial concentration (CFU/mL)
  - `pfu`: Phage concentration (PFU/mL)

### Monte Carlo Results
- `montecarlo_fijo_ori.xlsx`: Monte Carlo results for the experimental protocol (fixed timing and dose).
- `montecarlo_fijo.xlsx`: Monte Carlo results for the optimized strategy (fixed optimal timing and dose).
- `resultados_montecarlo2.csv`: Monte Carlo results where injection timing and dose were re-optimized for each simulation.

### Post-analysis
- `post_regression_analysis.py`: Python script for post-simulation statistical analysis and visualization.

## Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/marinacontreras/phage-bioreactor-optimization.git
   cd phage-bioreactor-optimization
2. Open and run `Código_completo_G3.ipynb` in Jupyter Notebook or Google Colab.

3. Use the .xlsx files as input data, and analyze the Monte Carlo outputs using the .csv and .xlsx result files.

4. (Optional) Run `post_regression_analysis.py` to perform visual/statistical analysis of the simulation results.

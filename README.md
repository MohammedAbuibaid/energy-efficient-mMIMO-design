# Energy-Efficient Massive MIMO Design: Optimal Number of Antennas Ensuring Guaranteed Bit Rate

This repository contains the resources and codebase for the paper titled "Energy-Efficient Massive MIMO Design: Optimal Number of Antennas Ensuring Guaranteed Bit Rate," published in the 2022 IEEE Future Networks World Forum (FNWF).

## Abstract
Our study explores the challenges and solutions in achieving energy efficiency in massive multi-input multi-output (mMIMO) antenna systems while ensuring a guaranteed bit rate (GBR). We introduce a novel algorithm based on symmetric game theory to optimize the energy efficiency (EE) of mMIMO systems in a multi-cell network, taking into account GBR constraints. Additionally, we propose a data traffic model that aligns various user equipment capabilities and mobile data applications with corresponding GBR levels.

## Key Contributions
- **Algorithm Development**: Developed an algorithm using the symmetric-game best response strategy to balance the GBR requirements with the highest possible energy efficiency in mMIMO designs.
- **Data Traffic Model**: Introduced a model that translates user bit rate requirements into GBR levels, considering a range of user equipment capabilities and applications.
- **Simulation Results**: Provided simulation results demonstrating the effectiveness of the proposed algorithm in achieving GBR requirements with minimal deviation from optimal EE targets.

## How to Use
Instructions on how to set up the environment, run simulations, and interpret results.
- To re-draw the figures, please execute the ```generateFigures.m```
- To re-produce the simulation results, please execute ```Run_Simulation_OptimumEE_BS_GDR.m``` 


## Citation
If you use the resources from this repository, please cite:

```bibtex
@INPROCEEDINGS{abuibaid2022mmimo,
  author={Abuibaid, Mohammed and St-Hilaire, Marc and Aldirmaz-Colak, Sultan and Eid, Imad},
  booktitle={2022 IEEE Future Networks World Forum (FNWF)}, 
  title={Energy-Efficient Massive MIMO Design: Optimal Number of Antennas Ensuring Guaranteed Bit Rate}, 
  year={2022},
  pages={328-333},
  doi={10.1109/FNWF55208.2022.00064}
}
```

# Overview 
This code implements the formulation in the following paper. The explanation and formulation can be found in the paper. 

[1] K. Ho, M. Isaji, M. Patel, and K. Garoust, "Optimizing Flexible Complex Systems with Coupled and Co-Evolving Subsystems under Operational Uncertainties" Systems Engineering (accepted).

The code was written with MATLAB R2022b and to be used with Gurobi 11.0.1. 

# Code
## Main functions:
runall.m: run all codes

TwoSubsys_NLP_Sec5.m: Used in the Paper Sec. 5.1; Two Subsystem NLP case

NSubsys_NLP_partitioning_Sec6.m: Used in the Paper Sec. 6.2; including (1)Deterministic; (2)Fully Flexible; (3)Disconnecting weak coupling

WaterEnergyFoodNexus.m: Used in the Paper Sec. 5.3; Water Energy-Food-Nexus use case

ComputationalExperiments.m: Used in the Paper Sec. 5.2 and Sec. 6.2; test algorithms over random problem instances

## Subfunctions:

solve.m: Main solve function

ThreeSubsystemSub.m: Subfunction solving the three-subsystem case

gurobi_linprog.m: Solving LP using Gurobi (This subfunction was developed with the assistance of AI tools, but all of the code has been reviewed and adapted.)

## Data:

outputSec5.mat: Result from TwoSubsys_NLP_Sec5.m

outputSec6.mat: Result from NSubsys_NLP_partitioning_Sec6.m

outputWaterEnergyFoodNexus.mat: Result from WaterEnergyFoodNexus.m

computaionalResults_0.1.mat/computaionalResults_1.0.mat: Result from ComputationalExperiments.m

# Acknowledgment
This material is based upon work supported by the National Science Foundation under Award No. 1942559. Any opinions, findings and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation. Artificial intelligence (AI) technologies (e.g., Microsoft Copilot) were used to assist with code development. All AI-generated content was reviewed, validated, and refined by the authors, who take full responsibility for the final publication.

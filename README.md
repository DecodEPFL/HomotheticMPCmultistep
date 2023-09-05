# Homothetic tube Model Predictive Control with multi-step predictors

This GitHub repository contains the accompanying code for the paper titled **"Homothetic Tube Model Predictive Control with Multi-Step Predictors."** The code is implemented in MATLAB and includes two main files:

## Main Files:

1. **main_offlineComputation.m**: This script takes the identified parameters and the dataset contained in the file **"IdentifiedModelandData.mat"** as inputs. It is responsible for generating all the offline computed sets and parameters needed for the MPC (Model Predictive Control) approach described in the paper. Specifically, it computes the set **𝒳₀**, the feedback gain **K**, matrix **P**, and all other offline quantities defined in equations (10)-(12). These computed quantities are saved in a file named **"offlineComputation.mat"** for later use.

2. **main_HomoMPC_multistep.m**: This script serves as the main controller. It takes as input the data stored in **"offlineComputation.mat"** and simulates the system by applying the control law generated by solving the MPC problem addressed in the paper.

## Getting Started:

To use this code, follow these steps:

1. Clone this repository to your local machine using the following command: 
```
git clone https://github.com/DecodEPFL/HomotheticMPCmultistep.git
```

2. Ensure you have MATLAB installed on your system.

3. Open MATLAB and navigate to the cloned repository's directory.

4. Run **main_offlineComputation.m** to compute the offline quantities and save them to **"offlineComputation.mat"**.

5. After the offline computations are completed, you can run **main_HomoMPC_multistep.m** to simulate the system using the MPC controller.

## Note:

- Make sure to have the required datasets and identified parameters in the appropriate format as mentioned in the paper.

- The code provided here is a companion to the paper and should be used in conjunction with the paper's instructions and explanations for a comprehensive understanding of the algorithm and methodology.

- Feel free to reach out Danilo Saccani (danilo.saccani@epfl.ch) or Johannes Köhler (jkoehle@ethz.ch) if you have any questions.

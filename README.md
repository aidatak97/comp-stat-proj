#  Final Project for Computational Statistics SS2021


The main notebook is called Simulation_Project_SS21 and contains the simulation study **Evaluating Performance of Bayesian Classification Trees for High-Dimentional Prediction**. You can find the code in code.R and simulated data in results.RData.

For my final project in the Computational Statistics course I did a simulation study with application on microarray data. In my term paper, I focused on the classification task for a high-dimensional setting in order to evaluate the performance of probabilistically motivated methods based on ensembles of decision trees ( Linero, 2018). These methods were proposed to address the high dimensionality as according to the literature random forests might have some limitations when the number of predictors is much larger than the number of observations. Namely, I evaluated the performance of BART, DART and RF methods on the real data on patients (Gravier et al., 2010) and simulated data based on the Homogenous Gene Expression Model and Multi-Hit Model of Cancer (implemented in the R package "Umpire"). The results reinforced the controversial conclusions existing in the literature. In the simplified cases for the simulated DNA data, BART and DART yielded the best performance. However, in the complex data cases with large networks as well as in the case with real data RF outperformed both BART and DART.

The structure of the paper is as follows:

- Introduction
- Description of BART and DART
- Description of Homogenous Gene Expression Model and Multi-Hit Model of Cancer
- Simulation Results
- Empirical Application
- Conclusion

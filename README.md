<p align="center"><img width="100%" src="img/Trajectron++.png"/></p>

# ROMPC: Reduced Order Model Predictive Control #
This repository provides an implementation of the reduced order model predictive control scheme described in (journal paper)

## Relevant Publications ##
J. Lorenzetti and M. Pavone, [“Error Bounds for Reduced Order Model Predictive Control,”](https://arxiv.org/pdf/1911.12349.pdf) in _Proc. IEEE Conf. on Decision and Control_, Jeju Island, Republic of Korea, 2020. (Submitted)

J. Lorenzetti, B. Landry, S. Singh, and M. Pavone, [“Reduced Order Model Predictive Control For Setpoint Tracking,”](https://arxiv.org/pdf/1811.06590.pdf) in _European Control Conference_, Naples, Italy, 2019.

## Requirements ##
[MATLAB]([https://www.mathworks.com/products/matlab.html](https://www.mathworks.com/products/matlab.html))
[MPT3](https://www.mpt3.org/)
[YALMIP](https://yalmip.github.io/)
LP/QP solver (e.g [CPLEX](https://www.ibm.com/analytics/cplex-optimizer) or [MOSEK](https://www.mosek.com/))



## Example Models ##
The following examples are included in this repository:
| Model                                     | Description                                                                                                                                                                                          |
|-------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| smallSynthetic                                      | A synthetic discrete-time example with dimension $n^f = 6$.       |
| largeSynthetic                    | A synthetic discrete-time example with dimension $n^f = 20$.       |
| distillationColumn               | A time-discretized model of a binary distillation column [1],[2] with dimension $n^f = 86$. |
| tubularReactor | A continuous-time model of a controlled chemical reaction process [3] with dimension $n^f = 600$. |
| heatflow | A time-discretized model for a distributed control heat flow problem with dimension $n^f = 3,481$. This model is a modified version of the HF2D9 model described in [4].       |
| supersonicDiffuser | A time-discretized CFD model for the active control of a supersonic diffuser [5],[6] with dimension $n^f = 11,730$. |
| aircraft | A continuous-time aircraft dynamics model [7],[8] (with an embedded CFD aerodynamics model) with dimension $n^f = 0$.        |

## References ##
[1] S. Skogestad and M. Morari, "Understanding the dynamic behavior of distillation columns", Ind. & Eng. Chem. Research, 27, 10, 1848-1862 (1988) 

[2] S. Skogestad and I. Postlethwaite, "Multivariable Feedback Control", Wiley (1996) 

[3] O. Agudelo and J. Espinosa, "Control of a tubular chemical reactor by means of POD and predictive control techniques", European Control Conference (2007)

[4] F. Leibfritz, "COMPleib: Constrained matrix optimization problem library" (2006)

[5] K. Willcox and G. Lassaux, "Model Reduction of an Actively Controlled Supersonic Diffuser", Dimension Reduction of Large-Scale Systems, 357-361 (2005)

[6] G. Lassaux, "High-Fidelity Reduced-Order Aerodynamic Models: Application to Active Control of Engine Inlets", Master’s Thesis, Massachusetts Institute of Technology (2002)

[7] J. Lorenzetti, A. McClellan, C. Farhat, and M. Pavone, [“UAV Aircraft Carrier Landing Using CFD-Based Model Predictive Control,”](http://asl.stanford.edu/wp-content/papercite-data/pdf/Lorenzetti.McClellan.Farhat.Pavone.AIAA20.pdf) AIAA Scitech Forum, (2020)

[8] A. McClellan, J. Lorenzetti, M. Pavone, and C. Farhat, [“Projection-based Model Order Reduction for Flight Dynamics and Model Predictive Control,"](https://arc.aiaa.org/doi/abs/10.2514/6.2020-1190) AIAA Scitech Forum (2020)



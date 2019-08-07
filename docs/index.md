# On the perils of stabilizing prices when agents are learning

Replication instructions. 

<hr>

## Quick Start

1. Run ```setup.m```.
2. Run the file you are interested in, from the folder ```FiguresAndTables```.

That's it. 

<hr>

## The code

If you want to use this code for your work, here are some details about how it works. 

The code is divided in 5 folders:

* <b>```FiguresAndTables```</b>: includes files for replicating the graphs and tables in the paper

* <b>```Libraries```</b>: includes the latest version of ```CompEcon```, a set of additional computational routines (```libm```) and several utilities used in the code (```Utilities```)

* <b>```Results```</b>: where replication results are collected.

* <b>```SimulationsRoutines```</b>: this folder contains the files used to simulate the model, and the montecarlo routines. It also includes files for calculating the consumption equivalent welfare changes.

* <b>```SolveRoutines```</b>: contains the files for solving the model with collocation.

## Preliminary steps

The ```setup.m``` file adds the directories to the Matlab path. 

## Parameters for solving the model

We need to provide a set of parameters for the model, the collocation algorithm, and the benchmark parameters for the chosen functional forms.

First we define which type of gain we use. This is string that can be either decreasing (```'decr'```) or constant (```'const'```) gain:

```matlab
gaintype = 'const'; %'decr'; %
```

We then set the parameters of the model using the function ```benchmarkparameters``` using ```gaintype``` as input:

```matlab
p = benchmarkparameters(gaintype);
```

This function returns a structure ```p``` which contains all the parameters used in the model. [Here is the list of them](benchparam.md). If you want to overwrite one or more of these benchmark parameters, just add a line after the call of ```benchmarkparameters```. For example, if you want to change the value of the discount factor, you can write:

```matlab
p = benchmarkparameters(gaintype);
p.betta = 0.99;
```

Finally, we set up the parameters for the type of basis functions used in the collocation algorithm. This step uses CompEcon routines. The function ```projectionparameters``` takes two inputs. The first is the type of basis function, and the second (which is optional) is the order of the splines if we are using splines. If not set, the second argument is automatically assumed to be equal to ```[]``` (this implies that, in case of splines, the default is cubic). For example, for Chebichev polinomials we write:

```matlab
pp = projectionparameters('cheb');
```

For cubic splines, we can write:

```matlab
pp = projectionparameters('spli');
```

For splines of order 5, we use:

```matlab
pp = projectionparameters('spli', 5);
```

This function returns a sturcture which contains the two parameters used as inputs, called respectively ```pp.approxtype``` and ```pp.splineorder```.

The parameters' structures are taken as inputs in several files for solving and simulating the model.

## Solving the model

The function ```main_solver``` does all is needed.

```matlab
[solution, functionalspace, Grid, max_accuracy_test] = main_solver(benchmark_parameters,projection_parameters, gaintype);
```

It takes as inputs the two structures with parameters (for the model and for the collocation algorithm) and the gain type string.

It returns four results:

* ```solution```: this is a vector containing the coefficients of the basis functions that solve the problem

* ```functionalspace```: a structure containing the parameters defining the functional space for the basis functions (from CompEcon)

* ```Grid```: the grid over which the solution is computed

* ```max_accuracy_test```: accuracy test returning the maximum residuals error over a large grid

## Simulations

The simulation files are used for two purposes:

1. To create IRFs
2. To perform Montecarlo experiments for welfare analysis

For both tasks, we can use the same routines, by taking initial conditions and parameters conditional on the type of simulation we are performing.

### Simulation parameters

We need to setup the simulation parameters by using the function ```mcparameters```. The syntax is the following:

```matlab
sim_parameters = mcparameters(i,d,ps, seed0,x0,bg0, bp0,g0,sdl);
```

The arguments of this function are related to the simulation.

* ```i```:      number of draws
* ```d```:      number of repetitions of the draws (set to 1 for all purposes; this is legacy from previous versions of the code, will be eliminated in next version)
* ```ps```:     number of periods for the simulation
* ```seed```:   random number generator seed
* ```x0```:     output gap initial condition
* ```bg0```:    b_gap initial condition
* ```bp0```:    b_pi initial condition
* ```g0```:     gamma initial condition (essential for decreasing gain, could be anything when using constant gain since the simulation functions take the gain parameters from the benchmark parameters of the model)
* ```sdl```:    shutdownlearning (a parameter equal to 1 if we shut down learning, 0 otherwise; legacy for some experiments we did at the beginning of the project; set to 0 for all purposes)

The output is a structure with all parameters for the simulation/Montecarlo experiment, which can be easily passed to other functions for the simulations.

### Simulations/MC

Here is an example on how to run simulations:

```matlab
[gap_lag, pi , gap,b_pi,b_gap,...
    gamma_t,lambda1, welfare_cumul, welfare_x_cumul, ...
    welfare_pi_cumul, welfare_cumul_final, welfare_x_cumul_final , ...
    welfare_pi_cumul_final, costpushshock, welfare_inst, welfare_x_inst, ...
    welfare_pi_inst ] = simul(solution, functionalspace, benchmarkparameters, model, ...
    gaintype, shocktype, sim_parameters);
```

This is a bit involved, so we will decypher it in pieces. First, the inputs:

* ```solution```: vector of coefficients for the interpolation
* ```functionalspace```: functional space created with Miranda Fackler Compecon Toolbox
* ```benchmarkparameters```: structure containing benchmark parameters for the model
* ```model```: string, this is the model we are running, can be 'MMS', 'EHCOMM', 'EHDISCR', 'RECOMM', 'REDISCR'
* ```gaintype```: can be set to 'decr' for decreasing gain, or 'const' for constant gain
* ```shocktype```: string, can be set to 'irf' for impulse response function or 'series' for general Montecarlo simulation
* ```sim_parameters```: parameters for the simulation

The output is given by vectors (or matrices) containing the simulated series.

* ```gap_lag```: output gap lagged
* ```pi```: inflation
* ```gap```: output gap
* ```b_pi```: learning coefficient for inflation
* ```b_gap```: learning coefficient for output gap
* ```gamma_t```: learning parameter series (this will be constant if constant gain is chosen)
* ```lambda1```: Lagrange multiplier on the law of motion for learning coefficient on output gap
* ```welfare_cumul```: series of the welfare
* ```welfare_x_cumul```: series of the welfare from output gap
* ```welfare_pi_cumul```: series of the welfare from inflation
* ```welfare_cumul_final```: average series of the welfare
* ```welfare_x_cumul_final```: average series of the welfare from output gap
* ```welfare_pi_cumul_final```: average series of the welfare from inflation
* ```costpushshock```: the cost push shock
* ```welfare_inst```: instantanous undiscounted loss
* ```welfare_x_inst```: instantanous undiscounted loss from output gap
* ```welfare_pi_inst```: instantanous undiscounted loss from inflation

There are several functions that work in the same way (i.e. same syntax) for simulations:

* ```simul```
* ```simul_fast```: this is a faster version, which is mostly used for Montecarlo analysis.
<!-- * ```simul_fast_ALA```: this is used for the exercises related to <b>A</b>lternative <b>L</b>earning <b>A</b>lgorithms. -->

#### Montecarlo simulations for consumption equivalent welfare analysis

For this task, the function to use is ```mc_cons_equiv```:

```matlab
mcce = mc_cons_equiv(solution, functionalspace, benchmarkparameters, model, gaintype, shocktype, sim_parameters)
```

The output is the consumption equivalent welfare measure, calculated by calling the function ```consumption_equivalent``` after a Montecarlo simulation performed with ```simul_fast```.

<!-- Notice that there is also a version of this function for alternative learning algorithms (```mc_cons_equiv_ALA```). -->

## Reproducing graphs and tables in the paper

The folder ```FiguresAndTables``` contains the files for replicating the paper. The names of the files should be self-explanatory. 


## If there is any problem:

Open an issue in the [Github repository](https://github.com/meleantonio/perils-stab-prices), or contact Antonio Mele at meleantonio@gmail.com

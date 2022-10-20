# Cryopeg-Carbon-Model
A model of the carbon cycle and bacterial growth in cryopeg brines. Aimed at constraining the energetic kinetics of microbial populations in subzero hypersaline brines.

## Reference
Kannan et al. _In prep._

## Installation
1. Install [Julia](https://julialang.org/downloads/)
2. Install [Python](https://www.python.org/downloads/), guarranteed to compile on Python 3.9
3. Create a Python virtual environment [How-to](https://docs.python.org/3/library/venv.html)
4. Install the requirements in your venv `pip install -r requirements.txt`
5. Launch Julia in the terminal and do `using Pkg; Pkg.add(["DifferentialEquations", "DiffEqSensitivity", "ForwardDiff", "GlobalSensitivity"])
6. Run main.py to reproduce presented analyses.

## Documentation
I have strived to write clear, readable and well commented code. Docstrings exist for all files, large functions and classes. 
Comments throughout the code should guide your understanding. Especially in the julia code.
The best documentation is the paper cited above.

## Architecture & Codebase Design
The model is coded in Julia and is contained in `model.jl`. This file contains the code that defines and solves the differential equations of the model, and performs the sensitivyt analysis.

The python code interfaces with the model to perform different analyses. It first defines Scenario objects that contain paramaters to run all analyses. The Scenario class also defines common paramters that are not meant to be changed (designated by the standard leading _).

An analysis object is then created, it is assigned a scenario, and executes all analyses. Individual analyses can be performed by calling their respective functions.

The main file and `plots.py` both participate in figure generation and value logging. 

`utils.py` contains miscellaneous classes and functions.

I will not that since passing objects between Python and Julia is not nice, model paramters are passed as *ordered arrays*. If you decide to extend or modify this code you must *ensure this order is respected* across Julia and Python code. The class structures were designed to maximize encapsulation and minimize the number of places modification would have to be carried over.

## Issues & Contribution
Should you have any questions or suggestions, please open an issue on this repository. 

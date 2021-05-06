# Trajectory Optimization based on Subdivison

## Dataset
Here is environment meshes data and initial trajectories data in our paper. 

- [Data](https://drive.google.com/file/d/1DM86tO0wUNef2G3BqX1U6s52vXGT5wuf/view?usp=sharing/)


## Installation

#### via CMake

Our code was originally developed on Ubuntu 16.04. We provide the commands for installing our project in Ubuntu 16.04:

- Clone the repository into your local machine:

```bash
git clone https://github.com/ruiqini/traj-opt-subdivision.git
```

- Compile the code using cmake (default in Release mode):

```bash
cd traj-opt-subdivision
mkdir build
cd build
cmake ..
make
```

## Usage

#### Input data

Extract compressed data file, there are three folders inside: `Config file/` `init/` `mesh/`. 
`Config file/3D.json` saves input parameters in paper.
`init/name_init.txt` saves initial trajectory of input environment mesh `name`.
`mesh/name` saves input environment mesh `name`, for example `name=bridge.stl`.

#### Run command
Move folders `Config file/` `init/` to `build/`, move all environment meshes to `build/`. Run commands:
```bash
cd build
mkdir result
debugPathPlanning3D 3D.json name
```
Result information will be saved in `result/`.


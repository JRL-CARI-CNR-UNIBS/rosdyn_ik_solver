# Compute IK solutions using RosDyn #



## Description
This package computes all the IK solutions for a TF using RosDyn. __No collision checking are provided__

Services are described [here](https://github.com/JRL-CARI-CNR-UNIBS/ik_solver_msgs).


## Usage

run the node _rosdyn_ik_ providing the following parameter
### required (local) params
```yaml
base_frame: world # base frame of the chain
flange_frame: tool0 # end frame of the chain
tool_frame: open_tip # destination frame of the IK (it should be rigid attached to flange_frame)
desired_solutions: 32 # number of desired solution

joint_names: [jnt1, jnt2, ...]  #[optional] desired joints order
```
Example in [this launcher](launch/test.launch)

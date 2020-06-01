# MPC-controller
## Description
Eigen-based basic MPC-controller that solves the receding horizon problem with quadratic cost, linear dynamics and constraints. No theoretical stability or recursive feasibility guarantees! Additional QP solver is used that can be downloaded form [here](https://www.cs.cmu.edu/~bstephe1/eiquadprog.hpp). Use +/- std::numeric_limits<double>::infinity() as a maximal/minimal value for unconstrained states and inputs. 
## Demo
  main.cpp contains an example on how to use the controller. I modelled an inverted pendulum on a slide that can move along the horizontal axis in 2D. The input is a horizontal force applied to the slide. The goal of the controller is to keep the pendulum at its upper equilibrium position. Horizontal displacement and input force magnitude are limited. The equations of motion are linearized around the upper equilibrium.
## Animation
  I simulated the behavior of the system under the controller, the two gif files contain the results. I also added random disturbance to the system in form of random angular accelerations of the pendulum, in order to model the effect of wind for example.

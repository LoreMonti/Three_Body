# Three-Body Problem Simulation

Numerical simulation of the three-body problem applied to the **Sun–Earth–Jupiter system**, with a comparison between fixed-step and adaptive numerical integration methods.

---

## Objective

This project studies the dynamical behavior of a gravitational three-body system, a classical problem known for its **chaotic nature and lack of a general analytical solution**.

The main goals are:
- simulate the orbital evolution of the bodies over time  
- compare different numerical integration methods  
- analyze the **energy stability of the system**

---

## Numerical Methods

### 🔹 4th Order Runge-Kutta (RK4)
- fixed time step  
- good accuracy and simplicity  
- **unstable under strong perturbations**

When Jupiter’s mass is artificially increased, the system becomes unstable and the Earth is ejected.

---

### 🔹 Runge-Kutta Cash-Karp (RKCK)
- adaptive time step  
- local error estimation  
- improved numerical stability  

The method uses two embedded approximations (4th and 5th order) to estimate the local error:

$\epsilon = \left\| y^{(5)} - y^{(4)} \right\|$

and dynamically adjust the time step.

Result: better energy conservation and improved stability of the system.

---

## Results

- **RK4:**
  - significant energy drift  
  - instability under perturbations  
  - Earth ejection from the system  

- **RKCK:**
  - excellent energy conservation  
  - stable simulations even with modified masses  
  - higher efficiency due to adaptive step size  

Energy plots clearly show the superiority of RKCK over RK4.

---

## N-Body Extension

The model was extended to an N-body system including:
- Sun  
- major planets (up to Pluto)  

The simulation runs over **250 years** to capture Pluto’s full orbit.

Results:
- excellent energy conservation  
- very low relative error  
- dynamic step-size adaptation  

---

## Code Structure

The C++ implementation includes:
- gravitational equations of motion  
- RK4 integrator  
- RKCK integrator with error control  
- total energy computation of the system  

The integration method can be selected by modifying:

```cpp
#define METHOD RKCK
// or
#define METHOD RK_4

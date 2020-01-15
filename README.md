# Engineering study: Flight command for Dassault Mirage 2000

This project is my work on a use case for **"Flight Commands"** as a 5th Year **Engineering student at ESTACA**

## Abstract

This aim of this project was to take the **Dassault Aviation Mirage 2000** as a use case for Fligt command scaling. We had to do a full longitudinal and transversal aerodynamic study of the aircraft before scaling flight commands and then run simulations on it. We had to look over the aircraft's specification and scale our flight commands according to the architecture and required maneuvers. This project made us apply all our theoretical knowledge upon aerodynamics and commands on a concrete use case to understand fully what scaling an aircraft's flight command is like. We also got to see how the aircraft architecture and flight commands requirement influence one another.

## Summary

### A. Longitudinal Study

- **I. Aerodynamic Model**
	- 1. *Czα* coefficient estimation
	- 2. Aircraft center and *Cmα* coefficient estimation at 53%
	- 3. Static margin calculation of centered aircraft at 53%
	- 4. Center of gravity determination for 4%/0%/-4% stability
	- 5. *Czδm* coefficient estimation
	- 6. *Cmδm* coefficient estimation

- **II Aircraft Balancing**
	- 1. Incidence and steering calculation for the rudder in flight
	- 2. Incidence and steering calculation for the rudder in approach
	- 3. Static margin position influence

- **III. Small movement around the balance point**
	- 1. Dampening, static margin and eigenvalue calculation for both flight points
	- 2. Load factor calculation at the center of gravity

### B. Transversal Study

- **I. Aerodynami Model**
	- 1. *Cyβ, Cnβ* and *Cnr* coefficients calculation
	- 2. Result analysis
	- 3. Aircraft Balancing

### C. Flight Command

- **I. Command Law Synthesis**
	- 1. *Ga* gain adjustment
	- 2. *Gq* gain adjustment

- **II. Simulation Results**
	- 1. Aircraft incidence command
	- 2. Aircraft pitch speed command
	- 3. Aircraft load factor command

### D. Conclusion

## Credits

This work was done with my former school colleague **Lucas PFLIEGER**

Our teachers help and the ressources of **Dassault Aviation**

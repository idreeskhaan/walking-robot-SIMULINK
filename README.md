# SIMULINK model for a walking robot with ZMP based trajectory solver for torso

## Mathematical Modelling
<br>
<li> Plan the foot trajectory</li>
<li> Calculate z3 using TMIPM (Two Mass Inverted Pendulum Model)</li>
<li> Get joint angles using inverse kinematics solution</li>
<li>Compute the new trajectories for m4, m5 and m6 using Forward Kinematics</li>
<li>Conmpute new z3 based on the new trajectoires</li>
<li>Compare with previous and keep iterating until convergence</li>

![Alogorithm Flowchart](/algorithm.PNG)

### Model Foot Trajectory
![Foot Trajectory](/foot-trajectory.PNG)

### ZMP model using TMIPM
![Torso trajectory model using TMIPM](/zmp-TMIPM.PNG)

### ZMP model using MMIPM
![Torso trajectory model using MMIPM](/zmp-MMIPM.PNG)

### SIMULINK model
![Simulink Model](/simulink-model.PNG)

### Simulated Torso Trajectory
![Torso Trajectory](/torso-traj.PNG)

<br>
<br>
For more info, read report.pdf.

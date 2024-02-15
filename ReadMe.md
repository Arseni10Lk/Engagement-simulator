# Engagement-simulator #

The simulation occurs in 2D only, assumes head-on collision and does not consider manoeuvrability limitations of a pursuer. 

The program simulates the engagement scenarios in 3 steps: 
* Setting initial condition
* Application of Pro-Nav laws and integration with RK4
* Post-processing output

In case program compiles with the following initial conditions:

Parameter | Value
-----|---------
Pro-Nav type | True
Target acceleration | 0
Heading error | -20°
Proportianility const N | 3
∠(Target Velocity, horizon) |  0
Initial position of pursuer | (0 , 10000) m
Initial position of target | (40000,10000) m
Pursuer speed | 4000 m/s
Target speed | 1000 m/s

Then we can obtain two general visualizations:
Full engagement:

https://github.com/Arseni10Lk/Engagement-simulator/assets/141524111/d9014651-1566-4b15-bf9b-2f820ea5f915

Engagement zoomed on the pursuer: 

https://github.com/Arseni10Lk/Engagement-simulator/assets/141524111/b66eb4b1-8300-4b37-958b-ceecf4dba833



And visualizations of particular parameter values:
Parameter | Graph
-----|---------
Line of sight angle change rate | <img src = "./Media/LOSrate.png" alt = "Rate of change of Line of Sight angle" width = 500>
Pursuer to target distance (on a log scale) | <img src = "./Media/RTM.png" alt = "Pursuer to target distance" width = 500>
Closing Velocity | <img src = "./Media/VC.png" alt = "Closing velocity" width = 500>
Pursuer acceleration | <img src = "./Media/aP.png" alt = "Pursuer acceleration" width = 500>

In order to see results for different initial conditions you have to change the first section with the initial conditions (obviously) and also the third section wth visualizations.

For example, by editing several lines of the code we can obtain animations for the same initial conditions, but with
Parameter | Value
-----|---------
Pro-Nav type | Pure
Target acceleration | 20g

And let the acceleration be perpendicular to the velocity vector.

Thus, the results of the siulation will look in the following way:

Full engagement:

https://github.com/Arseni10Lk/Engagement-simulator/assets/141524111/06d3c6ec-99d9-4fe2-a485-0c0af9ec847f

Engagement zoomed on the pursuer: 

https://github.com/Arseni10Lk/Engagement-simulator/assets/141524111/58784efc-ea69-43b9-a663-0f0741030e5a



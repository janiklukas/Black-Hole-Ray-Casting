# Black Hole Ray Casting

The purpose of this project is to simulate a Schwarzschild black hole surrounded by a light-emitting accretion disc using ray casting: A set of light rays corresponding to pixels on an observer's screen is initialized and then cast using the geodesic equation for lightlike particles.

## Ray Initialization

The trajectory of a given light ray is described by a state vector with $8$ components:

$$X(\lambda)=(t(\lambda),r(\lambda),\theta(\lambda),\varphi(\lambda);\dot t(\lambda),\dot r(\lambda),\dot\theta(\lambda),\dot\varphi(\lambda))$$

where $\lambda$ is the curve parameter.

Each light ray starts at the camera location given by $(r_0,\theta_0,0)$ in Schwarzschild coordinates. The screen center is located some distance $d$ closer to the singularity, i.e. at $(r_0-d,\theta_0,0)$. To determine the locations of the pixels on the screen, a Cartesian coordinate system with the camera on its $x$-axis is used. The initial velocity of a given ray can then be calculated in the screen coordinates as the difference between the camera location and that of the associated pixel. 

The initial light rays are then transformed from the screen coordinates to flat spherical coodinates compatible with the Schwarzschild coordinates. To achieve this, we need to first rotate by $\theta_0$ about the $y$-axis before using the standard transformation expressions.

## Ray Evolution

The time evolution of $X(\lambda)$ is governed by the following equations derived from the geodesic equation:

$$\dot X_1=X_5,\quad\dot X_2=X_6,\quad\dot X_3=X_7,\quad\dot X_4=X_8;$$

$$\dot X_5=\frac{1}{X_2(1-X_2)}X_5X_6,$$

$$\dot X_6=\frac{1-X_2}{2X_2^3}(X_5)^2 + \frac{1}{2X_2(X_2-1)}X_6^2 + (X_2-1)[X_7^2+\sin^2(X_3)\,X_8^2]$$

$$\dot X_7=-\frac{2}{X_2} X_6 X_7 + \sin(X_3) \cos(X_3)\,X_8^2$$

$$\dot X_8= -\frac{2}{X_2} X_6 X_8 - 2\cot(X_3)\,X_7X_8$$

where we have set $c=1$ and $R_\text{S}=1$. The equation for $\dot X_6$ can be simplified using the lightlike condition, yielding

$$\dot X_6=\frac{1}{2}(2X_2-3)[X_7^2+\sin^2(X_3)\,X_8^2].$$

In order to enable detecting collisions with the accretion disc, both the current and previous state vector are saved at each time step. Thus a light ray is described by a $2\times 8$ matrix at any given time. After each step, we check if the ray has crossed the $\theta=\pi/2$ plane that contains the disc between steps. In this case, an interval search is performed to find the precise crossing point and check if it lies within the accretion disc. If yes, the pixel corresponding to the ray is coloured orange.

If any ray comes too close to the event horizon or moves too far away from the black hole, the evolution is stopped and the associated pixel is coloured black.

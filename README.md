This code is for simulating satellite in 3D Pyvista condition. 

It reads tle data from web and convert into satellite location data.Starting from current time, it will show orbit of satellite for 6000 steps.
Satellite will be expressed as stl image of certain satellite and show its own orthornormal body frame. Also, it can point certain vectors from its origin for example, sun vectors or LOS vectors. 

Ground station is predetermined vectors and convert into certain cartesian location. If satellite pass the ground station and its distance between ground station and satellite is smaller than 3000km, I assumed that the communication systems are available and visualize it.

Attitude of satellite is determined by DCM method from quaternions and quaternions should be predetermined as input value. With those input values, satellite can visualize its attitude as rotation of orthornormal body frame.

sat_trajectory[ECEF].py

sat_trajectory[ECEF] is simulation in ECEF i.e., earth centered and earth fixed frame. In this frame, orgin of coordinates is identical to center of the earth and z axis is identical with rotation axis of earth, x aixs is pointing to zero-lattitude line, and y axis follows orthonormal coordinates directions.
Since earth is fixed in this simulation, satellite update its location with modifying its orbit with earth's rotation. So, if we run this simulation longer than 1 orbit, the result will be looks like 'referecen_figure_1~5'


sat_trajectory[J2000].py

sat_trajectory[J2000] is simulation in J2000 frame. J2K frame is ECI frame that earth centered but, inertial frame that orthonormal coordinates reflect earth's rotation. z axis is identical to rotation axis, x axis is pointing equinox of 2000, and y axis follows orthonormal coordinates directions.
In this simulation, 'rot_angle' parameter reflect earth rotation angle after 2000 year. 

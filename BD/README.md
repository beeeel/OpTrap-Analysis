Simulate particle trajectories for trapped beads using Brownian Dynamics.

First use `func_BDopts` to create an options structure. This function takes name-value pairs as inputs and returns a struct as output. Then call `func_simulateDataField` to run the simulation.

Note: Currently viscosity must be time-independent (i.e. Newtonian), however if you pester me I can update this to create a version which is more generally applicable.

Example:

`opts = func_BDopts('radius',3e-6,'temp',37,'kappaNm',[1 1 0.2] * 1e-6, 'Nt', 5e5, 'eta', 0.75e-3,'output','tracks');
tracks = func_simulateDataField(opts);`

After the simulation has run, the variable `tracks` will be a cell array where each cell contains a Nt x 4 array, where each row is [t x y z].

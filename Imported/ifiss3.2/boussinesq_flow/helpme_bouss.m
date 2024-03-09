%HELPME_BOUSS Boussinesq flow problem interactive help 
%   IFISS scriptfile: DJS;   3 May 2012.
% Copyright (c) 2012 D.J. Silvester, M.D. Mihajlovic.
fprintf(' \n');
fprintf('To setup and evolve a Boussinesq problem\n'); 
fprintf('simply run the driver: unsteady_bouss_testproblem\n');
fprintf('Nonzero temperature boundary conditions are set \n'); 
fprintf('in the user-defined function: /diffusion/specific_bc.m\n');
fprintf('Velocity boundary conditions are set \n'); 
fprintf('in the user-defined function: /stokes_flow/specific_flow.m\n');
fprintf(' \n');
fprintf('Checkpointing is performed every 200 timesteps (the parameter \n'); 
fprintf('<nnt> in the driver stabtrBouss) --- to continue time integration\n');
fprintf('from the latest checkpoint, simply run\n');
fprintf('restart_bouss or restart_step_bouss\n');
fprintf(' \n');
fprintf('To visualise the latest checkpoint solution\n'); 
fprintf('run the script file: ppbouss_checkpointdata\n\n');
fprintf('To test iterative solvers on the discrete system arising at \n'); 
fprintf('a specific time step, run the script file\n');
fprintf('snapshot_solvebouss\n');
fprintf(' \n');

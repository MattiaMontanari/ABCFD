%HELPME IFISS interactive help facility
%   IFISS scriptfile: DJS; 3 May 2012.
% Copyright (c) 2012 D.J. Silvester, H.C. Elman, A. Ramage (see readme.m)

fprintf(' \n');
fprintf(' IFISS\n')
fprintf(' Copyright (c) 2012 by D.J. Silvester, H.C. Elman and A. Ramage (see readme.m)\n')
fprintf(' \n');
fprintf(' To install the toolbox, the script-file\n');
fprintf(' gohome.m must be edited to reflect the correct path\n');
fprintf(' to the ifiss home directory on the installed computer: \n')
fprintf(' \n %s\n\n',pwd);
%currentdir=pwd; fprintf(['\n ',currentdir,' \n\n']);
fprintf(' This step only needs to be done once. Once IFISS is installed,\n');
fprintf(' for all subsequent uses the MATLAB search path must include\n');
fprintf(' the IFISS subdirectories. This can be done by running the\n');
fprintf(' script-file setpath.m\n');
fprintf(' \n');

fprintf(' (Type any character to continue.)')
pause;

fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
fprintf('\nTry running a software demo to get started \n'); help ifissdemos

fprintf('\n For further information on\n');
fprintf('    special features                                  enter 0\n');
fprintf('    solving a diffusion problem                             1\n');
fprintf('    solving a convection-diffusion problem                  2\n');
fprintf('    solving a Stokes flow problem                           3\n');
fprintf('    solving a Navier-Stokes flow problem                    4\n\n');
fprintf('    exploring multigrid solvers                             6\n');
fprintf('    exploring preconditioned Krylov subspace solvers        5\n\n');
fprintf('    general information on time-dependent problems         11\n');
fprintf('    solving a Boussinesq flow problem                      12\n');
hlp=default('    Help topic',-1);
if hlp==1, gohome; cd diffusion; helpme_diff;  gohome, cd timestepping; helpme_heat; 
elseif hlp==2, gohome; cd convection; helpme_cd; gohome; cd timestepping; helpme_unsteady_cd; 
elseif hlp==3, gohome; cd stokes_flow; helpme_stokes; helpme_stopping
elseif hlp==4, gohome; cd navier_flow; helpme_navier; gohome; cd timestepping; helpme_unsteady_navier;
elseif hlp==5, gohome; cd solvers; helpme_it;
elseif hlp==6, gohome; cd solvers; helpme_mg;
elseif hlp==11, gohome; cd timestepping; helpme_timestepping;
elseif hlp==12, gohome; cd boussinesq_flow; helpme_bouss;
elseif hlp==0, 
   gohome; helpme_ch4;
   if exist('djs','file') == 7,
      gohome; helpme_djs;
   end
   if exist('hce','file') == 7,
      gohome; helpme_hce;
   end
   if exist('ar','file') == 7,
      gohome; helpme_ar;
   end
end
gohome;
fprintf(' \n');



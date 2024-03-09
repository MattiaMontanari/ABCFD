% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SUB-ROUTINE TO SET COMSOL's PARAMETER FOR A    %
%  /----\ |  \|    |--  |   |   GIVEN PHASE. NAME AND VALUES OF EACH PARAMTER  %
% /      \|__/ \__ |    |__/    AD DEFINED AD-HOC                              %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%abCFD_parameters( ph ) adds to the current model a set of scalar
%   parameters: three for Diriclet BC, four for both Neumann and Robin BCs,
%   other two parameter for the ramp function that initilize all the
%   boundary conditinos

function abCFD_parameters( model, ph , init_T )
%% TIME STEP PARAMETER
model.param.set('step_t', [num2str( ph.step_t ) , '[s]']) ;
model.param.descr('step_t', 'Time Step for phase 1');

slope_t_step = abCFD_SlopeRamp( 1 , ph.step_t * 0.95);
model.param.set('step_t_slope', num2str( slope_t_step )  ) ;
model.param.descr('step_t_slope', ' Slope for time step _ phase 1');

%% INITIAL CONDITIONS
% constant initial temperature field
model.param.set('initial_T', [num2str( init_T ) , '[K]']) ;
model.param.descr('initial_T', 'Initial temperature');

%% BOUNDARY CONDITIONS' PARAMETERS

% Diriclet boundary conditions  on \Gamma_1
model.param.set('BC_Diriclet_cutoff', num2str( ph.g1.cutoff ));
model.param.descr('BC_Diriclet_cutoff', 'Cut off value');

model.param.set('BC_Diriclet_slope', num2str( ph.g1.slope ));
model.param.descr('BC_Diriclet_slope', 'Slope value');

model.param.set('BC_Diriclet_t_initial', num2str( ph.g1.t_initial ));
model.param.descr('BC_Diriclet_t_initial', 'Trigger time');

% Robin boundary conditions on \Gamma_2
model.param.set('BC_Robin_h', num2str( ph.g2.pi ) );
model.param.descr('BC_Robin_h', 'heat transfer coeff.');

model.param.set('BC_Robin_pi', num2str( ph.g2.pi ) );
model.param.descr('BC_Robin_pi', 'conductvity for Robin');

model.param.set('BC_Robin_cutoff', num2str( ph.g2.cutoff ));
model.param.descr('BC_Robin_cutoff', 'Cut off value');

model.param.set('BC_Robin_slope', num2str( ph.g2.slope ));
model.param.descr('BC_Robin_slope', 'Slope value');

model.param.set('BC_Robin_t_initial', num2str( ph.g2.t_initial ));
model.param.descr('BC_Robin_t_initial', 'Trigger time');

% Neumann boundary conditions  on \Gamma_4
model.param.set('BC_Neumann_pi', num2str( ph.g4.pi ) );
model.param.descr('BC_Neumann_pi', 'conductvity for Neumann');

model.param.set('BC_Neumann_cutoff', num2str( ph.g4.cutoff ));
model.param.descr('BC_Neumann_cutoff', 'Cut off value');

model.param.set('BC_Neumann_slope', num2str( ph.g4.slope ));
model.param.descr('BC_Neumann_slope', 'Slope value');

model.param.set('BC_Neumann_t_initial', num2str( ph.g4.t_initial ));
model.param.descr('BC_Neumann_t_initial', 'Trigger time');

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 2.0                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   2.0 - Update to better match Comsol's requirements              13/02/2013 %     
%   1.0 - IMPLEMENTATION FOR FEW COMSOL PARAMETERS                  01/02/2013 %     
% ---------------------------------------------------------------------------- %
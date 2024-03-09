% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SUB-ROUTINE TO SET COMSOL's PARAMETER FOR A    %
%  /----\ |  \|    |--  |   |   GIVEN PHASE. NAME AND VALUES OF EACH PARAMTER  %
% /      \|__/ \__ |    |__/    AD DEFINED AD-HOC                              %                     
%  * * * CALLS * * *                                                           %
%           i. load('geo.mat')                                                 %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %

function abCFD_parameters( model, ph )

% Load data
load( [char(pwd),'/Database/geo.mat'] );

%% MATERIAL AND GEO PROPERTIES
model.param.set('L', [num2str(L),'[m]'] );
model.param.descr('L', 'Cylinder dyameter');
model.param.set('Tmax', [num2str(ph.Tmax),'[K]']);
model.param.descr('Tmax', 'Max INLET temperature');
model.param.set('Tinit', [num2str(ph.Tinit),'[K]']);
model.param.descr('Tinit', 'Refence/Bulk temperature at infinity (cold)');
model.param.set('Tsur', [num2str(ph.Tsur),'[K]']);
model.param.descr('Tsur', 'Max surface estimated temperature');
model.param.set('rho', [num2str(ph.mat.rho),'[kg/m^3]']');
model.param.descr('rho', 'Density');
model.param.set('mud', [num2str(ph.mat.mud),'[N*s/m^2]']) ;
model.param.descr('mud', 'Dynamic viscosity');
model.param.set('k', [num2str( ph.mat.k ),'[W/(m*K)]']);
model.param.descr('k', 'Thermal conductivity');
model.param.set('Cp', [num2str( ph.mat.Cp ),'[kJ/(kg*K)]']);
model.param.descr('Cp', 'Heat capacity');
model.param.set('beta', [num2str( ph.mat.thex ),'[1/K]']);
model.param.descr('beta', 'Thermal expansion coefficient');

%% TIME STEP PARAMETER
model.param.set('step_t', [num2str( ph.step_t ) , '[s]']) ;
model.param.descr('step_t', 'Time Step for phase 1');

%% HEAT TRANSFER USEFUL PARAMETERS
% INITIAL CONDITIONS
model.param.set('initial_T', 'Tinit' ) ;
model.param.descr('initial_T', 'Initial temperature');

% BOUNDARY CONDITIONS' PARAMETERS
% Diriclet
model.param.set(['HT_DIR_','cutoff'], num2str( ph.ht.DIR.start.cutoff ));
model.param.descr(['HT_DIR_','cutoff'], 'Cut off value');

model.param.set( ['HT_DIR_','slope'], num2str( ph.ht.DIR.start.slope ));
model.param.descr( ['HT_DIR_','slope'], 'Slope value');

model.param.set( ['HT_DIR_','t_initial'] , num2str( ph.ht.DIR.start.t_initial ));
model.param.descr( ['HT_DIR_','t_initial'] , 'Trigger time');

% Robin
model.param.set( ['HT_ROB_','h'], num2str( ph.ht.ROB.kappa ) );
model.param.descr(['HT_ROB_','h'], 'heat transfer coeff.');

model.param.set(['HT_ROB_','pi'], num2str( ph.ht.ROB.pi ) );
model.param.descr( ['HT_ROB_','pi'] , 'conductvity for Robin');

model.param.set( ['HT_ROB_','cutoff'], num2str( ph.ht.ROB.start.cutoff ));
model.param.descr( ['HT_ROB_','cutoff'] , 'Cut off value');

model.param.set( ['HT_ROB_','slope'] , num2str( ph.ht.ROB.start.slope ));
model.param.descr( ['HT_ROB_','slope'] , 'Slope value');

model.param.set( ['HT_ROB_','t_initial'] , num2str( ph.ht.ROB.start.t_initial ));
model.param.descr( ['HT_ROB_','t_initial'] , 'Trigger time');

% Neumann
model.param.set( ['HT_NEU_','pi'] , num2str( ph.ht.NEU.pi ) );
model.param.descr( ['HT_NEU_','pi'] , 'conductvity for Neumann');

model.param.set( ['HT_NEU_','cutoff'] , num2str( ph.ht.NEU.start.cutoff ));
model.param.descr( ['HT_NEU_','cutoff'] , 'Cut off value');

model.param.set( ['HT_NEU_','slope'] , num2str( ph.ht.NEU.start.slope ));
model.param.descr( ['HT_NEU_','slope'] , 'Slope value');

model.param.set( ['HT_NEU_','t_initial']  , num2str( ph.ht.NEU.start.t_initial ));
model.param.descr( ['HT_NEU_','t_initial'] , 'Trigger time');

%% SINGLE PHASE FLOW USEFUL PARAMETERS
% INITIAL CONDITIONS
model.param.set('initial_V', [num2str( ph.fd.init ) , '[m/s]']) ;
model.param.descr('initial_V', 'Initial velocity field');

% BOUNDARY CONDITIONS' PARAMETERS
% Diriclet
model.param.set(['fd_DIR_','cutoff'], num2str( ph.fd.DIR.start.cutoff ));
model.param.descr(['fd_DIR_','cutoff'], 'Cut off value');

model.param.set( ['fd_DIR_','slope'], num2str( ph.fd.DIR.start.slope ));
model.param.descr( ['fd_DIR_','slope'], 'Slope value');

model.param.set( ['fd_DIR_','t_initial'] , num2str( ph.fd.DIR.start.t_initial ));
model.param.descr( ['fd_DIR_','t_initial'] , 'Trigger time');

% Robin
model.param.set( ['fd_ROB_','h'], num2str( ph.fd.ROB.kappa ) );
model.param.descr(['fd_ROB_','h'], 'heat transfer coeff.');

model.param.set(['fd_ROB_','pi'], num2str( ph.fd.ROB.pi ) );
model.param.descr( ['fd_ROB_','pi'] , 'conductvity for Robin');

model.param.set( ['fd_ROB_','cutoff'], num2str( ph.fd.ROB.start.cutoff ));
model.param.descr( ['fd_ROB_','cutoff'] , 'Cut off value');

model.param.set( ['fd_ROB_','slope'] , num2str( ph.fd.ROB.start.slope ));
model.param.descr( ['fd_ROB_','slope'] , 'Slope value');

model.param.set( ['fd_ROB_','t_initial'] , num2str( ph.fd.ROB.start.t_initial ));
model.param.descr( ['fd_ROB_','t_initial'] , 'Trigger time');

% Neumann
model.param.set( ['fd_NEU_','pi'] , num2str( ph.fd.NEU.pi ) );
model.param.descr( ['fd_NEU_','pi'] , 'conductvity for Neumann');

model.param.set( ['fd_NEU_','cutoff'] , num2str( ph.fd.NEU.start.cutoff ));
model.param.descr( ['fd_NEU_','cutoff'] , 'Cut off value');

model.param.set( ['fd_NEU_','slope'] , num2str( ph.fd.NEU.start.slope ));
model.param.descr( ['fd_NEU_','slope'] , 'Slope value');

model.param.set( ['fd_NEU_','t_initial']  , num2str( ph.fd.NEU.start.t_initial ));
model.param.descr( ['fd_NEU_','t_initial'] , 'Trigger time');

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.4                                 date:  August 2013            % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.4 - Changed Tbul with Tinit for a more defined description    17/08/2013 %
%       moreover the position of the geo.mat file has been changed             %
%   0.3 - Changed reference .mat file: from 'input' to 'geo'        09/04/2013 %     
%       The new file is better defined and contains information only on the    %
%       geometry of the problem. BC's etc now is included into ph struct.      %
%   2.0 - Update to better match Comsol's requirements              13/02/2013 %     
%   1.0 - IMPLEMENTATION FOR FEW COMSOL PARAMETERS                  01/02/2013 %     
% ---------------------------------------------------------------------------- %
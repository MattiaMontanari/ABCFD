% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SUB-ROUTINE TO SET COMSOL's PARAMETER FOR A    %
%  /----\ |  \|    |--  |   |   GIVEN PHASE. NAME AND VALUES OF EACH PARAMTER  %
% /      \|__/ \__ |    |__/    AD DEFINED AD-HOC                              %
%                                                                              %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%abCFD_parameters( ph ) adds to the current model a set of scalar parameters

function abCFD_functions( model )

%% RAMP FUNCTIONS FOR BOUNDARY CONDITIONS
name_Robin = 'RobinRampValue';
name_Diriclet = 'DiricletRampValue';
name_Neumann = 'NeumannRampValue';

% Working on \Gamma_1
DiricletRamp = model.func.create(name_Diriclet, 'Ramp');
DiricletRamp.name( name_Diriclet );
DiricletRamp.set('funcname', 'g1');
DiricletRamp.set('location', 'BC_Diriclet_t_initial');             
DiricletRamp.set('cutoff', 'BC_Diriclet_cutoff');
DiricletRamp.set('slope', 'BC_Diriclet_slope');

RobinRamp = model.func.create( name_Robin , 'Ramp');
RobinRamp.name( name_Robin );                       % Name of the function shown
RobinRamp.set('funcname', 'g2');                    % function's tag
RobinRamp.set('location', 'BC_Robin_t_initial');             
RobinRamp.set('cutoff', 'BC_Robin_cutoff');
RobinRamp.set('slope', 'BC_Robin_slope');
 
NeumannRamp = model.func.create( name_Neumann , 'Ramp');
NeumannRamp.name( name_Neumann );                       % Name of the function shown
NeumannRamp.set('funcname', 'g4');                    % function's tag
NeumannRamp.set('location', 'BC_Neumann_t_initial');             
NeumannRamp.set('cutoff', 'BC_Neumann_cutoff');
NeumannRamp.set('slope', 'BC_Neumann_slope');

%% RAMP FUNCTIONS FOR TIME INITIALIZATION

InitializationRamp = model.func.create( 'InitializeRamp' , 'Ramp');
InitializationRamp.name( 'InitializeRamp' );                       % Name of the function shown
InitializationRamp.set('funcname', 'g_init');                    % function's tag
InitializationRamp.set('location', '0');             
InitializationRamp.set('cutoff', '1');
InitializationRamp.set('slope', 'step_t_slope');

%% INTERPOLATION FUNCTIONS FOR BOUNDARY FLUX EVALUATIONS

% On \Gamma1
IntFun = model.func.create('InterpGamma1', 'Interpolation');
IntFun.model('mod1');
IntFun.set('source', 'file');
IntFun.set('filename', [pwd , '\myfile.txt'] );
IntFun.setIndex('funcs', 'IntGamma1', 0, 0);
IntFun.set('defvars', 'on');
IntFun.set('extrap', 'value');
IntFun.set('argunit', 'm');
IntFun.set('fununit', 'K');   
model.func('InterpGamma1').active(false);

% On \Gamma2
IntFun = model.func.create('InterpGamma2', 'Interpolation');
IntFun.model('mod1');
IntFun.set('source', 'file');
IntFun.set('filename', [pwd , '\myfile.txt'] );
IntFun.setIndex('funcs', 'IntGamma2', 0, 0);
IntFun.set('defvars', 'on');
IntFun.set('extrap', 'value');
IntFun.set('argunit', 'm');
IntFun.set('fununit', 'K');  
model.func('InterpGamma2').active(false);

% On \Gamma3
IntFun = model.func.create('InterpGamma3', 'Interpolation');
IntFun.model('mod1');
IntFun.set('source', 'file');
IntFun.set('filename', [pwd , '\myfile.txt'] );
IntFun.setIndex('funcs', 'IntGamma3', 0, 0);
IntFun.set('defvars', 'on');
IntFun.set('extrap', 'value');
IntFun.set('argunit', 'm');
IntFun.set('fununit', 'K');  
model.func('InterpGamma3').active(false);

% On \Gamma4
IntFun = model.func.create('InterpGamma4', 'Interpolation');
IntFun.model('mod1');
IntFun.set('source', 'file');
IntFun.set('filename', [pwd , '\myfile.txt'] );
IntFun.setIndex('funcs', 'IntGamma4', 0, 0);
IntFun.set('defvars', 'on');
IntFun.set('extrap', 'value');
IntFun.set('argunit', 'm');
IntFun.set('fununit', 'K');  
model.func('InterpGamma4').active(false);

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 3.0                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.3 - Added functions for heat flux monitor on buondaries       25/02/2013 %
%   2.0 - Update to better match Comsol's requirements              13/02/2013 %
%   1.0 - IMPLEMENTATION FOR FEW COMSOL PARAMETERS                  01/02/2013 %     
% ---------------------------------------------------------------------------- %
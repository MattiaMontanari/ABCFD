% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	CHANGES THE MODEL OBJECT TO CHANGE THE PHYSICS %
%  /----\ |  \|    |--  |   |   IT SIMPLIFIES THE FULL-PHYSICS MODEL TO A      %
% /      \|__/ \__ |    |__/    STOKES FLOW AND HEAT TRANSFER IN SOLIDS        %
%                               PROBLEM. NEGRECTS BUOYANCY EFFECT AND ANY HEAT %
% SOURCE. ENSURES THE RIGHT TIME STEPS WILL BE TAKEN                           %
% * * * REMIND * * IT DOES NOT RUN THE SOLVER, IT JUST CHANGES THE MODEL OBJ   %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%abCFD_mdlobj_ChangePhysics( SolTag, Steps, model, mud, physics )
% INPUT:
%       - physics   : Defines which physics should be solved
%                       if == 'Simple' solves Stokes + heat conduction
%                       if == 'Full' solves Navier Stokes + heat in fluids
%                       if == 'NavierStokes' solves Nav-Stks + heat conduction
%       - Steps     : exact time steps to be taken
%       - model     : comsol's structure
%       - mud       : dynamic viscosity which might be decrased
%       - SolTag    : Sovler tag to activated: 
%                       if == 1 use Transient Solver   
%                       if == 2 use Stationary solver
function abCFD_mdlobj_ChangePhysics( model,SolTag, Steps, physics,SinPhasFlow,HeatTransfe )

% Initialize   
SolvTag = cell( model.sol.tags );
SolvTag = SolvTag{ SolTag };
% RelTol = min( Steps(Steps > 0) ) /10;

%% MODIFY MODELOBJECT
if 1 == strcmp( physics, 'Simple')

    % Take strict time steps
    model.sol('TranSol').feature('t1').set('tstepsbdf', 'strict');
    % Decrease relative tollerance
    model.study('std1').feature('time').set('rtolactive', 'on');
    model.study('std1').feature('time').set('rtol', 1e-10 ); 

    % Stokes flow
    model.physics( SinPhasFlow ).prop('StokesFlowProp').set('StokesFlowProp', 1, '1');
        % Decrease the dynamic viscosity
    model.physics( SinPhasFlow ).feature('fp1').set('mu_mat', 1, 'userdef');
    model.physics( SinPhasFlow ).feature('fp1').set('mu', 1, num2str( 1 ) );
        % Deactivate buoyancy effect
%     model.physics( SinPhasFlow ).feature('buoyancy').active(false);
    % Heat Transfer in Solid
    model.physics( HeatTransfe ).feature('fluid1').active(true);
    model.physics( HeatTransfe ).feature('fluid1').set('minput_velocity_src', 1, 'userdef');
    model.physics( HeatTransfe ).feature('fluid1').set('minput_velocity', {'0' '0' '0'});
        % Deactivate heat sources
%     model.physics( HeatTransfe ).feature('heatsource').active(false);
    % Impose stepes to the solver
    model.study('std1').feature('time').set('tlist', Steps' );
    if SolTag == 1
        model.sol( SolvTag ).feature('t1').set('tstepsbdf', 'intermediate');
        model.sol( SolvTag ).feature('t1').set('control', 'user');
        model.sol( SolvTag ).feature('t1').set('tlist', Steps' );
        % Control output solution
        model.sol( SolvTag ).feature('t1').set('tout', 'tlist');
    end
    
elseif 1 == strcmp( physics, 'NavierStokes')
    model.physics( SinPhasFlow ).prop('StokesFlowProp').set('StokesFlowProp', 1, '0');

elseif 1 == strcmp( physics, 'Full')
    model.physics( SinPhasFlow ).prop('StokesFlowProp').set('StokesFlowProp', 1, '0');
    model.physics( SinPhasFlow ).feature('buoyancy').active(true);
    model.physics( HeatTransfe ).feature('heatsource').active(true);
    model.physics( HeatTransfe ).feature('fluid1').set('minput_velocity_src', 1, 'root.mod1.u');

end

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   0.1 - kick-off                                                  21/03/2013 %
% ---------------------------------------------------------------------------- %
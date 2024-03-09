% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \  	CHANGES THE MODELOBJECT TO COMPUTE THE         %
%  /----\ |  \|    |--  |   |   REFERENCE SOLUTION: STOKES FLOW + HEAT         %
% /      \|__/ \__ |    |__/    CONDUCTION.                                    %
%                               TRANSIENT SOLVER WITH STABILIZATION TECHNIQUES %
% WITH SAME BOUNDARY CONDITIONS AS THE FULL PHYSICAL MODEL                     %
%  * * * CALLS * * *                                                           %
%           i. abCFD_getU                                                      %
%          ii. abCFD_mdlobj_SelectSolStab                                      %
%         iii. abCFD_mdlobj_ChangePhysics                                      %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ K, D]= abCFD_assembl_MVE( model, assemb, mesh)
% INPUTS:
%   - model     : common phase structure
%   - mesh      : ordinary mesh structure
%   - TimeSteps : Output time values taken by the solver
% OUTPUT:
%   - R         : Reference solution, in a typical structure format
 
function [ R ] = abCFD_REF_solver( model, mesh, TimeSteps )

%% INITIALIZE VARIABLES

tiin = toc;
fprintf('Restart Comsol for Reference solution at %s  \n', datestr(now,15) )
 
%% STABILITY AND SOLVER SELECTION

abCFD_mdlobj_SelectSolStab( model, 1 , 'spf_FEM', 'ht_FEM', 'Stable' )

%% APPLY BOUNDARY CONDITIONS
 
%% CHANGE PHYSICS

abCFD_mdlobj_ChangePhysics( model, 1, TimeSteps, 'Simple','spf_FEM', 'ht_FEM' )

%% CHANGE SOLVER SETTING 

    % Physics to solve for
model.study('std1').feature('time').set('physselection', 'spf_FEM');
model.study('std1').feature('time').set('activate', {'spf_FEM' 'on' 'ht_FEM' 'on' 'spf_AUX' 'off' 'ht_AUX' 'off'});
    % convergence criteria
model.sol('TranSol').feature('t1').feature('fc1').set('ntermconst', 'tol');
model.sol('TranSol').feature('t1').feature('fc1').set('maxiter', '9');
    % Compile equations
model.sol('TranSol').feature('v1').set('control', 'time');
    % Time stepping control
model.sol('TranSol').feature('t1').set('control', 'time');

%% STUDY SETTING

    % relative tolerance
model.study('std1').feature('time').set('rtolactive', 'off');
    % Deactivate AUX physics'
model.study('std1').feature('time').set('useadvanceddisable', 'on');
model.study('std1').feature('time').set('disabledphysics', {'spf_AUX' 'ht_AUX'});


%% RUN THE SOLVER

model.sol('TranSol').runAll;

%% EXTRACT REFERENCE SOLUTION
[ R ] = abCFD_getU( model, 'TranSol', mesh );

%% EXIT

fprintf('Restart Comsol for Reference solution after %s  \n', datestr(toc-tiin,15) )

end
% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  April 2013             %
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick-off                                                  08/04/2013 %
% ---------------------------------------------------------------------------- %
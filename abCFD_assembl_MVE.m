% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \  	CHANGES THE MODELOBJECT TO COMPUTE THREE       %
%  /----\ |  \|    |--  |   |   COEFFICIENT MATRICES.                          %
% /      \|__/ \__ |    |__/    SOLVE A TRANSIENT SOLVER, WITHOUT STABILIZATION%
%                               TECHNIQUES, AND WITH SPECIFIC B.C.'s.          %
%  * * * CALLS * * *                                                           %
%           i. abCFD_MatrixAssembly                                            %
%          ii. abCFD_mdlobj_SelectSolStab                                      %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ M, V, E ] = abCFD_assembl_MVE( model, assemb, mesh)
% INPUTS:
%   - model     : common phase structure
%   - assemb    : Structure defining the matrices to be assembled
%   - mesh      : ordinary mesh structure
%   - TotDofs   : Total degrees of freedom 

function [M, V, E] = abCFD_assembl_MVE( model, assemb, mesh, TotDofs)

%% INITIALIZE VARIABLES
% Tags of unknown fields on domains
ta1 = mesh{ 2 }.ele{ 1 }.dof{ 1 }.uniqTAG ; % temperature
ta2 = mesh{ 2 }.ele{ 1 }.dof{ 2 }.uniqTAG ; % pressure (not used)
ta3 = mesh{ 2 }.ele{ 1 }.dof{ 3 }.uniqTAG ; % x-velocity
ta4 = mesh{ 2 }.ele{ 1 }.dof{ 4 }.uniqTAG ; % y-velocity
 
%% STABILITY AND TranSol
    % Remove stability and selec TranSol
abCFD_mdlobj_SelectSolStab( model, 1 , 'spf_AUX', 'ht_AUX', 'UnStable' )

%% APPLY BOUNDARY CONDITIONS
    % HT  - Deactivate diriclet and Neumann
model.physics('ht_AUX').feature('Diricl').active(false); 
model.physics('ht_AUX').feature('Neumann').active(true); %TRUE
model.physics('ht_AUX').feature('Robin').active(false);
    % SPF - Deactivate outlet
model.physics('spf_AUX').feature('out1').active(false);
    % SPF - All Diriclet BC, velocity field = [ 1 ,1 ]
model.physics('spf_AUX').feature('inl1').active(true);
model.physics('spf_AUX').feature('inl1').set('ComponentWise', 1, 'VelocityFieldComponentWise');
model.physics('spf_AUX').feature('inl1').set('u0', {'1' '1' '0'});
model.physics('spf_AUX').feature('inl1').selection.all;
    % PRESSURE
    model.physics('spf_AUX').feature.create('prpc1', 'PressurePointConstraint', 0);
model.physics('spf_AUX').feature('prpc1').selection.set([8]);


%% CHANGE PHYSICS
% model.physics('spf_AUX').prop('StokesFlowProp').set('StokesFlowProp', 1, '1');
% might be not necessary!
 
%% CHANGE SOLVER SETTING 
    % Convergence criteria
model.sol('TranSol').feature('t1').feature('fc1').set('ntermconst', 'iter');
    % Time steps
model.study('std1').feature('time').set('tlist', 'range(0,0.5,1)');
    % Physics to solve for
model.study('std1').feature('time').set('physselection', 'spf_FEM');
model.study('std1').feature('time').set('activate', {'spf_FEM' 'on' 'ht_FEM' 'off' 'spf_AUX' 'on' 'ht_AUX' 'on'});
model.study('std1').feature('time').set('physselection', 'spf_FEM');
model.study('std1').feature('time').set('activate', {'spf_FEM' 'off' 'ht_FEM' 'off' 'spf_AUX' 'on' 'ht_AUX' 'on'});

%% RUN THE SOLVER
model.sol('TranSol').runAll;

%% ASSEMBLE MATRICES
[ TrUnsFul_MVE ] = abCFD_MatrixAssembly( model, assemb, 1 );

%%  SAVE MATRICES: M & V & E
% Initialize
E = sparse( TotDofs , TotDofs );
M = sparse( TotDofs , TotDofs );
V = sparse( TotDofs , TotDofs );
% Extract
E( ta1, ta1) = TrUnsFul_MVE.D( ta1, ta1 );
% M( unique( [ta3,ta4],'stable' ),unique( [ta3,ta4],'stable' )) = ...
%     TrUnsFul_MVE.D( unique( [ta3,ta4],'stable' ),unique( [ta3,ta4],'stable' ) );
% V( ta1, unique( [ta3,ta4],'stable' )) = TrUnsFul_MVE.K( ta1, unique( [ta3,ta4],'stable' ) );
M( [ta3;ta4], [ta3;ta4] ) = TrUnsFul_MVE.D( [ta3;ta4] , [ta3;ta4] );
V( ta1, [ta3;ta4] ) = TrUnsFul_MVE.K( ta1, [ta3;ta4] );

% E = E ( [ta1;ta2;ta3;ta4] , [ta1;ta2;ta3;ta4] );
% M = M ( [ta1;ta2;ta3;ta4] , [ta1;ta2;ta3;ta4] );
% V = V ( [ta1;ta2;ta3;ta4] , [ta1;ta2;ta3;ta4] );

end
% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.3                                 date:  April 2013             %
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.3 - Added pressure point to avoid matrix singularity          11/04/2013 %
%       And changed ordering of output matrices
%   0.2 - Avoided the use of 'unique'                               09/04/2013 %
%   0.1 - kick-off                                                  08/04/2013 %
% ---------------------------------------------------------------------------- %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \  	CHANGES THE MODELOBJECT TO COMPUTE TWO         %
%  /----\ |  \|    |--  |   |   COEFFICIENT MATRICES.                          %
% /      \|__/ \__ |    |__/    SOLVE A STEADY SOLVER, WITHOUT STABILIZATION   %
%                               TECHNIQUES, AND WITH SPECIFIC B.C.'s.          %
%  * * * CALLS * * *                                                           %
%           i. abCFD_MatrixAssembly                                            %
%          ii. abCFD_mdlobj_SelectSolStab                                                 %
%         iii. abCFD_mdlobj_ChangePhysics                                                 %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ K, D]= abCFD_assembl_MVE( model, assemb, mesh)
% INPUTS:
%   - model     : common phase structure
%   - assemb    : Structure defining the matrices to be assembled
%   - mesh      : ordinary mesh structure
%   - TimeSteps : Output time values taken by the solver
%   - TotDofs   : Total degrees of freedom 
% OUTPUT:
%   - K         : Sparse matrix due Laplace operator in Energy equation 
%   - D         : Sparse matrix due 'Laplace operator' in Momentum equation 

function [ K, D] = abCFD_assembl_KD( model, assemb, mesh, TimeSteps, TotDofs)

%% INITIALIZE VARIABLES
% Tags of unknown fields on domains
ta1 = mesh{ 2 }.ele{ 1 }.dof{ 1 }.uniqTAG ; % temperature
ta2 = mesh{ 2 }.ele{ 1 }.dof{ 2 }.uniqTAG ; % pressure (not used)
ta3 = mesh{ 2 }.ele{ 1 }.dof{ 3 }.uniqTAG ; % x-velocity
ta4 = mesh{ 2 }.ele{ 1 }.dof{ 4 }.uniqTAG ; % y-velocity
 
%% STABILITY AND SOLVER SELECTION
    % Remove stability and selec StatSol
abCFD_mdlobj_SelectSolStab( model, 2 , 'spf_AUX', 'ht_AUX', 'UnStable' )

%% APPLY BOUNDARY CONDITIONS
    % Deactivate diriclet for heat 
model.physics('ht_AUX').feature('Diricl').active(true);
    % Deactivate diriclet for single phase flow
model.physics('spf_AUX').feature('inl1').active(true);

%% CHANGE PHYSICS

abCFD_mdlobj_ChangePhysics( model, 2, TimeSteps, 'Simple','spf_AUX', 'ht_AUX' )

%% CHANGE SOLVER SETTING 
    % Convergence criteria
model.sol('StatSol').feature('s1').feature('fc1').set('ntermauto', 'itertol');
    % Physics to solve for
model.study('std1').feature('stat').set('physselection', 'spf_FEM');
model.study('std1').feature('stat').set('activate', {'spf_FEM' 'on' 'ht_FEM' 'off' 'spf_AUX' 'on' 'ht_AUX' 'on'});
model.study('std1').feature('stat').set('physselection', 'spf_FEM');
model.study('std1').feature('stat').set('activate', {'spf_FEM' 'off' 'ht_FEM' 'off' 'spf_AUX' 'on' 'ht_AUX' 'on'});

%% RUN THE SOLVER
model.sol('StatSol').runAll;

%% ASSEMBLE MATRICES
[ StUnsSim_KD ] = abCFD_MatrixAssembly( model, assemb, 2 );

%%  SAVE MATRICES: M & V & E
% Initialize
K = sparse( TotDofs , TotDofs );
D = sparse( TotDofs , TotDofs );
% Extract
K( ta1, ta1) = StUnsSim_KD.K( ta1, ta1 );
% D( unique( [ta3,ta4],'stable' ),unique( [ta3,ta4],'stable' )) = ...
%     abCFD_round(1e-12, StUnsSim_KD.K( unique( [ta3,ta4],'stable' ),unique( [ta3,ta4],'stable' ) ));
D([ta3;ta4],[ta3;ta4] ) = abCFD_round( StUnsSim_KD.K( [ta3;ta4] ,[ta3;ta4] ),1e-12 );

% K = K ( [ta1;ta2;ta3;ta4] , [ta1;ta2;ta3;ta4] );
% D = D ( [ta1;ta2;ta3;ta4] , [ta1;ta2;ta3;ta4] );

end
% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  April 2013             %
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.2 - Changed ordering output matrices                          10/04/2013 %
%   0.1 - kick-off                                                  08/04/2013 %
% ---------------------------------------------------------------------------- %
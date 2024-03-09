% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \  	CHANGES THE MODELOBJECT TO COMPUTE ONE         %
%  /----\ |  \|    |--  |   |   COEFFICIENT MATRIX.                            %
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
%   - Cx        : Sparse matrix due non-linear term in Momentum equation 
%             computed with test constant flow field [u,v] = [ 1 , 1 ]    
%   - Cy        : Sparse matrix due non-linear term in Momentum equation 
%             computed with test constant flow field [u,v] = [ 2 , -.1 ]    
 
function [ Cx, Cy ] = abCFD_assembl_C( model, assemb, mesh, TimeSteps, TotDofs, D)

%% INITIALIZE VARIABLES
% Tags of unknown fields on domains
ta1 = mesh{ 2 }.ele{ 1 }.dof{ 1 }.uniqTAG ; % temperature (not used)
ta2 = mesh{ 2 }.ele{ 1 }.dof{ 2 }.uniqTAG ; % pressure (not used)
ta3 = mesh{ 2 }.ele{ 1 }.dof{ 3 }.uniqTAG ; % x-velocity
ta4 = mesh{ 2 }.ele{ 1 }.dof{ 4 }.uniqTAG ; % y-velocity
 
%% STABILITY AND SOLVER SELECTION
    % Remove stability and selec StatSol
abCFD_mdlobj_SelectSolStab( model, 2 , 'spf_AUX', 'ht_AUX', 'UnStable' )

%% CHANGE PHYSICS

abCFD_mdlobj_ChangePhysics( model, 2, TimeSteps, 'NavierStokes','spf_AUX', 'ht_AUX' )

%% CHANGE SOLVER SETTING 
    % Convergence criteria
model.sol('StatSol').feature('s1').feature('fc1').set('ntermauto', 'itertol');
    % Physics to solve for
model.study('std1').feature('stat').set('physselection', 'spf_FEM');
model.study('std1').feature('stat').set('activate', {'spf_FEM' 'on' 'ht_FEM' 'off' 'spf_AUX' 'on' 'ht_AUX' 'on'});
model.study('std1').feature('stat').set('physselection', 'spf_FEM');
model.study('std1').feature('stat').set('activate', {'spf_FEM' 'off' 'ht_FEM' 'off' 'spf_AUX' 'on' 'ht_AUX' 'on'});

%% APPLY BOUNDARY CONDITIONS - make C(u,v) => C( u ) only
    % Velocity field U = [ 1 , 0 ]
model.physics('spf_AUX').feature('inl1').set('u0', {'1' '1' '0'});
model.physics('spf_AUX').feature('init1').set('u2', {'1' '1' '0'});

%% RUN THE SOLVER

model.sol('StatSol').runAll;

%% ASSEMBLE MATRICES

[ StUnsNavstk_Cx ] = abCFD_MatrixAssembly( model, assemb, 2 );

%% APPLY BOUNDARY CONDITIONS - make C(u,v) => C( v ) only
    % Velocity field U = [ 0 , 1 ]
model.physics('spf_AUX').feature('inl1').set('u0', {'2' '-.1' '0'});
model.physics('spf_AUX').feature('init1').set('u2', {'2' '-.1' '0'});

%% RUN THE SOLVER

model.sol('StatSol').runAll;

%% ASSEMBLE MATRICES

[ StUnsNavstk_Cy ] = abCFD_MatrixAssembly( model, assemb, 2 );


%%  SAVE MATRICES: Cx & Cy
% Initialize
Cx = sparse( TotDofs , TotDofs );
Cy = sparse( TotDofs , TotDofs );
% Assemble
Cx( ta3 ,ta3 ) = StUnsNavstk_Cx.K( ta3 ,ta3 ) ;
Cy( ta3 ,ta3 ) = StUnsNavstk_Cy.K( ta3 ,ta3 ) ;

Cx( ta3 ,ta3 ) = abCFD_round(Cx( ta3 ,ta3 ) - D( ta3 ,ta3 ) , 1e-12);
Cy( ta3 ,ta3 ) = abCFD_round(Cy( ta3 ,ta3 ) - D( ta3 ,ta3 ) , 1e-12);


end
% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.3                                 date:  April 2013             %
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.3 - Assemble one block per time                               10/04/2013 %
%   0.2 - C is now provided as a 'single' block, that is: it is a   09/04/2013 %
%       a block C(u,v) which is valid for ta3 as well as ta4                   %
%   0.1 - kick-off                                                  08/04/2013 %
% ---------------------------------------------------------------------------- %
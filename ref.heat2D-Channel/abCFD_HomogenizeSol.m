% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \    SUB-ROUTINE TO COMPUTE AND EVENTUALLY PLOT     %
%  /----\ |  \|    |--  |   |   THE HOMOGENOUS SOLUTION FIELD OF U             %
% /      \|__/ \__ |    |__/    ELEMENT PROBLEMS ARE SET UP TO FOR THE HEAT    %                     
%                               EQUATION PROBLEM AND STOKES FLOW               %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[  ]=abCFD_HomogenizeSol( 

function [ R ] = abCFD_HomogenizeSol( U, mesh, ph1 )
%%  PROBLEM STATEMENT
% The main is to compute the solution vector of a multiphysics problem
%                    
%                   | u |
%   Um = Hm + Rm =  | v |
%                   | p |
%                   | T | @ t = tm
% Where:
%   Um is the imported solution time
%   Hm is the homogenous solution unknown
%   Rm is an abritrary reference solution which gives:
% 
%                   Hm = Um - Rm
% 
% Boundary conditions for the thermal problem
%          < - - - -  Um - - - - - > < - - - -  Hm - - - - - - > < - - Rm - - >
% \gammaD |  T = fD(x,t)            |   T = 0                   |  T = fD(x,t) |
% \gammaN | dT/dn = fN(x,t)         |   dT/dn = fN(x,t)         |  dT/dn = 0   |
% \gammaR | 1/a*dT/dn+c*T = b(x,t)  |   1/a*dT/dn+c*T = b(x,t)  |  T = 0       |

% Boundary conditions for the Stokes problem
% 
%  ? ? ? ? ? ?

%%      HEAT PROBLEM

% Find variable Index for Temperature field
Field = 'mod1.T';
IndexC = strfind( var.fieldNames , Field );
Index = find(not(cellfun('isempty', IndexC)));
% Number of unknown for T
dof = var.fieldNDofs( Index );

% Find nodes' tag for temperature field
Field = 'mod1_T';
if strcmp( char( mesh{2}.ele{1}.dof{1}.name ) , Field )
    Index = 1;
else
%     investigate if Index is 2 or 3 or 4
    return
end
allTags = mesh{2}.ele{1}.dof{Index}.TAG;

% Obaint Stiffness matrix from Comsol
K = assemb.K( allTags, allTags);

% Find \gammaD
DirGamma = ph1.ht.Diriclet;
% Initialize
% DirGamma = [];
for id_gamma = 1 : numel(DirGam) 
    DirTag{ id_gamma } = mesh{1}.ele{ DirGamma }.dof{ Index }.TAG; %#ok<AGROW>
end

% Construct null space constrain 
% Indepenedent nodes
NullInd = eye( dof );               
NullInd( :, unique( DirTag ) ) = [];
% Dependent nodes
IndTag = [ 1 : dof ]';
IndTag ( unique( DirTag ) ) = [];
NullDir = eye( dof );               
NullDir( :, unique( IndTag ) ) = [];

% Delete Stiffness matrix
Kind = NullInd' * K * NullInd;
KDir = NullDir' * K * NullDir;

% SOLVE LINEAR SYSTEM FOR BETA
%  [Kind] * [NullInd' * R] = [KDir] * [NullDir' * R]

% Initialize
R = zeros( dof, U.solinfo.sizesolvals );

for i_tm  = 1 : U.solinfo.sizesolvals

    R( DirTag , i_tm) =  U.T.d1( i_tm ,:)';
    R( IndTag , i_tm) = Kind \ ( - KDir * R( DirTag , i_tm) );
    
end

% % COMPUTE HOMOGENEOUS FIELD
% Uhom = bsxfun(@times,U,( 1 - beta.field) );
% % COMPUTE Reference FIELD
% Uref = bsxfun(@times,U,( beta.field) );



%% % COMPUTE SOLUTION FIELD OF BETA AND HOMOGENIZE SOLUTION
% % Thom = T - Tref = T - Beta * T = T ( 1 - Beta )
% %   Where Beta(x) is a scalar field constant in time. It is computed
% %   solving the Laplace equation:     d^2( Beta ) = 0
% %   with BC´s:  Beta = 1 on SIGMA_2  & Beta = 0 elsewhere
% %   Using the mesh from Comsol I can use the stiffness matrix:
% 
% % ANALYZE EDGES AND DEGREES OF FREEDOM
% beta.dof.dof = dof.number_V1;
% 
% % Find the Diriclet nodes
% beta.Edg.tag = [ 1 , 2 , 3 , 4 ];
% beta.Edg.num = size( beta.Edg.tag , 2);
% 
% beta.dof.diric = [ dof.gamma1.DofTags; dof.gamma2.DofTags; ...
%                    dof.gamma3.DofTags; dof.gamma4.DofTags ];
% beta.dof.diric = unique( beta.dof.diric ); % Tags of diriclet nodes for beta
% 
% beta.dof.ind = (1 : beta.dof.dof)';
% beta.dof.ind ( beta.dof.diric ) = [];
% 
% % BUILD THE NULL MATRICES
% beta.NullDir = sparse( beta.dof.diric , 1:(size( beta.dof.diric,1)) ,1,beta.dof.dof,(size(beta.dof.diric,1)));
% beta.Nullind = sparse( beta.dof.ind , 1: size(beta.dof.ind,1), 1 , beta.dof.dof , size(beta.dof.ind,1));
% % APPLY REDUCED STIFFNESS MATRIX
% beta.KKind = beta.Nullind' * KK * beta.Nullind;
% beta.KKDir = beta.Nullind' * KK * beta.NullDir;
% % SOLVE LINEAR SYSTEM FOR BETA
% %  [NullBeta]' * [KK] * [NullBeta] * [NullBeta]' * [Beta] = 0
% %           [KKind]        *       [Betaind]  =  [KKDir] * [BetaDir]
% 
% %Initialize
% beta.field = zeros(beta.dof.dof,1);   
% % Apply Boundary condition BETA = 1 on \Gamma1
% beta.field( dof.gamma1.DofTags , 1) = 1;
% % Solve for the unknowns nodes of Beta
% beta.field( beta.dof.ind , 1 ) = beta.KKind \ (- beta.KKDir * beta.field( beta.dof.diric , 1));
% % Apply the exponential to beta.field
% beta.field = ( beta.field ).^pow;
% % COMPUTE HOMOGENEOUS FIELD
% Uhom = bsxfun(@times,U,( 1 - beta.field) );
% % COMPUTE Reference FIELD
% Uref = bsxfun(@times,U,( beta.field) );
% 
% % EVENTUALLY PLOT HOMOGENEOUS FIELD
% abCFD_PlotExtMesh( model, A, dof.dofCoord', ele.DOFToNode',Uhom, elements)


end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick-off                                                  19/03/2013 %     
%  
% ---------------------------------------------------------------------------- %
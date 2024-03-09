% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SUB-ROUTINE TO SOLVE A LINEAR SYSTEM WITH A    %
%  /----\ |  \|    |--  |   |   GENERAL PURPOSE. WORKING ROUTINE FOR FEM       %
% /      \|__/ \__ |    |__/    THE ASSEMBLED LINEAR SYSTEM HAS A GENERAL      %
%                               FORMULATION SUITABLE FOR FEM AND POD.          %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ U ] = abCFD_HybSolver( assemb, dof, B, te, ph,U,init_T, DirTag) 
% output: U       - is a stucture with new fields: .fem with the solution of
%                   the linear system and, .timing list of time stepts
%                   taken by the solver
% inputs: B       - is the maxtrix which defines the size of the linear system 
%         te      - tetha coefficiente for the time-discretization scheme
%         init_T  - scalar which defines the homogeneous initial condition
%         DorTag  - DOF's tags prescribed by Diriclet's boundary conditions

function [ U ] = abCFD_FemSolver( assemb, dof, te, ph,U,init_T, DirTag)
%% Assembly the following linear system
% 
%                  ?!?!?!?!?

%% INITILIAZIE VARIABLES
% Number of time steps in current phase
nt = ph.nt; 
% Numer of ALL dofs for physical problem (in general DofNum ~= ee )
DofNum = dof.number_V1;
% Numer of PRESCRIBED (dir) dofs
DirNum = abs(diff(size(assemb.Null)));
% Numer of INDEPENDENT (ind) dofs
IndNum = min(size(assemb.Null));
% Unknown array
U.fem = ones( DofNum, nt );
% Unknown delted array
a = ones( IndNum, nt );

%% ASSEMBLY DELETED MATRICES

% Coefficient deleted matrices
D = assemb.Null' * assemb.D * assemb.Null;
K = assemb.Null' * assemb.K * assemb.Null;        % Pure Stiffness matrix
KR = assemb.Null' * assemb.KRobin * assemb.Null ; % Robin's Stiffness contribute
% F = assemb.Null' * 0;
SR = assemb.Null' * assemb.SRobin;
SN = assemb.Null' * assemb.SNeuma;

% Coefficient reduced-deleted matrices
A = (1/ph.step_t * D +   te   * ( K + KR ) ); 
Z = (1/ph.step_t * D - (1-te) * ( K + KR ) );

%% APPLY INITIAL AND BOUNDARY CONDITIONS

% Initial conditions for Unknown delted array
a(:,1) = a(:,1) .* init_T;
% Initial conditions for Unknown array
U.fem(:,1) = U.fem(:,1) .* init_T;
% COMSOL's INTERPOLATION for boundary conditions for Unknown array
U.fem(DirTag, :) = U.U( DirTag , : );

%% RETRIEVE DELETED CONTRIBUTES
Kdir = zeros( IndNum , DirNum );
Ddir = zeros( IndNum , DirNum );

for id_dir = 1 :DirNum      % For each prescribed dof
    
    K_col = assemb.K( : , DirTag( id_dir ) ) ;    % K's deleted column 
    K_col( DirTag , : ) = [];              % delete the row of the prescribed dof
    Kdir( : , id_dir) = K_col;            % add the contribute
    
    KR_col = assemb.KRobin( : , DirTag( id_dir ) );    % KR's deleted column 
    KR_col( DirTag , : ) = [];              % delete the row of the prescribed dof
    KRdir( : , id_dir) = KR_col;            % add the contribute
    
    D_col = assemb.D( : , DirTag( id_dir ) );    % D's deleted column 
    D_col( DirTag , : ) = [];              % delete the row of the prescribed dof
    Ddir( : , id_dir) = D_col;            % add the contribute
    
end
% These are not correct:  they have to right to exist!!!
% SRdir = assemb.SRobin(DirTag);
% SNdir = assemb.SNeuma(DirTag);

% Coefficient reduced-retrieved matrices
Adir = 1/ph.step_t * Ddir +   te  * ( Kdir + KRdir );
% Zdir = 1/ph.step_t * B' * Ddir + (1-te) * B' * ( Kdir + KRdir );

%% COMPUTE SOLUTION IN TIME
%   The linear system to be solved for a(t) being t [ tinit + ph.step_t , ph.t_final ]
%   runs for nt+1-1 time steps because a is computed starting from the 
%   1sh physical solution. 

U.timing = zeros(1,nt);
for id_t = 2: 1 : nt + 1 - 1
    % Time instant computed within the id_t-th iteration
    U.timing = U.timing + ph.step_t;
    % Assemble the right hand side vector
%     f = BZB * a(: , id_t - 1) - BAB * U.ref(: , id_t) + BZB * U.ref(:, id_t -1 ) +...
%         B' * ( F );
    f = Z * a(: , id_t - 1) + te.*(SR + SN) + (1-te).* (SR + SN);
    DirCon = -1* Adir * U.fem(DirTag, id_t) ;%+ Zdir *U.fem(DirTag, id_t - 1) ;
    a(: , id_t ) = A \ (f + DirCon);
    
end

%% MAP DELETEDO SOLUTION a(t) INTO U(t)
% All dof's tag
IndTag = 1 : DofNum;
% find independent dof's Tags
IndTag( DirTag ) = 0;
IndTag( IndTag == 0) = [];
% Substitute independent dof's into solution vector
U.fem( IndTag, : ) = a;

end
% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.2                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.2 - Improved routinef or a more general purpose. Tested for   21/02/2013 %     
%         different elements and time-discretization scheme. but ONLY for      %
%         problems with unitary physical quantities!                           %
%   0.1 - kick-off                                                  14/02/2013 %     
%  
% ---------------------------------------------------------------------------- %


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %                 WHAT FOLLOWS IS WRONG                           % %% % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%   BAB * a(t+1) = BZB * a( t ) + BAB * find(t+1) + BQOZ * find( t ) +
%                + BQOAdir * fdir(t+1) + BQOZdir * fdir( t ) + 
%                + BQ2A * g2(t+1) + BQ2Z * g2( t ) + BQ2Adir * g2dir(t+1) +
%                + BQ2Zdir * g2dir( t ) + BQ3A * g3(t+1) + BQ3Z * g3( t ) +
%                - BAind * Tref(t+1) - BZind * Tref( t ) +
%                - BAdir * Trefdir(t+1) - BZdir * Trefdir ( t ); 
%     where: 
% 
%   BAB     = (1/de*BDBOind -   te   * BKBOind +   te   * be * BDB2ind)
%   BZB     = (1/dt*BDBOind - (1-te) * BKBOind - (1-te) * be * BDB2ind)
%   BQOA    = (Bind'*QOind) *   te
%   BQOZ    = (Bind'*QOind) * (1-te)
%   BQOAdir = (Bind'*QOdir) *   te
%   BQOZdir = (Bind'*QOdir) * (1-te)
%   BQ2A    = (Bind'*Q2ind) *   te
%   BQ2Z    = (Bind'*Q2ind) * (1-te)
%   BQ2Adir = (Bind'*Q2dir) *   te
%   BQ2Zdir = (Bind'*Q2dir) * (1-te)
%   BQ3A    = (Bind'*Q3ind) *   te
%   BQ3Z    = (Bind'*Q3ind) * (1-te)
%   BAind   = (1/dt*Bind'*DOind + Bind'*KOind *  te   + Bind'*D2ind* be *   te
%   BZind   = -1/dt*Bind'*DOind + Bind'*KOind *(1-te) + Bind'*D2ind* be * (1-te)
%   BAdir   = (1/dt*Bind'*DOdir + Bind'*KOdir *  te   + Bind'*D2dir* be *   te
%   BZdir   = -1/dt*Bind'*DOdir + Bind'*KOdir *(1-te) + Bind'*D2ind* be * (1-te)
% 
%     where: 
%                                           |                 size
%   BDBOind = Bind' * DOind * Bind          |  [ee x ind][ind x ind ][ind x ee]
%   BKBOind = Bind' * KOind * Bind          |  [ee x ind][ind x ind ][ind x ee]
%   BDB2ind = Bind' * B2ind * Bind      
%   QOind   = Bind' * QOind
%   QOdir   = Bind' * QOdir
%   Q2ind   = Bind' * Q2ind
%   Q2dir   = Bind' * Q2dir
%   Q3ind   = Bind' * Q3ind
%   BDOind  = Bind' * DOind
%   BDOdir  = Bind' * DOdir
%   BkOind  = Bind' * KOind
%   BKOdir  = Bind' * KOdir
% 
%     where:
%   DOind = Int[ Nn * Nm | Omega]
%   KOind = Int[ Nn,j * Nm,j | Omega]
%   QOind = Int[ 1/rho Nn | Omega ]
%   Q2ind = Sur[ alpha/gamma * Nn | Sigma_2]
%   Q3ind = Sur[ alpha * Nn | Sigma_3]
%   D2ind = Sur[ alpha/gamma * Nn * Nm | Sigma_2 ]
%   
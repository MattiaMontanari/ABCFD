% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SOLVES LINEAR POD SYSTEM                       %
%  /----\ |  \|    |--  |   |                                                  %
% /      \|__/ \__ |    |__/                                                   %
%                                                                              %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %


function [ POD ]=abCFD_POD_solver( U,H,R,B,D_T,D_uv,K_T,K_uv,K_Tuv,Cx,Cy,ph,tSteps,ee,var,mesh)
% Assembly the following linear system
% 
%   BAB * a(t+1) = BZB * a( t ) + B'* (- A * Tref(t+1) + Z * Tref( t ) ) +
%          + B'*( F + te *(S2(t+1) + S4(t+1)) + (1-te) * (S2( t ) + S4( t )))
% where:
%       BAB     = (1/dt*BDB +   te   * BKKB )  
%       BZB     = (1/dt*BDB - (1-te) * BKKB )
%       BDB     = B' * D * B;
%       BKB     = B' * K * B;

%% INITIALIZE VARIABLES

% Velocity tags
ta3 = ( mesh{ 2 }.ele{ 1 }.dof{ 3 }.uniqTAG );
ta4 = ( mesh{ 2 }.ele{ 1 }.dof{ 4 }.uniqTAG );
% Tetha scheme
te = ph.tetha;
% Number of time steps in current phase
nt = numel( tSteps ); 
% Numer of ALL dofs
DofNum = sum( var.fieldNDofs );
% Unknown array
Q = zeros( DofNum, nt );
% Unknown reduced array
a = zeros( ee, nt );

%% INITIALIZE ALGEBRAIC SYSTEM
%While in FEM you should have deleted matriced to take into accout the
%essential boundary conditions, with POD you have:
%   1. to solve for a homogeneous field which has homogeneous BCs'
%   2. your modes are homogeneous on some boundaries
% For these reasons you won't have to compute deleted matrices

% Damping matrix
D = D_T + D_uv;
% Stiffness matrix ( C contribute not included here)
K = K_T + K_uv + K_Tuv;
% Convective matrix
C_0 = sparse( DofNum , DofNum );
C_1 = sparse( DofNum , DofNum );
% Robin's Stiffness contribute
% KR = assemb.KRobin;
% Robin's boundary contribute
% SR = assemb.SRobin;
% Neumann's boundary contribute
% SN = assemb.SNeuma;
% Modes matrices

%% DELETE PRESSURE ENTRIES
% From solution arrays..
% U.d1( mesh{2}.ele{1}.dof{2}.uniqTAG , : ) = [];
% H.d1( mesh{2}.ele{1}.dof{2}.uniqTAG , : ) = [];
% R.d1( mesh{2}.ele{1}.dof{2}.uniqTAG , : ) = [];
% % ..and from matrices
% B( mesh{2}.ele{1}.dof{2}.uniqTAG , : ) = []; % Delet only row
% D( mesh{2}.ele{1}.dof{2}.uniqTAG , : ) = []; D( : , mesh{2}.ele{1}.dof{2}.uniqTAG  ) = [];
% K( mesh{2}.ele{1}.dof{2}.uniqTAG , : ) = []; K( : , mesh{2}.ele{1}.dof{2}.uniqTAG  ) = [];
% C( mesh{2}.ele{1}.dof{2}.uniqTAG , : ) = []; C( : , mesh{2}.ele{1}.dof{2}.uniqTAG  ) = [];

%% REDUCED SYSTEM MATRICES

BKB = B' * ( K ) * B;
BDB = B' * D * B;

%% APPLY INITIAL CONDITIONS
%Since initial conditions are known, you can start your computations from
%time step #2.
a(:,1) = BDB \ (B' * D * H.d1(:,2) ) ;

%% COMPUTE SOLUTION IN TIME
%   runs for nt+1-2 time steps because a is computed starting from the 
%   2nd physical solution. 

for id_t = 2 : 1 : nt + 1 - 2
    % Current time step
    dt = tSteps( id_t , 1 ) - tSteps( id_t-1 , 1 );
    % Compute CONVECTIVE MATRIX
    Cpod = Cx(ta3,ta3) * diag(B(ta3,:) * a(: , id_t-1 )) + Cy(ta3,ta3) * diag(B(ta4,:) * a(: , id_t-1 ));
    C_0( ta3 , ta3 ) = Cpod;
    C_0( ta4 , ta4 ) = Cpod;
    Cref_0 = Cx(ta3,ta3) * diag( R.d1(ta3 , id_t - 1 ) ) + Cy(ta3,ta3) * diag(R.d1(ta4 , id_t - 1 ));
    C_0( ta3 , ta3 ) = C_0( ta3 , ta3 ) + Cref_0;
    C_0( ta4 , ta4 ) = C_0( ta4 , ta4 ) + Cref_0;
    Cref_1 = Cx(ta3,ta3) * diag( R.d1(ta3 , id_t ) ) + Cy(ta3,ta3) * diag(R.d1(ta4 , id_t));
    C_1( ta3 , ta3 ) = Cpod + Cref_1;
    C_1( ta4 , ta4 ) = Cpod + Cref_1;

    
    % Compute reduced matrices
    BAB	= ( 1/dt.*BDB +   te   * BKB +   te   * B' * C_1 * B); 
 
    BZB = ( 1/dt.*BDB - (1-te) * BKB - (1-te) * B' * C_0 * B); 
  
    BR  = B' * ( - te * C_1 * R.d1(: , id_t) - (1-te) * C_0 * R.d1(: , id_t - 1) );
    % Right hand side
    f = BZB * a(: , id_t - 1) + BR;
    % Compute reduced vector
    a(: , id_t ) = BAB \ f ; clf; bar(a(:,id_t))
    % Recontrstuct physical solution
    Q(: , id_t ) = B * a(: , id_t );
end

POD = R.d1 + Q;
POD( : , 1 ) = U.d1( : , 1 );
% POD( ([mesh{2}.ele{1}.dof{1}.uniqTAG ; ...
%        mesh{2}.ele{1}.dof{2}.uniqTAG...
%        mesh{2}.ele{1}.dof{3}.uniqTAG ;...
%        mesh{2}.ele{1}.dof{4}.uniqTAG]),...
%       1:nt) = R.d1 + Q;
% POD( ([mesh{2}.ele{1}.dof{1}.uniqTAG ; mesh{2}.ele{1}.dof{2}.uniqTAG...
%     mesh{2}.ele{1}.dof{3}.uniqTAG ;...
%     mesh{2}.ele{1}.dof{4}.uniqTAG]),...
%       1) = U.d1( : , 1);

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  APRIL 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick off                                                  05/04/2013 %
% ---------------------------------------------------------------------------- %

% 
%     %% COMPUTE SOLUTION IN TIME
%     %   runs for nt+1-2 time steps because a is computed starting from the 
%     %   2nd physical solution. 
% 
%     for id_t = 2 : 1 : nt + 1 - 2
%         % Current time step
%         dt = tSteps( id_t , 1 ) - tSteps( id_t-1 , 1 );
% 
%         % U-velocity @ t = tn-1
%         Qu = sparse( Q( ta3, id_t - 1 ));
%         % V-velocity @ t = tn-1
%         Qv = sparse( Q( ta4, id_t - 1 ));
% 
%         % Compute CONVECTIVE MATRIX
%         C( ta3 , ta3 ) = Cx( ta3 , ta3 ) * diag( Qu ) + Cy( ta3 , ta3 ) * diag( Qv ); 
%         C( ta4 , ta4 ) = C( ta3 , ta3 );
%         BCB = B' * ( C ) * B;
%         % Compute reduced matrices
%         BAB	= ( 1/dt.*BDB +   te   * (BCB + BKB) ); 
%     %     BAB	= ( 1/dt.*BDB +   te   *    BKB    ); 
%         BZB = ( 1/dt.*BDB - (1-te) * (BCB + BKB) );
%     %     BZB = ( 1/dt.*BDB - (1-te) *    BKB    );
%         A   = ( 1/dt.* D  +   te   * ( K)  );
%         Z   = ( 1/dt.* D  - (1-te) * ( K )  );
%         % Right hand side
%         f = BZB * a(: , id_t - 1) + B' * ( - A * R.d1(: , id_t) + Z *R.d1(: , id_t - 1));
%         % Compute reduced vector
%         a(: , id_t ) = BAB \ ( f );
%         % Recontrstuct physical solution
%         Q(: , id_t ) = B * a(: , id_t );
%     end
% 
%     POD = R.d1 + Q;
%     POD( : , 1 ) = U.d1( : , 1 );
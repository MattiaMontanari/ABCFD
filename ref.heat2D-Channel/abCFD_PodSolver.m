% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SUB-ROUTINE TO SOLVE A LINEAR SYSTEM WITH A    %
%  /----\ |  \|    |--  |   |   SPECIFIC PURPOSE. WORKING ROUTINE FOR POD      %
% /      \|__/ \__ |    |__/                                                   %
%                                                                              %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ U ] = abCFD_HybSolver( assemb, dof, B, te, ph,U,init_T, DirTag) 
% output: U       - is a stucture with new fields: .fem with the solution of
%                   the linear system and, .timing list of time stepts
%                   taken by the solver
% inputs: B       - is the maxtrix which defines the size of the linear system 
%         te      - tetha coefficiente for the time-discretization scheme
%         init_T  - scalar which defines the homogeneous initial condition
%         DorTag  - DOF's tags prescribed by Diriclet's boundary conditions
% ee     - Size of linear system

function [ U ] = abCFD_PodSolver( assemb, dof, B, te, ph,U, ee)
%% Assembly the following linear system
% 
%   BAB * a(t+1) = BZB * a( t ) + B'* (- A * Tref(t+1) + Z * Tref( t ) ) +
%          + B'*( F + te *(S2(t+1) + S4(t+1)) + (1-te) * (S2( t ) + S4( t )))
% where:
%       BAB     = (1/dt*BDB +   te   * BKKB )  
%       BZB     = (1/dt*BDB - (1-te) * BKKB )
%       BDB     = B' * D * B;
%       BKB     = B' * K * B;
% and:                                                  |          size
%        B =    {modes}[ bnm , bmn ]                    |       [ind x  ee ]
%        K =    Int[ k*Nn,j * Nm,j | Omega]             |       [ind x ind ]
%        D =    Int[ rho*cp*Nn * Nm | Omega]            |       [ind x ind ]
%        F =    Int[ Nn*Q | Omega ]                     |       [ind x  1  ]
%       S2 =    Sur[ 1/pi * Nn *g2 | \Gamma_2 ]         |       [ind x  1  ]
%       S4 =    Sur[ 1/pi * Nn *g4 | \Gamma_2 ]         |       [ind x  1  ]

%% INITILIAZIE VARIABLES
% Number of time steps in current phase
nt = ph.nt; 
% Numer of ALL dofs for physical problem (in general DofNum ~= ee )
DofNum = dof.number_V1;
% Unknown array
Q = sparse( ones( DofNum, nt ) );
% Unknown delted array
a = sparse( ones( ee, nt ) );

%% ASSEMBLY DELETED MATRICES

% Coefficient deleted matrices
% B = assemb.Null' * B * assemb.Null;
D = assemb.D;
K = assemb.K;        % Pure Stiffness matrix
KR = assemb.KRobin; % Robin's Stiffness contribute
% F = assemb.Null' * 0;
SR = assemb.SRobin;
SN = assemb.SNeuma;

% Coefficient reduced-deleted matrices
BKB = B' * (K + KR) * B;
BDB = B' * D * B;
BAB	= (1/ph.step_t *BDB +   te   *    BKB    ); 
BZB = (1/ph.step_t *BDB - (1-te) *    BKB    );
A   = (1/ph.step_t * D  +   te   * (K + KR)  );
Z   = (1/ph.step_t * D  - (1-te) * (K + KR)  );

%% APPLY INITIAL CONDITIONS

a(:,1) = BDB \ (B' * D * U.hom(:,1)) ;

%% COMPUTE SOLUTION IN TIME
%   The linear system to be solved for a(t) being t [ tinit + ph.step_t , ph.t_final ]
%   runs for nt+1-1 time steps because a is computed starting from the 
%   1sh physical solution. 

% U.timing = zeros(1,nt);
for id_t = 2: 1 : nt + 1 - 1
    % Time instant computed within the id_t-th iteration
%     U.timing = U.timing + ph.step_t;

%   BAB * a(t+1) = BZB * a( t ) + B'* (- A * Tref(t+1) + Z * Tref( t ) ) +
%          + B'*( F + te *(S2(t+1) + S4(t+1)) + (1-te) * (S2( t ) + S4( t )))


    % Assemble the right hand side vector
    f = BZB * a(: , id_t - 1) + B' * ( te*(SR + SN) + (1-te)*(SR + SN) ...
                           - A * U.ref(: , id_t) + Z *U.ref(: , id_t - 1));
    % Compute reduced vector
    a(: , id_t ) = BAB \ ( f );
    % Recontrstuct physical solutio
    Q(: , id_t ) = B * a(: , id_t );
    
end

U.pod = U.ref + Q;
U.pod( : , 1 ) = U.U( : , 1);

end
 
 
% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick-off                                                  21/02/2013 %     
%  
% ---------------------------------------------------------------------------- %
 
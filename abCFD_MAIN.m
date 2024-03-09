% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \   	MAIN ROUTINE TO SOLVE A MULTIPHYSICS HEAT AND  %
%  /----\ |  \|    |--  |   |   MASS TRANSFER OF A FLOW PAST A CYLINDER.       %
% /      \|__/ \__ |    |__/    THE FULL ORDER MODEL IS SOLVED IN COMSOL AND   %
%                               WHILE A REDUCED MODEL IS BUILT VIA PROPER      %
% ORTHOGONAL DECOMPOSITION P.O.D. TECHNIQUE. THE EMPIRICAL EIGENFUNCTIONS ARE  %
% COMPUTED WITH THE SNAPSHOT METHOD. BUOYANCY EFFECT IS INCLUDED.              %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %

%%   =======================================================================  %%
%%  - - - - - - - - - - - -   START   - - - - - - - - - - - - - - - - - - - - %%
%%   =======================================================================  %%

% Starter
abCFD_starter();
clear all; 
ModelVersion = '1_0_Re100';

%%   =======================================================================  %%
%%  - - - - - - - - - - -  G E O   F I L E  - - - - - - - - - - - - - - - - - %%
%%   =======================================================================  %%
% is the complete geometry specification

%   GEOMETRY DEFINITIONS -----------------------------------------|--- UNIT -- %
L = 0.2e+3;         % Characteristic lenght: cylinder diametre    |     [m]
x = [  0  , 4*L ];  % Domains' X-coordinate bounds                |     [m]
y = [  0  , 10*L ]; % Domains' Y-coordinate bounds                |     [m]
c = [ 2*L , 2*L ];  % Cylinder centre's coordinates               |     [m]
%   LOADS
% Here should be defined the number of boundaries, edges and domain. To each
%geometry entity a proper load has to be applied.
%   PHYSICS
% Here should be defined physics are activated in COMSOL
%   MESH
% Here details about the mesh should be speficied (type,size, etc.)

%%   =======================================================================  %%
%%  - - - - - - - - - - - -  S E T T I N G S  - - - - - - - - - - - - - - - - %%
%%   =======================================================================  %%

%%   DISCRETIZATION PROPERTIES  -----------------------------------|-- UNIT -- %
ph1.COMSOL.elements = 'T1T1';  % 'Qx','Tx','Mx' with x = 1,2,3,4,5 |   [-]
ph1.COMSOL.el_grade = 5 ;       % 1 = extrafine , 9 = coarse       |   [-]
tetha = 1;          % Time discretization factor                   |   [-]

%%  F.A.Q. DEFINE ANSWERES --------------------------------------------------- %
A_PlotByComso = 'NY';       % ref. abCFD_PlotsByComsol 
A_SaveLoad = 'L';           % ref. abCFD_WalkinComsol
A_PlotModes = 'TPUV';       % ref
A_spy_FEMatrix = 'KDC';     % ref
A_spy_UNImatrix = 'KDC';    % ref

%%   =======================================================================  %%
%%  - - - - - - - - - - - - -  P H A S E S  - - - - - - - - - - - - - - - - - %%
%%   =======================================================================  %%

%% PHASE # 1 ----------------------------------------------------------------- %

%   MATERIAL PROPERTIES ------------------------------------------|--- UNIT -- %
ph1.mat.rho = 1;        % Density                                 |  [kg/m^3]
ph1.mat.mud = 1e-0;     % Dynamic viscosity                       |   [Pa*s]
ph1.mat.k = 1 ;         % Thermal conductivity                    |  [W/(m*K)]
ph1.mat.Cp = 1;         % Heat Capacity                         |  [J/(kg*K)]
ph1.mat.gamma = 1.0;    % Ration specific heats                   |     [-]
ph1.mat.thex = 1;       % Thermal expansion coefficient           |    [1/K]
%   OPERATING PARAMETERS -----------------------------------------|--- UNIT -- %
ph1.Umax = 1.0;         % Max inlet velocity                      |    [m/s]
ph1.Tmax = 1;           % Bulk (inlet) temperature                |     [K]
ph1.Tinit = 1;          % Initial temperature                     |     [K]
ph1.Tsur = 1;           % Estimated maximum temperature           |     [K]
ph1.g = 9.81;           % Gravity                                 |   [m/s^2]
%   DISCRETIZATION PROPERTIES  -----------------------------------|--- UNIT -- %
ph1.step_t = 50;     % Time step                                   |     [s]
ph1.t_final = 10000;% Computation final time                      |     [s]
ph1.tetha = 1;      % Time discretization factor                  |     [1]

%   BOUNDARY CONDITIONS ------------------------------------------------------ %
% Heat transfer
ph1.ht.Diriclet     =  2 ;      % Edge on which applies the Diriclet BC
ph1.ht.DIR.iwrite   =  [ num2str( ph1.Tinit ) , '*' ];
ph1.ht.Robin        =  [] ;     % Edge on which applies the Robin BC
ph1.ht.ROB.iwrite   =  '0*';
ph1.ht.Neumann      =  5 : 8;   % Edges on which applies the Neumann BC
ph1.ht.NEU.iwrite   =  '1*';
ph1.ht.Q            =  '0';
ph1.ht.init         = ph1.Tinit;
% Fluid flow
ph1.fd.Diriclet     =  2 ;      % Edge on which applies the Diriclet BC
ph1.fd.DIR.iwrite   =  [ num2str( ph1.Umax ) , '*4*s*(1-s)*' ];
ph1.fd.Robin        =  [] ;     % Edge on which applies the Robin BC
ph1.fd.ROB.iwrite   =  '1*';
ph1.fd.Neumann      =  3 ;      % Edge on which applies the Neumann BC
ph1.fd.NEU.iwrite   =  '1*';
                      % The following source term isn't complete
ph1.fd.Q            = [ num2str(-ph1.g*ph1.mat.thex*ph1.mat.rho),...
                                        '*0*( 0 -',num2str( ph1.Tinit ),')'];
                      % the following init.condit. must use Stoke's flow solution
ph1.fd.init         =  0;
ph1.fd.Pressure     =  3;       % Edge on which applies null pressure

%% PHASE # 2 ----------------------------------------------------------------- %
% Actually just a change of velocity should be included.

%% AUXILIAR PHASE  ----------------------------------------------------------- %

% MATERIAL PROPERTIES
phAUX.mat.rho = ph1.mat.rho;    % Density                           |  [kg/m^3]
phAUX.mat.mud = ph1.mat.mud;    % Dynamic viscosity                 |   [Pa*s]
phAUX.mat.k = ph1.mat.k ;       % Thermal conductivity              |  [W/(m*K)]
phAUX.mat.Cp = ph1.mat.Cp;      % Heat Capacity                     |  [J/(kg*K)]
phAUX.mat.gamma = ph1.mat.gamma;% Ration specific heats             |     [1]
phAUX.mat.thex = ph1.mat.thex;  % Thermal expansion coefficient     |    [1/K]

%   DISCRETIZATION PROPERTIES  
phAUX.step_t = 0.1;
phAUX.t_final = phAUX.step_t * 2 ;
phAUX.COMSOL.elements = ph1.COMSOL.elements;
phAUX.COMSOL.el_grade = ph1.COMSOL.el_grade;

%   INITIAL CONDITIONS
phAUX.ht.init = 1;
phAUX.fd.init = 1;
%   BOUNDARY CONDITIONS
% Heat transport
    phAUX.ht.Diriclet = 2 ;     phAUX.ht.DIR.iwrite = '0';
    phAUX.ht.Robin = [];        phAUX.ht.ROB.iwrite = '0';
    phAUX.ht.Neumann = 1 ;      phAUX.ht.NEU.iwrite = '1';
% MAss transport
    phAUX.fd.Diriclet = 2 ;     phAUX.fd.DIR.iwrite = '0';
    phAUX.fd.Neumann = 1:8;     phAUX.fd.NEU.iwrite = '0';
    phAUX.fd.Pressure = 1;      phAUX.fd.ROB.iwrite = '0';


%%   =======================================================================  %%
%% - - - - - - - - - - -  S O L V E   F E M  - - - - - - - - - - - - - - - -  %%
%%   =======================================================================  %%

%% P H A S E  # 1
% Compute dimensionless parameters
[ nuv ] = abCFD_pi( L, ph1);
% Initialize bounadry conditions
[ ph1.ht ] = abCFD_BCgeneralPhysics( ph1.ht , ph1.step_t );
[ ph1.fd ] = abCFD_BCgeneralPhysics( ph1.fd , ph1.step_t );

% Inspect elements
[el_type , T_order, V_order, P_order ] = abCFD_elements( ph1.COMSOL.elements );
% Save inputs
save( [char(pwd),'/Database/geo.mat'],'x','L','y','c')

% Set up and solve a model in Comsol, or load an old model
[ model, assemb ] = abCFD_WalkinComsol( ph1, A_SaveLoad, ModelVersion );

% Import Comsol's plots
% abCFD_PlotsByComsol( model , A_PlotByComso )  

%%   =======================================================================  %%
%% - - - - - - - - - - -  E X T R A C T   D A T A  - - - - - - - - - - - - -  %%
%%   =======================================================================  %%

%% EXTRACT DATA ON GEOMETRY, MESH and EXTENDED MESH -------------------------- %
[ SmeshInfo, var, mesh ] = abCFD_GeoXmeshInfo( model );

%% EXTRACT REFERENCT SOLUTION ------------------------------------------------ %

[ R ] = abCFD_getU( model, 'Stokes', mesh );

%% EXTRACT FULL MODEL SOLUTION ----------------------------------------------- %

[ U ] = abCFD_getU( model, 'NavStok', mesh );

%%   =======================================================================  %%
%% - - - - - -  R E D U C E D   O R D E R   M O D E L  - - - - - - - - - - -  %%
%%   =======================================================================  %%

eeMin = input('Define minimum number of modes:  ');
eeMax = input('Define maximum number of modes:  ');
abCFD_ReducedOrderModel( var, mesh, R, U, [eeMin eeMax] , model, ph1 );

abCFD_ReducedOrderModel( var, mesh, R, U, [7 9] , model, ph1, [22000 ,24930], [23000, 1 ,1])

test_ReducedOrderModel( var, mesh, R, U, [7 9] , model, ph1, [22000 ,24930], [23000, 1 ,1])
test_ReducedOrderModel( var, mesh, R, U, [4 35] , model, ph1)
abCFD_ReducedOrderModel( var, mesh, R, U, [30 55 ] , model, ph1)

% The Reduced Order Model has to fulfil the canonical form:
%   
% dX_k/dt = Q_klm * X_l * X_m + K_kl * X_l + P_k + S_kl * Y_k + V_klm X_l * Y_m 
%   |        |                   |            |     |            |     |     |
%   |        |                   V            |     |            |     |     V
%   |        |  2nd order STIFFNESS operator  |     |            |     | variables
%   |        |                                |     |            |     |    set
%   |        V                                |     |            |     V    
%   |   3rd order QUADRATIC operator          |     |            |  Set of    
%   |                                         |     |            | principal   
%   |                                         |     |            | variables     
%   V                                         |     |            V
%  First order time derivative of leading     |     |    3rd order coupling
%  variables set                              |     |    operator
%                                             |     V         
%                                             | 2nd order coupling operator
%                                             V
%                             1st order operator of explicilty known factors 
% 
         
%%   =======================================================================  %%
%% - - - - - - - - - - - - - -  P L O T S  - - - - - - - - - - - - - - - - -  %%
%%   =======================================================================  %%
%     % Plot velocity modes
    abCFD_plot_3DUinTime( model, mesh, var, 2, 1, [3,4], cat(3,B,B) ,1 )
    abCFD_plot_3DUinTime( model, mesh, var, 2, 1, [3,4], cat(3,B,B) ,2 )
    abCFD_plot_3DUinTime( model, mesh, var, 2, 1, [3,4], cat(3,B,B) ,3 )
%     % Plot temperature modes
%     abCFD_plot_3DUinTime( model, mesh, var, 2, 1, 1, cat(2,B) ,1)
%     abCFD_plot_3DUinTime( model, mesh, var, 2, 1, 1, cat(2,B) ,2)
%     abCFD_plot_3DUinTime( model, mesh, var, 2, 1, 1, cat(2,B) ,3)

% SOLUTION IN TIME
    % temperature
    figure;
abCFD_plot_3DUinTime( model, mesh, var, 2, 1, 1, U , U.solinfo.solvals )
    % Velocity
    figure;
abCFD_plot_3DUinTime( model, mesh, var, 2, 1, [3,4], cat(3,U.d1,U.d1) , U.solinfo.solvals )
% NODES' TAGS ON MESH
abCFD_plot_MeshWithTags( SmeshInfo, 1:4 , mesh, var );
% SHAPE FUNCTIONS - requires polyfitn.m
abCFD_plot_ShpGlobalElem( mesh ,var, 3, 2, 1, 'tri', T_order)
% SPY MATRICES
abCFD_plot_MatrixPattern( TrUnsFulStok_DKuv, var, mesh, A_spy_FEMatrix )


% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.5                                 date:   MAY  2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.5 - Rewritten.                                                04/05/2013 %
%   0.4 - All the subroutines, up to abCFD_GeoXmeshInfo incluede,   27/03/2013 %
%           have been reviewed and checked                                     %
%   0.3 - Added section 'extract data'                              22/03/2013 %
%   0.2 - Importan changes in required inputs                       06/03/2013 %
%   0.1 - kick-off                                                  27/02/2013 %
% ---------------------------------------------------------------------------- %
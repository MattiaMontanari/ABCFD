% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \   	MAIN ROUTINE TO SOLVE A MULTIPHYSICS HEAT AND  %
%  /----\ |  \|    |--  |   |   MASS TRANSFER WITH COMSOL AND BY MEANS OF      %
% /      \|__/ \__ |    |__/    THE PROPER ORTHOGONAL DECOMPOSITION P.O.D.     %
%                               TECHNIQUE. THE EMPIRICAL EIGENFUNCTIONS ARE    %
% COMPUTED WITH A DIRECT METHOD COMPUTES. THE GEOMETRY IS SAME AS VERTICAL     %
% VON KARMAN VORTEX STREE TEST CASE. BODY FORCES ARE INCLUDED AS BUOYANCY      %
% EFFECT. CONVECTIVE ACTION COOLS DOWN THE CYLINDER                            %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %

%%   =======================================================================  %%
%%  - - - - - - - - - - - -   START   - - - - - - - - - - - - - - - - - - - - %%
%%   =======================================================================  %%
% Starter
abCFD_starter();
clear all; 
ModelVersion = '0_0_2';

%%   =======================================================================  %%
%%  - - - - - - - - - - -  G E O   F I L E  - - - - - - - - - - - - - - - - - %%
%%   =======================================================================  %%

%   GEOMETRY DEFINITIONS -----------------------------------------|--- UNIT -- %
L = 0.2e+3;           % Characteristic lenght                       |     [m]
x = [  0  , 4*L ];  % X-coordinate bounds                         |     [m]
y = [  0  , 10*L ]; % Y-coordinate bounds                          |     [m]
c = [ 2*L , 2*L ];  % Circle center's coord                       |     [m]
%   LOADS
% Here should be defined the number of boundaries, edges and domain. To each
%geometry entity a proper load has to be applied.
%   PHYSICS
% Here should be defined physics are activated in COMSOL
%   MESH
% Here details about the mesh should be speficied (type,size, etc.)

%%   =======================================================================  %%
%%  - - - - - - - - - - - - -  P H A S E S  - - - - - - - - - - - - - - - - - %%
%%   =======================================================================  %%

%% PHASE # 1 ----------------------------------------------------------------- %

%   MATERIAL PROPERTIES ------------------------------------------|--- UNIT -- %
ph1.mat.rho = 1;        % Density                                 |  [kg/m^3]
ph1.mat.mud = 1e-0;     % Dynamic viscosity                       |   [Pa*s]
ph1.mat.k = 1 ;         % Thermal conductivity                    |  [W/(m*K)]
ph1.mat.Cp = 500;         % Heat Capacity                         |  [J/(kg*K)]
ph1.mat.gamma = 1.0;    % Ration specific heats                   |     [1]
ph1.mat.thex = 1;       % Thermal expansion coefficient           |    [1/K]
%   OPERATING PARAMETERS -----------------------------------------|--- UNIT -- %
ph1.Umax = 1.6*2;         % Max inlet velocity                      |    [m/s]
ph1.Tmax = 320;         % Inlet temperature                       |     [K]
ph1.Tbul = 315;         % Bulk temperature                        |     [K]
ph1.Tsur = 350;         % Estimated maximum temperature           |     [K]
ph1.g = 9.81;           % Gravity                                 |   [m/s^2]
%   DISCRETIZATION PROPERTIES  -----------------------------------|--- UNIT -- %
ph1.step_t = 10%0.01;  % Time step                                   |     [s]
ph1.t_final = 2000%6;    % Computation final time                      |     [s]
ph1.tetha = 1;          % Time discretizatino factor                  |     [1]

%   BOUNDARY CONDITIONS ------------------------------------------------------ %
ph1.ht.Diriclet =  [ 2 ] ;  %#ok<NBRAK>
ph1.ht.Robin    =  [] ; 
ph1.ht.Neumann  =  [5 : 8]; %#ok<NBRAK>
ph1.fd.Diriclet =  2 ; 
ph1.fd.Robin    =  [] ; 
ph1.fd.Neumann  =  3 ;
ph1.fd.Pressure =  3;
%% PHASE # 2 ----------------------------------------------------------------- %

%%   =======================================================================  %%
%%  - - - - - - - - - - - -  S E T T I N G S  - - - - - - - - - - - - - - - - %%
%%   =======================================================================  %%

%%   DISCRETIZATION PROPERTIES  -----------------------------------|--- UNIT -- %
ph1.COMSOL.elements = 'T1T1';  % 'Qx','Tx','Mx' with x = 1,2,3,4,5           |     [1]
ph1.COMSOL.el_grade = 7;       % 1 = extrafine , 9 = extremely coarse        |     [1]
ee = 20;            % Eigenmodes                                  |     [1]


%%  F.A.Q. DEFINE ANSWERES --------------------------------------------------- %
A_PlotByComso = 'NY';       % ref. abCFD_PlotsByComsol 
A_SaveLoad = 'L';          % ref. abCFD_WalkinComsol
A_PlotModes = 'TPUV';       % ref
A_spy_FEMatrix = 'KDC';     % ref
A_spy_UNImatrix = 'KDC';    % ref

%%   =======================================================================  %%
%%  - - - - - - - - - - - -  SETTING  - - - - - - - - - - - - - - - - - - - - %%
%%   =======================================================================  %%

% Compute convective cooling coefficient
[ nuv ] = abCFD_pi( L, ph1);
% Define bounadry conditinos
[ ph1 ] = abCFD_BCtoComsol( ph1 );
% Inspect elements
[el_type , T_order, V_order, P_order ] = abCFD_elements( ph1.COMSOL.elements );
% Save inputs
save( 'geo.mat','x','L','y','c')
% save('input.mat')

%%   =======================================================================  %%
%% - - - - - - - - - - -  S O L V E   F E M  - - - - - - - - - - - - - - - -  %%
%%   =======================================================================  %%
% Set up and solve a model in Comsol, or load an old model
[ model, assemb ] = abCFD_WalkinComsol( ph1, A_SaveLoad, ModelVersion );

% Import Comsol's plots
abCFD_PlotsByComsol( model , A_PlotByComso )  

%%   =======================================================================  %%
%% - - - - - - - - - - -  E X T R A C T   D A T A  - - - - - - - - - - - - -  %%
%%   =======================================================================  %%

%% EXTRACT DATA ON GEOMETRY, MESH and EXTENDED MESH -------------------------- %
[ SmeshInfo, var, mesh ] = abCFD_GeoXmeshInfo( model );

%% EXTRACT FEM SOLUTION ------------------------------------------------------ %
% Transient Stable Full physics
% [ U ] = abCFD_eval( model, T_order, V_order, P_order, mesh );
[ U ] = abCFD_getU( model, 'TranSol', mesh );

%%   =======================================================================  %%
%% - - - - - - - - -  A U X I L I A R   M O D E L  - - - - - - - - - - - - -  %%
%%   =======================================================================  %%
 
%% INITIALIZE AUXILIAR PHASE  ------------------------------------------------ %

% MATERIAL PROPERTIES
phAUX.mat.rho = ph1.mat.rho;    % Density                           |  [kg/m^3]
phAUX.mat.mud = ph1.mat.mud;    % Dynamic viscosity                 |   [Pa*s]
phAUX.mat.k = ph1.mat.k ;       % Thermal conductivity              |  [W/(m*K)]
phAUX.mat.Cp = ph1.mat.Cp;      % Heat Capacity                     |  [J/(kg*K)]
phAUX.mat.gamma = ph1.mat.gamma;% Ration specific heats             |     [1]
phAUX.mat.thex = ph1.mat.thex;  % Thermal expansion coefficient     |    [1/K]

%   DISCRETIZATION PROPERTIES  
phAUX.step_t = 0.1;
phAUX.t_final = 1;
phAUX.COMSOL.elements = ph1.COMSOL.elements;
phAUX.COMSOL.el_grade = ph1.COMSOL.el_grade;

%   INITIAL CONDITIONS
phAUX.ht.init = 1;
phAUX.fd.init = 1;
%   BOUNDARY CONDITIONS
% Heat transport
    phAUX.ht.Diriclet = [2];     phAUX.ht.DIR.iwrite = '0';
    phAUX.ht.Robin = [];        phAUX.ht.ROB.iwrite = '0';
    phAUX.ht.Neumann = [1];      phAUX.ht.NEU.iwrite = '1';
% MAss transport
    phAUX.fd.Diriclet = [2];     phAUX.fd.DIR.iwrite = '0';
    phAUX.fd.Neumann = 1:8;     phAUX.fd.NEU.iwrite = '0';
    phAUX.fd.Pressure = 1;      phAUX.fd.ROB.iwrite = '0';

%% INIZIALITE AUXLIAR MODEL -------------------------------------------------- %

[ model ] = abCFD_mdlobj_AUX( phAUX, model);

%% ASSEMBLE:  M & V & E  ----------------------------------------------------- %

[M, V, E] = abCFD_assembl_MVE( model, assemb, mesh, sum(var.fieldNDofs));

%% ASSEMBLE:  K & D  --------------------------------------------------------- %

[ K, D] = abCFD_assembl_KD( model, assemb, mesh, U.solinfo.solvals, sum(var.fieldNDofs));
% D=ph1.mat.rho .*ph1.mat.mud;

ta1 = ( mesh{ 2 }.ele{ 1 }.dof{ 1 }.uniqTAG ); ta3 = ( mesh{ 2 }.ele{ 1 }.dof{ 3 }.uniqTAG );
ta2 = ( mesh{ 2 }.ele{ 1 }.dof{ 2 }.uniqTAG ); ta4 = ( mesh{ 2 }.ele{ 1 }.dof{4 }.uniqTAG );

%% ASSEMBLE:  C  ------------------------------------------------------------- %
% Compute test convection matrices
[ C11, C20 ] = abCFD_assembl_C( model, assemb, mesh, U.solinfo.solvals, sum(var.fieldNDofs), D);
Q=zeros(sum(var.fieldNDofs),1);
Q( [ta3, ta4],:) = 1;
[ C,Kdel,Edel] = abCFD_AssembleNavier( var, [ mesh{2}.ele{1}.dof{3}.TAG;mesh{2}.ele{1}.dof{3}.TAG], Q, mesh);
%%   =======================================================================  %%
%% - - - - - - - - -  H O M O G E N I Z A T I O N  - - - - - - - - - - - - -  %%
%%   =======================================================================  %%

%% COMPUTE REFERENCE SOLUTION ------------------------------------------------ %

[ R ] = abCFD_REF_solver( model, mesh, U.solinfo.solvals );

%% COMPUTE HOMOGENEOUS SOLUTION
[ H ] = abCFD_Uhom( U, R, mesh);

%%   =======================================================================  %%
%% - - - - - - - - - - -  S O L V E   P O D  - - - - - - - - - - - - - - - -  %%
%%   =======================================================================  %%
for ee = 5 : 10
%% Compute modes
[ H ] = abCFD_Modes( H, ee, 'eig' );
% Assemble modes
B = zeros( sum(var.fieldNDofs) ,ee );
B( mesh{2}.ele{1}.dof{1}.uniqTAG, : ) = H.T.B;
B( mesh{2}.ele{1}.dof{3}.uniqTAG, : ) = H.u.B;
B( mesh{2}.ele{1}.dof{4}.uniqTAG, : ) = H.v.B;

% B = [ H.T.B ; zeros( var.fieldNDofs( 2 ),ee); H.u.B ; H.v.B];

%% ROM SOLVER
[ POD ] = abCFD_POD_solver( U.d1,H.d1,R.d1,B,E,M,K,D.*1e-0,V,ph1,U.solinfo.solvals,ee,var,mesh); 
clf; pause(1);
abCFD_plot_3DUinTime( model, mesh, var, 2, 1, [3,4], cat(3,POD,POD) , U.solinfo.solvals )
end
%%   =======================================================================  %%
%% - - - - - - - - - - - - - -  P L O T S  - - - - - - - - - - - - - - - - -  %%
%%   =======================================================================  %%
% SOLUTION IN TIME
    % temperature
abCFD_plot_3DUinTime( model, mesh, var, 2, 1, 1, U , U.solinfo.solvals )
    % Velocity
abCFD_plot_3DUinTime( model, mesh, var, 2, 1, [3,4], cat(3,U.d1,U.d1) , U.solinfo.solvals )
% NODES' TAGS ON MESH
abCFD_plot_MeshWithTags( SmeshInfo, 1:4 , mesh, var );
% SHAPE FUNCTIONS - requires polyfitn.m
abCFD_plot_ShpGlobalElem( mesh ,var, 3, 2, 1, 'tri', T_order)
% SPY MATRICES
abCFD_plot_MatrixPattern( TrUnsFul_MVE, var, mesh, A_spy_FEMatrix )
% whole linear system's matrices pattern
%Damping
figure(); hold on;
spy(E([ta1;ta2;ta3;ta4],[ta1;ta2;ta3;ta4]))
spy(M([ta1;ta2;ta3;ta4],[ta1;ta2;ta3;ta4]),'g')
%Stiffness
figure(); hold on;
spy(K([ta1;ta2;ta3;ta4],[ta1;ta2;ta3;ta4]),'c')
spy(D([ta1;ta2;ta3;ta4],[ta1;ta2;ta3;ta4]),'r')
spy(V([ta1;ta2;ta3;ta4],[ta1;ta2;ta3;ta4]),'g')
spy(C([ta1;ta2;ta3;ta4],[ta1;ta2;ta3;ta4]),'wo')
 

% PLOT TIME STEPS
% compute time steps
% figure(); hold on; grid on;
% for id_t = 2 : U.solinfo.sizesolvals
%     semilogy(1./( U.solinfo.solvals(id_t-1) -U.solinfo.solvals(id_t)))
% end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.3                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.4 - All the subroutines, up to abCFD_GeoXmeshInfo incluede,   27/03/2013 %
%           have been reviewed and checked                                     %
%   0.3 - Added section 'extract data'                              22/03/2013 %
%   0.2 - Importan changes in required inputs                       06/03/2013 %
%   0.1 - kick-off                                                  27/02/2013 %
%                                                                              %
% ---------------------------------------------------------------------------- %


 
% % Display relative error between Comsol and abCFD's FEM solution
% disp('relative error between Comsol and abCFD''s FEM solution')
% mean(mean( sqrt(abs(U.fem.^2-U.U.^2)) ))
% disp('relative error between Comsol and abCFD''s POD solution')
% mean(mean( sqrt(abs(U.pod.^2-U.U.^2)) ))
% 
% %% POSTPROCESS RESULTS
% 
% % Identify elements on \gammaX
% abCFD_NorHT( model, ph1, dof, edg, ele, U.U', A_PlotGratT )
% % Compute heat fluxes
% 
% % Tangengial flux to \Gamma4
% abCFD_TanHT( U.pod, ph1, dof, edg, edg.gamma4, 'nY' );
% 
% 
% %% PLOT RESULTS
% % Evenutally plot element in local coordinates
% abCFD_PlotElement( A_PlotElem , ele.localCoor )
% 
% % Eventually plot Simple Mesh
% abCFD_PlotSplMesh( A_PlotSplMesh , SmeshInfo.data.vertex' , ...
%     double(model.mesh('mesh1').getElem( dof.types_V1{2} ) + 1)', elements  )
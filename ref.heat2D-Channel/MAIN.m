% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    % 
%   /  \  |__  /   |    |  \   	MAIN ROUTINE TO SOLVE A HEAT TRANSIENT LINEAR  %
%  /----\ |  \|    |--  |   |   PROBLEM WITH COMSOL AND Proper Orthogonal      %
% /      \|__/ \__ |    |__/    Decomposition METHOD. ROBIN,NEUMANN & DIRICLET %
%                               BOUNDARY CONDITIONS ARE APPLIED. HEAT          %
% SOURCE IS INCLUDED. A DIRECT METHOD IS USED TO COMPUTE THE HOMOGENEOUS       %
% EIGENFUNCTIONS, WHICH ARE OBTAINED WITH THE BETA FIELD. GEOMETRY IS SAME     %
% AS HORIZONTAL VON KARMAN VORTEX STREE GEOMETRY.                              %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %

clear all;  close all; colordef black; %Sclc;    

%% INPUT
ModelVersion = '0_2_0';
%   PHYSICAL PROBLEM SET UP -------------------------------------------------- %
rho = 1;            % Material density      [kg/m^3]
k = 1 ;             % Thermal conductivity  [W/(m*K)]
Cp = 1;             % Heat Capacity         [J/(kg*K)]
ph1.BC.Diriclet =  1 ; % Geo's edge tag with Diriclet BC            %#ok<NBRAK>
ph1.BC.Robin    =  2 ; % Geo's edge tag with Robin BC               %#ok<NBRAK>
ph1.BC.Neumann  =  4  ; % Geo's edge tag with Neumann BC             %#ok<NBRAK>
%   DISCRETIZATION PROPERTIES  ----------------------------------------------- %
elements = 'T1';    % options: 'Qx','Tx','Mx' with x = 1,2,3,4,5            [-]
el_grade = 4;       % with 1 = extrafine  and  9 = extremely coarse         [-]
ee = 2;             % Number of considered eigenmodes                       [-]
tetha = 1;          % time discretizatino factor fot tetha-family scheme    [-]
init_T = 1;         % homogeneous initial temperature                       [K]
ph1.step_t = 0.1;   % Time step                                             [s]
ph1.t_final = 1;    % Computation final time                                [s]
% %   BOUNDARY CONDITIONS ---------------------------------------------------- %
[ ph1 , el_type , el_order ] = abCFD_BCtoComsol( ph1, elements );

%%  F.A.Q. DEFINE ANSWERES
A_PlotByComso = 'Y';   	% See abCFD_PlotsByComsol 
A_BuondaryCondCheck = 'Y';  % See abCFD_PlotsByComsol 
A_PlotElem =    'nY';       % See abCFD_PlotElement
A_PlotSplMesh = 'Y';        % See abCFD_PlotSplMesh
A_PlotExtFEM = 'VF';        % See abCFD_PlotExtMesh
A_PlotExtPOD = 'nY';        % See abCFD_PlotExtMesh
A_PlotExtCOM = 'nY';        % See abCFD_PlotExtMesh
A_PlotModes =   'nY';       % Plot modes
A_PlotThom =    'nY';       % Plot homogeneous solution
A_Save =        'nY';       % Save or not the model being solved
A_spyMatrix =   'nY';       % Spy coeff. matrices
A_PlotGratT =   'nY';       % Plot gradient of particular scalar field
%   Save inputs
save( 'input.mat' )

%% COMSOL COMPUTATIONS ------------------------------------------------------- %
% Load or Run new model
A_load = input('Load or solve new model? [ {L} / Solve ] ', 's');
if strcmp( A_load, 'Solve') == 1
    [ model, assemb ] = abCFD_modelobject( ph1, A_BuondaryCondCheck, A_Save );
    save('input_assemb.mat', 'assemb');
else
    disp(['LOAD MODEL: ',cd,'heat_mode14_',ModelVersion,'.mph']);
    [ model , PathModel ] = mphload( ['heat_mode14_',ModelVersion,'.mph'] );    
    load('input_assemb.mat', 'assemb');
end

% Plot from Comsol
abCFD_PlotsByComsol( model , A_PlotByComso )

%% EXTRACT DATA ON GEOMETRY, 'geometric' MESH and EXTENDED MESH -------------- %
[ SmeshInfo, dof, ele , edg ] = abCFD_GeoXmeshInfo( model );

%% EXTRACT COEFFICENT AND SOLUTIONS ARRAYS
% Coefficient matrices
[ assemb ,FC_spy ] = abCFD_MatrixAssembly(assemb, model , dof, A_spyMatrix);

% Temperature field
solinfo = mphsolinfo(model,'nu','on');
ph1.dof = solinfo.size;
ph1.nt =  solinfo.sizesolvals;
U.U = mphgetu( model ,'Solnum', 1 : ph1.nt );
U.u = mpheval( model , {'T'}, 'solnum', [1:ph1.nt] , 'smooth','none',...
        'refine' , num2str( el_order ) ,'matherr' , 'on');
    
%% ASSEMBLE BOUNDARY CONDITIONS' COEFFICENT ARRAYS
[ assemb ] = abCFD_AssembleAllBCs( edg, dof, ele, ph1 ,assemb, 'n');
%% ASSEMBLE source contribute
% [ assemb.SourcContribute ] = abCFD_AssembleEleCntrb( edg.gamma2, dof, ele);

%% COMPUTE REFERENCE SOLUTION AND HOMOGENEOUS MODES
% Compute reference solution
[beta , U.ref, U.hom] = abCFD_HomogenizeSol( model, A_PlotThom, dof, ele, ...
            assemb.K - assemb.KRobin , U.U, elements, 3);
% Compute modes
[ eigval , B ] = abCFD_ComputeModes( U.hom,model,ee,dof, ele,elements,'eig',A_PlotModes);

%% REDUCED ORDER MODEL - ASSEMBLY LINEAR SYSTEM              ref. notes

% Solve FEM
[ U ] = abCFD_FemSolver( assemb, dof,  tetha, ...
                            ph1,U,init_T,dof.gamma1.DofTags);
% Solve POD
[ U ] = abCFD_PodSolver( assemb, dof, B, tetha, ph1, U, ee);
% Display relative error between Comsol and abCFD's FEM solution
disp('relative error between Comsol and abCFD''s FEM solution')
mean(mean( sqrt(abs(U.fem.^2-U.U.^2)) ))
disp('relative error between Comsol and abCFD''s POD solution')
mean(mean( sqrt(abs(U.pod.^2-U.U.^2)) ))

%% POSTPROCESS RESULTS

% Identify elements on \gammaX
abCFD_NorHT( model, ph1, dof, edg, ele, U.U', A_PlotGratT )
% Compute heat fluxes

% Tangengial flux to \Gamma4
abCFD_TanHT( U.pod, ph1, dof, edg, edg.gamma4, 'nY' );


%% PLOT RESULTS
% Evenutally plot element in local coordinates
abCFD_PlotElement( A_PlotElem , ele.localCoor )

% Eventually plot Simple Mesh
abCFD_PlotSplMesh( A_PlotSplMesh , SmeshInfo.data.vertex' , ...
    double(model.mesh('mesh1').getElem( dof.types_V1{2} ) + 1)', elements  )

% Eventually Plot U.U in time
abCFD_PlotExtMesh(model, A_PlotExtCOM, dof.dofCoord', ele.DOFToNode',...
                                U.U, elements)
% Eventually Plot U.fem in time
abCFD_PlotExtMesh(model, A_PlotExtFEM, dof.dofCoord', ele.DOFToNode',...
                                U.fem, elements)
                            
% Eventually Plot U.pod in time
abCFD_PlotExtMesh(model, A_PlotExtPOD, dof.dofCoord', ele.DOFToNode',...
                                U.pod, elements)

%     %% EVALUATION SOLUTION FIELDS
% 
%     % %function MPHEVAL
%     NodalEvaluated = mpheval( model , {'T','ht.tfluxx','ht.tfluxy','x','y' } ,...
%         'smooth', 'none' , 'solnum' , 'all', 'refine' , num2str( el_order ) ,...
%         'matherr' , 'on') ;
% 
% NodalEvaluated = mpheval( model , {'T', 'ht.tfluxx','ht.tfluxy' } ,...
%         'solnum' , [1:solinfo.sizesolvals] , 'refine' , num2str( el_order ) ,...
%         'matherr' , 'on');

% ht.ntflux

% 
%  NodalEvaluated = mpheval( model , {'T', 'ht.tfluxx','ht.tfluxy' ,'ht.ntflux'} ,...
%         'solnum' , [1:solinfo.sizesolvals] , 'smooth','none','refine' , num2str( el_order ) ,...
%         'matherr' , 'on','edim',1,'selection',2);



% The six kinds of fields included concern:
%  << expr >>    [ ] : names of the evaluated expressinos
%  << d1,..dn >>     : Solution fields.columns correspond to node poing coordinates in columns in 'p'
%                [x] :      .d1 --> 'T' temperature field
%                [ ] :      .d2 --> 'ht.tfluxx'
%                [ ] :      .d3 --> 'ht.tfluxy
%                [ ] :      .d4 --> 'x'
%                [ ] :      .d5 --> 'y'
%  <<  p  >>     [ ] : [ X-Ycood.   x     ??        ]
%                      max(.p), min(.p) match the geometry verteces
%  <<  t  >>     [ ] : [    3       x   Simpl.el.Nu ]

%  <<  ve  >>    [ ] : [  ??        x       1     ]
%                      max(.ve) = Free.el.Nu   |     min(.ve) = 1
%  << unit >>    [ ] : list of unit for each expression.


% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.2                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.2 - THE CODE IS STILL UNDER CONSTRUCTION AND SINCE TODAY IS   08/02/2013 %
%       ABOUT 140 lines. THE HOMOGENEOUS TEMPERATURE FIELD IS CALCULATED       %
%       AND PLOTTED ON DIFFERENT KINDS MESH TOPOLOGIES. SO FAR ELEMENTS T1,    %
%       T2, T3, Q1, Q2 HAVE BEEN TESTED. LINEAR ELEMENTS ARE THE BEST          %
%       NEXT IMPROVEMENTS SHOULD IMPROVE THE DATA MANAGEMENT. ONE MAIN         %
%       STRUCTURE SHOULD INCLUDE ALL RELEVANT DATA TO WORK ON THE EXTENDEDMESH %
%   0.1 - kick-off                                                  05/02/2013 %
%                                                                              %
% ---------------------------------------------------------------------------- %

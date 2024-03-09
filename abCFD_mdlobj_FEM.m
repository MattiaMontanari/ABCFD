% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \  	MODELOBJECT FOR 2D HEAT AND MASS TRANSFER      %
%  /----\ |  \|    |--  |   |   PROBLEM IN A VERTICAL CHANNEL WITH A CYLINDER. %
% /      \|__/ \__ |    |__/    EIGHT BOUNDARY CONDITIONS PLUS HEAT SOURCE AND %
%                               BODY FORCES ARE TO BE APPLIED                  %
%  * * * CALLS * * *                                                           %
%           i. abCFD_parameters                                                %
%          ii. abCFD_functions                                                 %
%         iii. abCFD_variables                                                 %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ model, assemb ] = abCFD_mdlobj_FEM( ph, {NewModelID} )
% INPUTS:
%   - ph           : common phase structure
%   - {NewModelID} : defines a new model ID to be used, different from the 
%                       one loaded in 'input.mat'                     

function [ model, assemb ] = abCFD_mdlobj_FEM( ph,ModelVersion, ModelName, varargin )

%% INITIALIZE ---------------------------------------------------------------- %
% Load input
load( [char(pwd),'/Database/geo.mat'] );
% Elements
[el_type , T_order, V_order, P_order ] = abCFD_elements( ph.COMSOL.elements );

% check input
if numel( varargin ) ~= 0
   ModelVersion = varargin{ 1 };
end

% Initizilize the matrices to be extracted
% matrices = ['uscale';'K     ';'L     ';'M     ';'N     ';'D     ';'E     ';...
%             'NF    ';'NP    ';'MP    ';'MLB   ';'MUB   ';'Kc    ';'Lc    ';...
%             'Dc    ';'Ec    ';'Null  ';'Nullf ';'ud    '];
assemb.matrices = {'K','N','D','E','NF','NP','Kc','Dc','Ec','Null','Nullf'};
assemb.vector = {'L','M','MP','MLB','MUB','Lc','ud','uscale'};


import com.comsol.model.*
import com.comsol.model.util.*
% Model generalities
model = ModelUtil.create( ['Model',ModelVersion] );
model.modelPath( pwd ) ;

ModelUtil.showProgress( true );
model.version( ModelVersion );
model.author( 'Montanari Mattia') ;

% Create model, geometry, mesh and physic
model.modelNode.create( ModelName );
geom = model.geom.create( ['geom',ModelName], 2);
mesh = model.mesh.create( ['mesh',ModelName], ['geom',ModelName] ); 

% Create physics
singlephase = model.physics.create('spf_FEM', 'LaminarFlow', ['geom',ModelName] );
heatransfer = model.physics.create('ht_FEM', 'HeatTransfer', ['geom',ModelName] );

% Create stadionary study
Stdy_Steady = model.study.create('Study_Std');
% Create Stationary step and activate physics
StatStep = Stdy_Steady.feature.create('stat', 'Stationary');
StatStep.activate('spf_FEM', true);
StatStep.activate('ht_FEM', true);

% Create transient study
Stdy_Transi = model.study.create('Study_Trn');
% Create Transient step and activate physics
TranStep = Stdy_Transi.feature.create('time', 'Transient');
TranStep.activate('spf_FEM', true);
TranStep.activate('ht_FEM', true);

% Set interpolation function order
singlephase.prop('ShapeProperty').set('order_fluid', 1, num2str( V_order ));


%% Single fluid flow phase --------------------------------------------------- %
% Apply incompressibility
singlephase.prop('CompressibilityProperty').set('Compressibility', 1, 'Incompressible');
    
% Set up time ranges
TranStep.set('tlist', ...
    ['range(0,200,5000)',',','range(5000,', num2str(ph.step_t),',',num2str(ph.t_final),')' ]);
  
% transtudy.set('tlist', ...
%     ['range(0,', num2str(ph.step_t),',',num2str(ph.t_final)...
%     ,') range(',num2str( ph.t_final+ph.step_t),',',...
%       num2str(ph.step_t/10),',', num2str(3*ph.t_final),')']);

%% Parameters' definitions --------------------------------------------------- %
abCFD_parameters( model, ph );

%% Functions' definitions ---------------------------------------------------- %
% abCFD_functions( model, ModelName );

%% Variables definition ------------------------------------------------------ %
abCFD_variables( model, ModelName);

%% Geo file ------------------------------------------------------------------ %
% Rectangle
geom.feature.create('r1', 'Rectangle');
geom.feature('r1').set('type', 'solid');
geom.feature('r1').set('base', 'corner');
geom.feature('r1').set('pos', [x(1) y(1)] );
geom.feature('r1').set('lx', num2str( x(2)-x(1) ) );
geom.feature('r1').set('ly', num2str( y(2)-y(1) ));
% Circle
geom.feature.create('c1', 'Circle');
geom.feature('c1').set('type', 'solid');
geom.feature('c1').set('base', 'center');
geom.feature('c1').set('pos', c );
geom.feature('c1').set('r', num2str( L/2 ) );
% Difference
geom.feature.create('dif1', 'Difference');
geom.feature('dif1').selection('input').set({'r1'});
geom.feature('dif1').selection('input2').set({'c1'});
geom.runAll;
geom.run;

%% Selections ---------------------------------------------------------------- %
 

%% Material ------------------------------------------------------------------ %
material = model.material.create('mat1');
material.propertyGroup('def').set('thermalconductivity', ph.mat.k );
material.propertyGroup('def').set('density', ph.mat.rho );
material.propertyGroup('def').set('heatcapacity', ph.mat.Cp );
material.propertyGroup('def').set('dynamicviscosity', ph.mat.mud );
material.propertyGroup('def').set('ratioofspecificheat', ph.mat.gamma );


%%  Heat Transfer in solids -------------------------------------------------- %

heatransfer.feature.create('fluid1', 'FluidHeatTransferModel', 2);
heatransfer.feature('fluid1').selection.set([1]);
heatransfer.feature('fluid1').set('minput_velocity_src', 1, ['root.',ModelName,'.u']);

% Heat transfer in fluids 
heatransfer.prop('ShapeProperty').set('order_temperature', 1, num2str( T_order ));

%% Apply boundary contidions ------------------------------------------------- %

% Initial boundary conditions 
heatransfer.feature('init1').set('T', 1, 'initial_T');
singlephase.feature('init1').set('u', [ ph.fd.init , ph.fd.init , 0 ]);
singlephase.feature('init1').set('p', 1, 'rho*g_const*(L-y)');

% HEAT SOURCE 
heatransfer.feature.create('heatsource', 'HeatSource', 2);
heatransfer.feature('heatsource').selection.set( [] );
heatransfer.feature('heatsource').set('Q', 1, ph.ht.Q  );

% BODY FORCE
singlephase.feature.create('buoyancy', 'VolumeForce', 2);
singlephase.feature('buoyancy').selection.set( 1 );
singlephase.feature('buoyancy').set('F', {'0' ph.fd.Q '0'});

% INLET    - SINGLE PHASE FLOW
singlephase.feature.create('inl1', 'Inlet', 1);
singlephase.feature('inl1').selection.set( ph.fd.Diriclet );
singlephase.feature('inl1').set('U0in', 1, ...
    [ ph.fd.DIR.iwrite,'1'] );%'HT_DiricletRamp(t[1/s])'] );
% OUTLET   - SINGLE PHASE FLOW
singlephase.feature.create('out1', 'Outlet', 1);
singlephase.feature('out1').selection.set( ph.fd.Neumann );
singlephase.feature('out1').set('BoundaryCondition', 1, 'NormalStress');
% singlephase.feature('out1').set('BoundaryCondition', 1, ...
%     'PressureNoViscousStress');

% DIRICLET - HEAT TRANSFER
heatransfer.feature.create('Diricl', 'TemperatureBoundary', 1);
heatransfer.feature('Diricl').selection.set([ ph.ht.Diriclet ]);  
heatransfer.feature('Diricl').name('Temperature');
heatransfer.feature('Diricl').set(...
    'T0', 1, [ ph.ht.DIR.iwrite , '1'] );%'HT_DiricletRamp(t[1/s])+initial_T'] );

% ROBIN    - HEAT TRANSFER
heatransfer.feature.create('Robin', 'ConvectiveCooling', 1);
heatransfer.feature('Robin').selection.set( ph.ht.Robin );  
heatransfer.feature('Robin').name('Convection');
heatransfer.feature('Robin').set('h', 1, 'HT_ROB_h');
heatransfer.feature('Robin').set('Text', 1, ...
    [ ph.ht.ROB.iwrite ,'1'] );% 'HT_RobinRamp(t[1/s])'] );

% NEUMANN  - HEAT TRANSFER
heatransfer.feature.create('Neumann', 'HeatFluxBoundary', 1);
heatransfer.feature('Neumann').selection.set([ ph.ht.Neumann ]);  
heatransfer.feature('Neumann').name('HeatFlux');
heatransfer.feature('Neumann').set('q0', 1, ...
	[ ph.ht.NEU.iwrite ,'1'] );% 'HT_NeumannRamp(t[1/s])'] );
% heatransfer.feature('Neumann').set('HeatFluxType', 1, 'GeneralInwardHeatFlux');


%% Mesh setting -------------------------------------------------------------- %
mesh.automatic(false);
mesh.feature.remove('size1');
% mesh.feature.remove('cr1');
mesh.feature.remove('bl1');
mesh.feature.remove('ftri1');
mesh.feature.remove('ftri2');

if 1 == strcmp( el_type , 'Free_Quad' )
% *---*  Free quadrilateral mesh  
    mesh.feature.create('msh', 'FreeQuad');
    mesh.feature('size').set('hauto', num2str( ph.COMSOL.el_grade ));
    
elseif 1 == strcmp( el_type , 'Free_Tria' )
% *---*  Free triangles
    mesh.feature.create('msh', 'FreeTri');
    
    % TRIVIAL MESH
    trivialmesh = 'false';
    if 1 == strcmp( trivialmesh , 'active') ;
        
        mesh.feature('msh').feature.create('dis1', 'Distribution');
        mesh.feature('msh').feature('dis1').selection.all;
        mesh.feature('msh').feature('dis1').set('numelem', '1');
        mesh.feature('size').set('custom', 'on');
        mesh.feature('size').set('hmax', '5');
        mesh.feature('size').set('hmin', '5');
        mesh.run;
    else
        mesh.feature('size').set('hauto', num2str( ph.COMSOL.el_grade ));
    end
%     mesh.feature('size').set('hgrad', '1');


elseif 1 == strcmp( el_type , 'Map_Quad' );
% *---*  Structured grid. This requires the domains to be joint
    geom.feature('c1').active(false);
    geom.feature('dif1').active(false);
    geom.run;
    % Apply element distribution
    mesh.feature.create('msh', 'Map');
    mesh.feature('msh').feature.create('dis1', 'Distribution');
    mesh.feature('msh').feature('dis1').set('numelem', '2');
    mesh.feature('msh').feature('dis1').selection.set([3 2]);
    mesh.feature('msh').feature.create('dis2', 'Distribution');
    mesh.feature('msh').feature('dis2').set('numelem', '2');
    mesh.feature('msh').feature('dis2').selection.set([1 4]);
    mesh.feature('size').set('hauto', num2str( el_grade ));
    
end
 
mesh.run;

%% STUDY --------------------------------------------------------------------- %
    
%% STUDY CONFIGURATION ------------------------------------------------------ %

% Initialize Stationary solver for Stokes flow
StokesStatNode = 'Stokes';
StokStat_eq = 'StokStat_eq';
StokStat_var = 'StokStat_var';
StokStat_sol = 'StokStat_sol';
StokStat_direct = 'StokStat_direct';
StokStat_couple = 'StokStat_couple';
% Create Stationary solver for Stokes flow
StokesStatNode = model.sol.create( StokesStatNode );
StokesStatNode.study('Study_Std');
StokesStatNode.feature.create( StokStat_eq , 'StudyStep');
StokesStatNode.feature( StokStat_eq ).set('study', 'Study_Std');
StokesStatNode.feature( StokStat_eq ).set('studystep', 'stat');
StokesStatNode.feature.create( StokStat_var , 'Variables');
StokesStatNode.feature( StokStat_var ).set('control', 'stat');
StokesStatNode.feature( StokStat_var ).set('scalemethod', 'none');
StokesStatSolver = StokesStatNode.feature.create( StokStat_sol , 'Stationary');
% StokesStatSolver.feature.create('seDef', 'Segregated');
StokesStatSolver.feature.create( StokStat_couple , 'FullyCoupled');
StokesStatSolver.feature( StokStat_couple ).set('initstep', 0.01);
StokesStatSolver.feature( StokStat_couple ).set('minstep', 1.0E-6);
StokesStatSolver.feature( StokStat_couple ).set('dtech', 'auto');
StokesStatSolver.feature( StokStat_couple ).set('maxiter', 50);
StokesStatSolver.feature.create( StokStat_direct , 'Direct');
StokesStatSolver.feature( StokStat_direct ).set('linsolver', 'pardiso');
StokesStatSolver.feature( StokStat_couple ).set('linsolver',  StokStat_direct );
StokesStatSolver.feature( StokStat_couple ).set('initstep', 0.01);
StokesStatSolver.feature( StokStat_couple ).set('minstep', 1.0E-6);
StokesStatSolver.feature( StokStat_couple ).set('dtech', 'auto');
StokesStatSolver.feature( StokStat_couple ).set('maxiter', 50);
StokesStatSolver.feature('aDef').set('convinfo', 'detailed');
StokesStatSolver.feature.remove('fcDef');
StokesStatSolver.feature.remove('seDef');
StokesStatNode.feature.create('StoredStokes', 'StoreSolution'); % Store solution
StokesStatNode.feature('StoredStokes').name('StoredStokes');    % Give it a name
StokesStatNode.attach('Study_Std');

% Initialize Transient solver for Navier-Stokes flow
NavStoksTranNode = 'NavStok';
NavSto_eq = 'NavStok_eq';
NavSto_var = 'NavStok_var';
NavSto_sol = 'NavStok_sol';
NavSto_direct = 'NavStok_direct';
NavSto_couple = 'NavStok_couple';
% Create Transient solver for Navier-Stokes flow
NavStksTranNode = model.sol.create( NavStoksTranNode );
NavStksTranNode.study('Study_Trn');
NavStksTranNode.feature.create( NavSto_eq , 'StudyStep');
NavStksTranNode.feature( NavSto_eq ).set('study', 'Study_Trn');
NavStksTranNode.feature( NavSto_eq ).set('studystep', 'time');
NavStksTranNode.feature.create( NavSto_var , 'Variables');
NavStksTranNode.feature( NavSto_var ).set('control', 'time');
NavStksTranNode.feature( NavSto_var ).set('scalemethod', 'none');
NavStksTranSolver = NavStksTranNode.feature.create( NavSto_sol , 'Time');
NavStksTranSolver.set('tlist', 'range(0,0.2,3.4) range(3.5,0.02,7)');
NavStksTranSolver.set('plot', 'off');
NavStksTranSolver.set('plotfreq', 'tout');
NavStksTranSolver.set('probesel', 'all');
NavStksTranSolver.set('probes', {});
NavStksTranSolver.set('probefreq', 'tsteps');
NavStksTranSolver.set('atolglobalmethod', 'scaled');
NavStksTranSolver.set('atolglobal', 0.0010);
NavStksTranSolver.set('estrat', 'exclude');
NavStksTranSolver.set('maxorder', 2);% Time discret-maximum order
NavStksTranSolver.set('control', 'time');
NavStksTranSolver.set('eventtol', '0.001');
NavStksTranSolver.set('tout', 'tsteps');           % all timesteps taken by solver
NavStksTranSolver.set('tstepsbdf', 'intermediate');% Take intermediate time steps  
NavStksTranSolver.feature('aDef').set('convinfo', 'detailed');
NavStksTranSolver.set('atolglobalmethod', 'unscaled');
NavStksTranSolver.feature.create( NavSto_couple , 'FullyCoupled');
NavStksTranSolver.feature( NavSto_couple ).set('jtech', 'once');
NavStksTranSolver.feature( NavSto_couple ).set('maxiter', 5);
NavStksTranSolver.feature.create( NavSto_direct , 'Direct');
NavStksTranSolver.feature( NavSto_direct ).set('linsolver', 'pardiso');
NavStksTranSolver.feature( NavSto_couple ).set('linsolver',  NavSto_direct );
NavStksTranSolver.feature( NavSto_couple ).set('jtech', 'once');
NavStksTranSolver.feature( NavSto_couple ).set('maxiter', 6);
NavStksTranSolver.feature.remove('fcDef');
NavStksTranNode.detach;     % Do not create default plot 
NavStksTranNode.feature( NavSto_var ).set('initmethod', 'sol'); % Select the previous
NavStksTranNode.feature( NavSto_var ).set('initsol', 'Stokes');
NavStksTranNode.feature( NavSto_var ).set('initsoluse', 'sol1');%  computation as initial condition
NavStksTranNode.feature.create('StoredStokes', 'StoreSolution'); % Store solution
NavStksTranNode.feature('StoredStokes').name('StoredNavStks');  % Give it a name
NavStksTranNode.attach('Study_Trn');

%% ASSEMBLY NODES ------------------------------------------------------------ %
% Create assembly feature for transient solver
assembl = StokesStatNode.feature.create('asmbl','Assemble');
for id_m = 1 : size( assemb.matrices , 2 )
    assembl.set(  assemb.matrices{id_m}    , 'on'  );
end
for id_v = 1 : size(assemb.vector , 2 )
    assembl.set(  assemb.vector{id_v}    , 'on'  );
end
% Avoid scaling in assembly
assembl.feature('aDef').set('rowscale', 'off');

% Create assembly feature for transient solver
assembl = NavStksTranNode.feature.create('asmbl','Assemble');
for id_m = 1 : size( assemb.matrices , 2 )
    assembl.set(  assemb.matrices{id_m}    , 'on'  );
end
for id_v = 1 : size(assemb.vector , 2 )
    assembl.set(  assemb.vector{id_v}    , 'on'  );
end
model.param.set('timestep','1');
model.param.set('t', '0.1');
% Avoid scaling in assembly
assembl.feature('aDef').set('rowscale', 'off');


fprintf('Modelobject setup complete after %f sec. \n Start Comsol computation at %s \n',...
toc, datestr(now,15) )
tii = toc;


%% AD HOC     STUDY AND SOLVER SETTING

% NEGLECT INERTIAL TERM - STOKES
model.physics('spf_FEM').prop('StokesFlowProp').set('StokesFlowProp', 1, '1');
% NEGLECT BUOYANCY
model.physics('spf_FEM').feature('buoyancy').active(false);
% STATIONARY
model.physics('spf_FEM').prop('EquationForm').set('form', 1, 'Stationary');
model.physics('ht_FEM').prop('EquationForm').set('form', 1, 'Stationary');
% EVENTUALLY RUN STOKES PROBLEM
model.sol('Stokes').runAll;
disp('Stokes Flow Computed')

tii = toc;
fprintf('\n - + - + Solving the full model - + - + \n')
% CONSIDER INERTIAL TERM - NAVIER STOKES
model.physics('spf_FEM').prop('StokesFlowProp').set('StokesFlowProp', 1, '0');
% CONSIDER BUOYANCY
model.physics('spf_FEM').feature('buoyancy').active(true);
% TRANSIENT 
model.physics('spf_FEM').prop('EquationForm').set('form', 1, 'Transient');
model.physics('ht_FEM').prop('EquationForm').set('form', 1, 'Transient');
% EVENTUALLY RUN NAVIER - STOKES PROBLEM
model.sol('NavStok').runAll;




% Computation Completed
fprintf('Navier-Stokes simulation completed in %f minutes \n', (toc-tii)/60 )

%% Results ------------------------------------------------------------------- %

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.5                                 date:   May  2013             %
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.5 - Big changes in the model architecture. A Stokes flow is   04/05/2013 %
%       solved both as reference field for POD and initia condition for the    %
%       subsequent FEM solution of Navier-Stokes. All these two solution are   %
%       stored into the model, each one into its respective study node.        %
%   0.4 - No longer base on exteral .mat files. All the useful data 09/04/2013 %
%       are included into the ph structure                                     %
%   0.3 - General rivision, added few comments                      27/03/2013 %
%   0.2 - Added varargin and a condition on 'run' command           20/03/2013 %
%   0.1 - kick-off                                                  27/02/2013 %
% ---------------------------------------------------------------------------- %
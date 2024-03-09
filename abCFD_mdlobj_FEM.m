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

function [ model, assemb ] = abCFD_mdlobj_FEM( ph,ModelVersion, varargin )

%% INITIALIZE ---------------------------------------------------------------- %
% Load input
load( 'geo.mat' );
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
ModelUtil.showProgress( true );
model = ModelUtil.create( ['Model',ModelVersion] );
model.modelPath( pwd ) ;
model.version( ModelVersion );
model.author( 'Montanari Mattia') ;

% Create model, geometry, mesh and physic
model.modelNode.create('mod1');
geom = model.geom.create('geom1', 2);
mesh = model.mesh.create('mesh1', 'geom1'); 

%% Single fluid flow phase --------------------------------------------------- %
singlephase = model.physics.create('spf_FEM', 'LaminarFlow', 'geom1');
singlephase.prop('CompressibilityProperty').set('Compressibility', 1, 'Incompressible');

    
% Create and define physics of the study node
study = model.study.create('std1');
transtudy = study.feature.create('time', 'Transient');
    
transtudy.activate('spf_FEM', true);
% Set up time ranges
transtudy.set('tlist', ...
    ['range(0,', num2str(ph.step_t),',',num2str(ph.t_final),')' ]);
  
% transtudy.set('tlist', ...
%     ['range(0,', num2str(ph.step_t),',',num2str(ph.t_final)...
%     ,') range(',num2str( ph.t_final+ph.step_t),',',...
%       num2str(ph.step_t/10),',', num2str(3*ph.t_final),')']);

%% Parameters' definitions --------------------------------------------------- %
abCFD_parameters( model, ph );

%% Functions' definitions ---------------------------------------------------- %
abCFD_functions( model );

%% Variables definition ------------------------------------------------------ %
abCFD_variables( model );

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

%% About the Physics --------------------------------------------------------- %
singlephase.prop('EquationForm').set('form', 1, 'Transient');
singlephase.prop('ShapeProperty').set('order_fluid', 1, num2str( V_order ));

%%  Heat Transfer in solids -------------------------------------------------- %
heatransfer = model.physics.create('ht_FEM', 'HeatTransfer', 'geom1');
heatransfer.feature.create('fluid1', 'FluidHeatTransferModel', 2);
heatransfer.feature('fluid1').selection.set([1]);
heatransfer.feature('fluid1').set('minput_velocity_src', 1, 'root.mod1.u');

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
    [ ph.fd.DIR.iwrite,'HT_DiricletRamp(t[1/s])'] );
% OUTLET   - SINGLE PHASE FLOW
singlephase.feature.create('out1', 'Outlet', 1);
singlephase.feature('out1').selection.set( ph.fd.Neumann );
% singlephase.feature('out1').set('BoundaryCondition', 1, ...
%     'PressureNoViscousStress');

% DIRICLET - HEAT TRANSFER
heatransfer.feature.create('Diricl', 'TemperatureBoundary', 1);
heatransfer.feature('Diricl').selection.set([ ph.ht.Diriclet ]);  
heatransfer.feature('Diricl').name('Temperature');
heatransfer.feature('Diricl').set(...
    'T0', 1, [ ph.ht.DIR.iwrite , 'HT_DiricletRamp(t[1/s])+initial_T'] );

% ROBIN    - HEAT TRANSFER
heatransfer.feature.create('Robin', 'ConvectiveCooling', 1);
heatransfer.feature('Robin').selection.set( ph.ht.Robin );  
heatransfer.feature('Robin').name('Convection');
heatransfer.feature('Robin').set('h', 1, 'HT_ROB_h');
heatransfer.feature('Robin').set('Text', 1, ...
    [ ph.ht.ROB.iwrite , 'HT_RobinRamp(t[1/s])'] );

% NEUMANN  - HEAT TRANSFER
heatransfer.feature.create('Neumann', 'HeatFluxBoundary', 1);
heatransfer.feature('Neumann').selection.set([ ph.ht.Neumann ]);  
heatransfer.feature('Neumann').name('HeatFlux');
heatransfer.feature('Neumann').set('q0', 1, ...
	[ ph.ht.NEU.iwrite , 'HT_NeumannRamp(t[1/s])'] );
% heatransfer.feature('Neumann').set('HeatFluxType', 1, 'GeneralInwardHeatFlux');


%% Mesh setting -------------------------------------------------------------- %
mesh.automatic(false);
model.mesh('mesh1').feature.remove('size1');
% model.mesh('mesh1').feature.remove('cr1');
model.mesh('mesh1').feature.remove('bl1');
model.mesh('mesh1').feature.remove('ftri1');
model.mesh('mesh1').feature.remove('ftri2');

if 1 == strcmp( el_type , 'Free_Quad' )
% *---*  Free quadrilateral mesh  
    mesh.feature.create('msh', 'FreeQuad');
    mesh.feature('size').set('hauto', num2str( el_grade ));
    
elseif 1 == strcmp( el_type , 'Free_Tria' )
% *---*  Free triangles
    mesh.feature.create('msh', 'FreeTri');
    
    % TRIVIAL MESH
    trivialmesh = 'false';
    if 1 == strcmp( trivialmesh , 'active') ;
        
        model.mesh('mesh1').feature('msh').feature.create('dis1', 'Distribution');
        model.mesh('mesh1').feature('msh').feature('dis1').selection.all;
        model.mesh('mesh1').feature('msh').feature('dis1').set('numelem', '1');
        model.mesh('mesh1').feature('size').set('custom', 'on');
        model.mesh('mesh1').feature('size').set('hmax', '5');
        model.mesh('mesh1').feature('size').set('hmin', '5');
        model.mesh('mesh1').run;
    else
        mesh.feature('size').set('hauto', num2str( ph.COMSOL.el_grade ));
    end
%     model.mesh('mesh1').feature('size').set('hgrad', '1');


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
    
%% SOLVER CONFIGURATION ------------------------------------------------------ %

% Create Transient solver
TransSol = 'TranSol';
    
    TransSolNode = model.sol.create( TransSol );
    TransSolNode.study('std1');
    TransSolNode.feature.create('st1', 'StudyStep');
    TransSolNode.feature('st1').set('study', 'std1');
    TransSolNode.feature('st1').set('studystep', 'time');

    % Set dependent variable(s)
    TransSolNode.feature.create('v1', 'Variables');
    TransSolNode.feature('v1').set('control', 'time');

    % Do not create default plots
    TransSolNode.detach;

% Create Stationary solver

StatSol = 'StatSol';
model.study('std1').feature.create('stat', 'Stationary');
StatSolNode = model.sol.create( StatSol );
StatSolNode.study('std1');
StatSolNode.feature.create('st1', 'StudyStep');
StatSolNode.feature('st1').set('study', 'std1');
StatSolNode.feature('st1').set('studystep', 'stat');
StatSolNode.feature.create('v1', 'Variables');
StatSolNode.feature('v1').set('control', 'stat');
StatSolNode.feature.create('s1', 'Stationary');
StatSolNode.feature('s1').feature.create('seDef', 'Segregated');
StatSolNode.feature('s1').feature.create('fc1', 'FullyCoupled');
StatSolNode.feature('s1').feature('fc1').set('initstep', 0.01);
StatSolNode.feature('s1').feature('fc1').set('minstep', 1.0E-6);
StatSolNode.feature('s1').feature('fc1').set('dtech', 'auto');
StatSolNode.feature('s1').feature('fc1').set('maxiter', 50);
StatSolNode.feature('s1').feature.create('d1', 'Direct');
StatSolNode.feature('s1').feature('d1').set('linsolver', 'pardiso');
StatSolNode.feature('s1').feature('fc1').set('linsolver', 'd1');
StatSolNode.feature('s1').feature('fc1').set('initstep', 0.01);
StatSolNode.feature('s1').feature('fc1').set('minstep', 1.0E-6);
StatSolNode.feature('s1').feature('fc1').set('dtech', 'auto');
StatSolNode.feature('s1').feature('fc1').set('maxiter', 50);
StatSolNode.feature('s1').feature.remove('fcDef');
StatSolNode.feature('s1').feature.remove('seDef');
StatSolNode.attach('std1');

%% TIME DEPENDENT SOLVER ----------------------------------------------------- %

TimeSolver = TransSolNode.feature.create('t1', 'Time');
TimeSolver.set('maxorder', '1'); % Time discret-ORDER
TimeSolver.set('control', 'time');
TimeSolver.set('plot', 'off');
TimeSolver.set('plotfreq', 'tout');
TimeSolver.set('probesel', 'all');
TimeSolver.set('probes', {});
TimeSolver.set('probefreq', 'tsteps');
TimeSolver.set('atolglobalmethod', 'scaled');
TimeSolver.set('atolglobal', 0.0010);
TimeSolver.set('estrat', 'exclude');
TimeSolver.set('maxorder', 2);
TimeSolver.set('control', 'time');

% Time stepping 
% Control tollerance. MIND: should be smaller than the result itself!
TimeSolver.set('eventtol', '0.001');
TimeSolver.set('maxorder', '1');            % Time-discretiazion order
TimeSolver.set('tout', 'tsteps');           % all timesteps taken by solver
TimeSolver.set('tstepsbdf', 'intermediate');% Take intermediate time steps

    se = 'NonSegregated';
    
if strcmp( se, 'Segregated');
    model.sol('Solver1').feature('t1').feature.create('se1', 'Segregated');
    model.sol('Solver1').runAll;
    disp('         ****  ---------------- ****')
    disp('         **** SEGREGATED SOLVER ****')
    disp('         ****  ---------------- ****')
else
    TimeSolver.feature.create('fc1', 'FullyCoupled');
    TimeSolver.feature('fc1').set('jtech', 'once');
    TimeSolver.feature('fc1').set('maxiter', 5);
    TimeSolver.feature.create('d1', 'Direct');
    TimeSolver.feature('d1').set('linsolver', 'pardiso');
    TimeSolver.feature('fc1').set('linsolver', 'd1');
    TimeSolver.feature('fc1').set('jtech', 'once');
    TimeSolver.feature('fc1').set('maxiter', 6 );
    TimeSolver.feature.remove('fcDef');
    TimeSolver.feature('aDef').set('convinfo', 'detailed');
end
% Avoid scaling
TransSolNode.feature('v1').set('scalemethod', 'none');
TimeSolver.set('atolglobalmethod', 'unscaled');
 
TransSolNode.attach('std1');
 

% Time-dependent solver    ****   ADVANCED *****  ---------------------------- %
% Control matrix singularity of mass matrix
TimeSolver.set('masssingular', 'yes');
% Constant consisten initialization of zero-filled diagonal entries in 
TimeSolver.set('consistent', 'off');

% AVOID PREOUDERING AND ROW EQUILIBRATION

model.sol('TranSol').feature('t1').feature('aDef').set('rowscale', 'off');
model.sol('TranSol').feature('t1').feature('d1').set('pardrreorder', 'off');
 
model.sol('StatSol').feature('s1').feature('d1').set('pardrreorder', 'off');
model.sol('StatSol').feature('s1').feature('aDef').set('rowscale', 'off');

%% ASSEMBLY NODES ------------------------------------------------------------ %
% Create assembly feature for transient solver
assembl = StatSolNode.feature.create('asmbl','Assemble');
for id_m = 1 : size( assemb.matrices , 2 )
    assembl.set(  assemb.matrices{id_m}    , 'on'  );
end
for id_v = 1 : size(assemb.vector , 2 )
    assembl.set(  assemb.vector{id_v}    , 'on'  );
end
% model.param.set('timestep','1');
% model.param.set('t', '0.1');

% Create assembly feature for transient solver
assembl = TransSolNode.feature.create('asmbl','Assemble');
for id_m = 1 : size( assemb.matrices , 2 )
    assembl.set(  assemb.matrices{id_m}    , 'on'  );
end
for id_v = 1 : size(assemb.vector , 2 )
    assembl.set(  assemb.vector{id_v}    , 'on'  );
end
model.param.set('timestep','1');
model.param.set('t', '0.1');
% ADVANCED
model.sol('TranSol').feature('asmbl').feature('aDef').set('rowscale', 'off');
model.sol('StatSol').feature('asmbl').feature('aDef').set('rowscale', 'off');


%  EVENTUALLY RUN
if numel(varargin) == 0
    fprintf('Modelobject setup complete after %f sec. \n Start Comsol computation at %s \n',...
    toc, datestr(now,15) )
    tii = toc;
    % NOT run this solver
    %     StatSolNode.attach('std1');
    % Deactivate this solver
    model.study('std1').feature('stat').active(false);
    % Run this solver
    TransSolNode.runAll;
else
    
    % FOR OPTIMIZE ROUTIN HERE I SHOULD HAVE ASSEMBLY FEATURE
    
    return
end

% Computation Completed
fprintf('Comsol simulation completed in %f minutes \n', (toc-tii)/60 )

%% Results ------------------------------------------------------------------- %
if 'Y' == 'teniamoliperdopo'
tiin = toc;
% PLOT TEMPERATURE SURFACE
% PlotGr = model.result.create('temperature', 'PlotGroup2D');
% PlotGr.name('Temperature');
% PlotGr.set('data', 'dset1');
% PlotGr.feature.create('surf1', 'Surface');
% PlotGr.feature('surf1').name('Surface');
% PlotGr.feature('surf1').set('colortable', 'ThermalLight');
% PlotGr.feature('surf1').set('data', 'parent');
% PlotGr.run;
% PLOT VELOCITY
PlotGr = model.result.create('Velocity', 2);
PlotGr.set('data', 'dset1');
PlotGr.feature.create('surf1', 'Surface');
PlotGr.feature('surf1').set('expr', {'spf.U'});
PlotGr.set('frametype', 'spatial');
PlotGr.name('Velocity (spf)');

PlotGr.feature.create('str1', 'Streamline');
PlotGr.feature('str1').set('expr', {'u' 'v'});
PlotGr.feature('str1').set('descr', 'Velocity field');
PlotGr.feature('str1').set('posmethod', 'magnitude');

PlotGr.run;
% PLOT PRESSURE
% PlotGr = model.result.create('Pressure', 2);
% PlotGr.set('data', 'dset1');
% PlotGr.feature.create('con1', 'Contour');
% PlotGr.feature('con1').set('expr', {'p'});
% PlotGr.set('frametype', 'spatial');
% PlotGr.name('Pressure (spf)');
% PlotGr.feature('con1').set('number', 40);
% PlotGr.run;


fprintf('Plots setup and in %f sec. \n', (toc-tiin) )



	% initialize values for time dependent variables:
    nt = ph.t_final/ph.step_t +1 ;
    range = ['range(2,', num2str(round( (nt-1)/7 ) + 1),',',num2str( nt-1 ),')'];
    
    % \Gamma_2 Normal heat flux in time (on Robin boundary)
    PlotGr_3 = model.result.create('pg3', 'PlotGroup1D');
    PlotGr_3.setIndex('looplevelinput', 'manualindices', 0);
    PlotGr_3.setIndex('looplevelindices', range , 0);
    PlotGr_3.set('titletype', 'none');
    PlotGr_3.name('Gamma2'); 
    PlotGr_3.set('xlabelactive', 'on');
    PlotGr_3.set('ylabelactive', 'on');
    PlotGr_3.set('xlabel', 'x-coordinate (m)');
    PlotGr_3.set('ylabel', 'Total normal heat flux (W/m<sup>2</sup>)');
    % plot flux at several instants
    LineInTime = PlotGr_3.feature.create('lngr1', 'LineGraph');%  plots in time
    LineInTime.selection.set([2]); %#ok<NBRAK>
    LineInTime.set('expr', 'ht.ntflux');
    LineInTime.set('xdataexpr', 'x');
    LineInTime.set('legend', 'on');
    LineInTime.set('legendmethod', 'automatic');
%     LineInTime.set('refine',  elements(2) );
    LineInTime.set('smooth', 'none');
    LineInTime.set('resolution', 'norefine');
    LineInTime.name('Normal flux frames');
    % plot flux at initial and final time
    LineOneEnd = PlotGr_3.feature.duplicate('lngr2', 'lngr1');
    LineOneEnd.set('data', 'dset1');
	LineOneEnd.setIndex('looplevelinput', 'manual', 0);
	LineOneEnd.setIndex('looplevel', ['1,',num2str(nt)] , 0);
    LineOneEnd.set('linewidth', '2');
	LineOneEnd.set('linemarker', 'cycle');    
    LineOneEnd.set('legendmethod', 'manual');
    LineOneEnd.setIndex('legends', 't = 0', 0);
    LineOneEnd.setIndex('legends', 't = LAST', 1); 
    LineOneEnd.name('Normal flux @t=0 and tfinal');
    LineOneEnd.set('markerpos', 'datapoints');
    PlotGr_3.run; 
    
    % \Gamma1 Temperature in time ( on Diriclet boundary )
    PlotGr_4 =  model.result.duplicate('pg4', 'pg3');
    PlotGr_4.name('Gamma1'); 
    PlotGr_4.feature('lngr1').selection.set([1]); 
    PlotGr_4.feature('lngr1').set('xdataexpr', 'y');
    PlotGr_4.feature('lngr1').set('expr', 'T');
    PlotGr_4.feature('lngr1').name('Temperature frames');
    PlotGr_4.feature('lngr1').set('resolution', 'norefine');
    PlotGr_4.feature('lngr1').name('Temperature frames');
    PlotGr_4.feature('lngr2').selection.set([1]); 
    PlotGr_4.feature('lngr2').set('xdataexpr', 'y');
    PlotGr_4.feature('lngr2').set('expr', 'T');
    PlotGr_4.feature('lngr2').set('markerpos', 'datapoints');
    PlotGr_4.feature('lngr2').set('resolution', 'norefine');
    PlotGr_4.feature('lngr2').name('Temperature @t=0 and tfinal');
    PlotGr_4.set('xlabelactive', 'on');
    PlotGr_4.set('xlabel', 'y-coordinate(m)');
    PlotGr_4.set('ylabelactive', 'on');
    PlotGr_4.set('ylabel', 'Temperature (K)');
    PlotGr_4.run; 
    
    % \Gamma_4 Normal heat flux in time (on Neumann boundary)
    PlotGr_5 =  model.result.duplicate('pg5', 'pg3');
    PlotGr_5.name('Gamma4'); 
    PlotGr_5.feature('lngr1').selection.set([4]); 
    PlotGr_5.feature('lngr1').set('xdataexpr', 'y');
    PlotGr_5.feature('lngr1').set('expr', 'ht.ntflux');
    PlotGr_5.feature('lngr1').name('Normal flux frames');
    PlotGr_5.feature('lngr2').selection.set([4]); 
    PlotGr_5.feature('lngr2').set('xdataexpr', 'y');
    PlotGr_5.feature('lngr2').set('expr', 'ht.ntflux');
    PlotGr_5.feature('lngr2').name('Normal flux @t=0 and tfinal');
    PlotGr_5.feature('lngr2').set('markerpos', 'datapoints');
    PlotGr_5.set('xlabel', 'y-coordinate(m)');
    
    % \Gamma3 Normal heat flux in time (on Insulated boundary)
    PlotGr_6 =  model.result.duplicate('pg6', 'pg3');
    PlotGr_6.name('Gamma3'); 
    PlotGr_6.feature('lngr1').selection.set([3]); 
    PlotGr_6.feature('lngr1').set('xdataexpr', 'x');
    PlotGr_6.feature('lngr1').set('expr', 'ht.ntflux');
    PlotGr_6.feature('lngr1').name('Normal flux frames');
    PlotGr_6.feature('lngr2').selection.set([3]); 
    PlotGr_6.feature('lngr2').set('xdataexpr', 'x');
    PlotGr_6.feature('lngr2').set('expr', 'ht.ntflux');
    PlotGr_6.feature('lngr2').name('Normal flux @t=0 and tfinal');
    
%% ADD PLOTS FOR EVALUATION INTERPOLATION FUNCTION

    % \Gamma_X Normal heat flux 
    PlotGr = model.result.create('eval_interp_hT', 'PlotGroup1D');
    PlotGr.setIndex('looplevelinput', 'first', 0);
    PlotGr.set('titletype', 'none');
    PlotGr.name('GradGammaX'); 
    PlotGr.set('xlabelactive', 'on');
    PlotGr.set('ylabelactive', 'on');
    PlotGr.set('xlabel', 'x-coordinate (m)');
    PlotGr.set('ylabel', 'Total normal heat flux (W/m<sup>2</sup>)');
    LineUnique = PlotGr.feature.create('lngr1', 'LineGraph'); 
    LineUnique.selection.set( 1 );
    LineUnique.set( 'legend', 'off');
    LineUnique.set('smooth', 'internal');
    LineUnique.set('resolution', 'fine');
    LineUnique.set('linewidth', '2');
    LineUnique.set('linemarker', 'none');
    LineInTime.name('Normal flux frames');
    LineUnique.active(false);
end

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.4                                 date:  April 2013             %
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.4 - No longer base on exteral .mat files. All the useful data 09/04/2013 %
%       are included into the ph structure                                     %
%   0.3 - General rivision, added few comments                      27/03/2013 %
%   0.2 - Added varargin and a condition on 'run' command           20/03/2013 %
%   0.1 - kick-off                                                  27/02/2013 %
% ---------------------------------------------------------------------------- %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \  	MODELOBJECT FOR 2D HEAT TRANSFER PROBLEM ON A  %
%  /----\ |  \|    |--  |   |   HORIZONTAL CHANNEL WITH A CYLINDER. ALL KIND   %
% /      \|__/ \__ |    |__/    OF BOUNDARY CONDITIONS PLUS HEAT SOURCE ARE    %
%                               APPLIED. 
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ model, matrices ] = abCFD_modelobject( ph , A_BuondaryCondCheck ) contains 
%   the whole set up of the model. All the  results and the plot group are
%   also defined. The 'model' variable is exported and ready to tranfer results.

function [ model, assemb ] = abCFD_modelobject( ph , A_BuondaryCondCheck, A_Save )
% Initizilize the extracted matrices
% matrices = ['uscale';'K     ';'L     ';'M     ';'N     ';'D     ';'E     ';...
%             'NF    ';'NP    ';'MP    ';'MLB   ';'MUB   ';'Kc    ';'Lc    ';...
%             'Dc    ';'Ec    ';'Null  ';'Nullf ';'ud    '];
assemb.matrices = {'K','N','D','E','NF','NP','Kc','Dc','Ec','Null','Nullf'};
assemb.vector = {'L','M','MP','MLB','MUB','Lc','ud','uscale'};

% Load input
load( 'input.mat' );

import com.comsol.model.*
import com.comsol.model.util.*
% Model generalities
ModelUtil.showProgress(true);
model = ModelUtil.create( ['Model',ModelVersion] );
model.modelPath('C:\Users\XMAMON\Desktop\20130204\heat');
model.version( ModelVersion );
model.author('Montanari Mattia');

% Create model, geometry, mesh and physic
model.modelNode.create('mod1');
geom = model.geom.create('geom1', 2);
mesh = model.mesh.create('mesh1', 'geom1');
heatransfer = model.physics.create('ht', 'HeatTransfer', 'geom1');
% Create and define the study node (holds nodes about: overall setting and
% solver sequence)
study = model.study.create('std1');
transtudy = study.feature.create('time', 'Transient');
transtudy.activate('ht', true);
transtudy.set('tlist', ...
    ['range(0,', num2str(ph.step_t),',', num2str(ph.t_final),')']);

%% Parameters' definitions
abCFD_parameters( model, ph, init_T );
% Diriclet boundary conditions
%% Functions' definitions
abCFD_functions( model );
%% Variables definition
abCFD_variables( model );

%% Geo file
geom.feature.create('r1', 'Rectangle');
geom.feature('r1').set('type', 'solid');
geom.feature('r1').set('base', 'corner');
geom.feature('r1').set('pos', {'0' '0'});
geom.feature('r1').set('lx', '0.4');
geom.feature('r1').set('ly', '2.2');
geom.feature.create('c1', 'Circle');
geom.feature('c1').set('type', 'solid');
geom.feature('c1').set('base', 'center');
geom.feature('c1').set('pos', {'0.2' '0.2'});
geom.feature('c1').set('r', '0.05');
geom.run;

%% Selections
model.selection.create('sel1');
model.selection('sel1').geom(2);
model.selection('sel1').name('Cylinder');
model.selection('sel1').set([2]); %#ok<NBRAK>

%% Material
model.material.create('mat1');
model.material('mat1').propertyGroup('def').set('thermalconductivity', {num2str( k )});
model.material('mat1').propertyGroup('def').set('density', {num2str( rho )});
model.material('mat1').propertyGroup('def').set('heatcapacity', {num2str( Cp )});

%% About the Physics
heatransfer.prop('EquationForm').set('form', 1, 'Transient');
heatransfer.prop('ShapeProperty').set('order_temperature', 1, num2str( el_order ));

%% Apply boundary contidions

    % Initial boundary conditions
heatransfer.feature('init1').set('T', 1, 'initial_T');
    % HEAT SOURCE
heatransfer.feature.create('hs1', 'HeatSource', 2);
heatransfer.feature('hs1').selection.named('sel1');
heatransfer.feature('hs1').selection.set([2]); %#ok<NBRAK> % Domain select
heatransfer.feature('hs1').set('Q', 1, num2str( ph.Q ));

    % DIRICLET ON \Gamma_1
heatransfer.feature.create('Diricl', 'TemperatureBoundary', 1);
heatransfer.feature('Diricl').selection.set([ ph.BC.Diriclet ]);  
heatransfer.feature('Diricl').name('DiricletBC');
heatransfer.feature('Diricl').set(...
    'T0', 1, [ph.g1.iwrite,'g_init(t[1/s])+initial_T'] );
%                 [g2.iwrite,'g2(t[1/s])*g_init(t[1/s])+initial_T'] );

    % ROBIN ON \Gamma_2 
heatransfer.feature.create('Robin', 'ConvectiveCooling', 1);
heatransfer.feature('Robin').selection.set([ ph.BC.Robin ]);  
heatransfer.feature('Robin').set('h', 1, 'BC_Robin_pi');
heatransfer.feature('Robin').set('Text', 1, [ph.g2.iwrite,'g_init(t[1/s])']);
heatransfer.feature('Robin').name('RobinBC');
heatransfer.feature('Robin').version( ModelVersion );

    % NEUMANN ON \Gamma_4
heatransfer.feature.create('Neumann', 'HeatFluxBoundary', 1);
heatransfer.feature('Neumann').selection.set([ ph.BC.Neumann ]);  
heatransfer.feature('Neumann').set('q0', 1, [ph.g4.iwrite,'g_init(t[1/s])']);
heatransfer.feature('Neumann').set('HeatFluxType', 1, 'GeneralInwardHeatFlux');

%% Mesh setting
mesh.automatic(false);
mesh.feature.remove('ftri1');
if 1 == strcmp( el_type , 'Free_Quad' )
% *---*  Free quadrilateral mesh  
    mesh.feature.create('msh', 'FreeQuad');
elseif 1 == strcmp( el_type , 'Free_Tria' )
% *---*  Free triangles
    mesh.feature.create('msh', 'FreeTri');
elseif 1 == strcmp( el_type , 'Map_Quad' );
% *---*  Structured grid. This requires the domains to be joint
    geom.feature.create('uni1', 'Union');
    geom.feature('uni1').selection('input').set({'r1' 'c1'});
    geom.feature('uni1').set('intbnd', 'off');
    geom.run;
    % Apply element distribution
    mesh.feature('msh').feature.create('dis1', 'Distribution');
    mesh.feature('msh').feature('dis1').set('numelem', '4');
    mesh.feature('msh').feature('dis1').selection.set([3 2]);
    mesh.feature('msh').feature.create('dis2', 'Distribution');
    mesh.feature('msh').feature('dis2').set('numelem', '10');
    mesh.feature('msh').feature('dis2').selection.set([1 4]);
    mesh.feature.create('msh', 'Map');
end
mesh.feature('size').set('hauto', num2str( el_grade ));
mesh.feature('msh').feature.create('dis1', 'Distribution');
mesh.feature('msh').feature('dis1').set('numelem', num2str( 11 - el_grade) );
mesh.feature('msh').feature('dis1').selection.set([ 3 ]); %#ok<NBRAK>
mesh.run;

%% General SOLVER SEQUENCE configuration
	% It requires the nodes: Study Step, depende variables, Solver
    
firstSolver = 'Solver1';
    
GlobalSolver_1 = model.sol.create( firstSolver );
GlobalSolver_1.study('std1');
GlobalSolver_1.feature.create('st1', 'StudyStep');
GlobalSolver_1.feature('st1').set('study', 'std1');
GlobalSolver_1.feature('st1').set('studystep', 'time');

% Set dependent variable(s)
GlobalSolver_1.feature.create('v1', 'Variables');
GlobalSolver_1.feature('v1').set('control', 'time');
% Avoid scaling
GlobalSolver_1.feature('v1').set('scalemethod', 'none');
GlobalSolver_1.feature('v1').feature('mod1_T').set('scalemethod', 'none');
% Set up Time-dependent solver
TimeSolver = GlobalSolver_1.feature.create('t1', 'Time');
TimeSolver.set('plot', 'off');
TimeSolver.set('plotfreq', 'tout');
TimeSolver.set('probesel', 'all');
TimeSolver.set('probes', {});
TimeSolver.set('probefreq', 'tsteps');
TimeSolver.set('atolglobalmethod', 'scaled');
TimeSolver.set('atolglobal', 0.0010);
TimeSolver.set('maxorder', 2);
TimeSolver.set('control', 'time');
TimeSolver.feature.create('fc1', 'FullyCoupled');
TimeSolver.feature('fc1').set('jtech', 'once');
TimeSolver.feature('fc1').set('maxiter', 5);
TimeSolver.feature.create('d1', 'Direct');
TimeSolver.feature('d1').set('linsolver', 'pardiso');
TimeSolver.feature('fc1').set('linsolver', 'd1');
TimeSolver.feature('fc1').set('jtech', 'once');
TimeSolver.feature('fc1').set('maxiter', 5);
TimeSolver.feature.remove('fcDef');
TimeSolver.feature('aDef').set('convinfo', 'detailed');

% Avoid scaling
TimeSolver.set('atolglobalmethod', 'unscaled');
TimeSolver.set('fieldselection', 'mod1_T');
TimeSolver.set('atolmethod', {'mod1_T' 'unscaled'});
TimeSolver.set('fieldselection', 'mod1_T');
TimeSolver.set('atoludotactive', {'mod1_T' 'on'});
TimeSolver.set('ewtrescale', 'off');

GlobalSolver_1.attach('std1');

% Create assembly feature
assembl = GlobalSolver_1.feature.create('asmbl','Assemble');
for id_m = 1 : size( assemb.matrices , 2 )
    assembl.set(  assemb.matrices{id_m}    , 'on'  );
end
for id_v = 1 : size(assemb.vector , 2 )
    assembl.set(  assemb.vector{id_v}    , 'on'  );
end

%  RUN
GlobalSolver_1.runAll;

disp('***************** COMSOL Transient done ******************')
disp(['Comsol Computation finisched at ', num2str(datestr(now,'HH:MM:SS')) ])
disp('**********************************************************')

%% Results
% PLOT TEMPERATURE SURFACE
PlotGr_1 = model.result.create('pg1', 'PlotGroup2D');
PlotGr_1.name('Temperature');
PlotGr_1.set('data', 'dset1');
PlotGr_1.feature.create('surf1', 'Surface');
PlotGr_1.feature('surf1').name('Surface');
PlotGr_1.feature('surf1').set('colortable', 'ThermalLight');
PlotGr_1.feature('surf1').set('data', 'parent');
PlotGr_1.run;
% PLOT ISOTHERMAL CONTOURS
PlotGr_2 = model.result.create('pg2', 'PlotGroup2D');
PlotGr_2.name('Isothermal contours');
PlotGr_2.set('data', 'dset1');
PlotGr_2.feature.create('con1', 'Contour');
PlotGr_2.feature('con1').name('Contour');
PlotGr_2.feature('con1').set('colortable', 'ThermalLight');
PlotGr_2.feature('con1').set('data', 'parent');
PlotGr_2.feature.create('arws1', 'ArrowSurface');
PlotGr_2.feature('arws1').name('Arrow surface');
PlotGr_2.feature('arws1').set('unit', {'' ''});
PlotGr_2.feature('arws1').set('color', 'gray');
PlotGr_2.feature('arws1').set('data', 'parent');
PlotGr_2.run;

switch A_BuondaryCondCheck
    case 'Y'

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
    
    
    % Eventually save model
   	if strcmp( A_Save, 'Y') == 1
        mphsave(model,['heat_mode14_',ModelVersion,'.mph'])
        disp(['NEW MODEL SAVED IN: ',cd,'heat_mode14_',ModelVersion,'.mph']);
    end
    
end

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

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.5                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.5 - Added a new plot for check the gradient of interolation   26/02/2013 %
%         functions.
%   0.4 - The definition of B.C.-type per geo's edge is made more   20/02/2013 %
%         easily by loading variables. However the rest of the code is not
%         adapted to this parametrization
%   0.3 - Added feature 'dis1' on edge 3 to avoid having mesh too   19/02/2013 %
%         coarse. If there's only one element on top edge other functions      %
%         will fail!                                                           %
%   0.2 - Better definition of the boundaries' plots. Each plot is  11/02/2013 %
%         showing the relative prescribed quantity. Moreover the variable
%         'iwrite' has been introduced to have a better control of the BC's
%         directly from the MAIN.m. An array 'matrices' is added as output
%   0.1 - kick-off                                                  04/02/2013 %
%  
% ---------------------------------------------------------------------------- %
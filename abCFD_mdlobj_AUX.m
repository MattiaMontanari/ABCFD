% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \  	MODELOBJECT FOR 2D AUXILIAR MODEL.             %
%  /----\ |  \|    |--  |   |   GEOMETRY AND MESH ARE IMPORTED FROM A PARENT   %
% /      \|__/ \__ |    |__/    MODEL. PHYSICS AND BOUNDARY CONDITIONS ARE     %
%                               USER-DEFINED                                   %
%  * * * CALLS * * *                                                           %
%           i. abCFD_parameters                                                %
%          ii. abCFD_functions                                                 %
%         iii. abCFD_variables                                                 %
%          iv. load( 'input.mat' )
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ model, assemb ] = abCFD_mdlobj_AUX( ph, {NewModelID} )
% INPUTS:
%   - ph           : common phase structure
%   - {NewModelID} : defines a new model ID to be used, different from the 
%                       one loaded in 'input.mat'                     

function [ model ] = abCFD_mdlobj_AUX( ph, model, varargin )

%% INITIALIZE ---------------------------------------------------------------- %

[el_type , T_order, V_order, P_order ] = abCFD_elements( ph.COMSOL.elements );
% Existing Mesh
FEMmesh = char(model.mesh.tags);
% check input
if numel( varargin ) ~= 0
   ModelVersion = varargin{ 1 };
end

import com.comsol.model.*
import com.comsol.model.util.*

model.modelNode.create('mod2');
model.geom.create('geom2', 2);

model.mesh.create('mesh2', 'geom2');

%% GEO AND MESH SETTING ------------------------------------------------------ %

model.mesh('mesh2').feature.create('imp1', 'Import');
model.mesh('mesh2').feature('imp1').set('source', 'sequence');
model.mesh('mesh2').feature('imp1').set('sequence', FEMmesh );
model.mesh('mesh2').feature('imp1').importData;


%% Single fluid flow phase --------------------------------------------------- %
singlephase = model.physics.create('spf_AUX', 'LaminarFlow', 'geom2');
singlephase.prop('CompressibilityProperty').set('Compressibility', 1, 'Incompressible');

%% Material ------------------------------------------------------------------ %
material = model.material.create('mat2');
material.propertyGroup('def').set('thermalconductivity', ph.mat.k );
material.propertyGroup('def').set('density', ph.mat.rho );
material.propertyGroup('def').set('heatcapacity', ph.mat.Cp );
material.propertyGroup('def').set('dynamicviscosity', ph.mat.mud );
material.propertyGroup('def').set('ratioofspecificheat', ph.mat.gamma );

%% About the Physics --------------------------------------------------------- %

% Discretization order
singlephase.prop('ShapeProperty').set('order_fluid', 1, num2str( V_order ));

%% Heat Transfer
heatransfer = model.physics.create('ht_AUX', 'HeatTransfer', 'geom2');
heatransfer.feature.create('fluid1', 'FluidHeatTransferModel', 2);
heatransfer.feature('fluid1').selection.set([1]);
heatransfer.feature('fluid1').set('minput_velocity_src', 1, 'root.mod2.u2');
heatransfer.prop('ShapeProperty').set('order_temperature', 1, num2str( T_order ));

%% Apply boundary contidions ------------------------------------------------- %

% Initial boundary conditions 
heatransfer.feature('init1').set('T2', 1, ph.ht.init );
singlephase.feature('init1').set('u2', [ ph.fd.init , ph.fd.init , 0 ]);
singlephase.feature('init1').set('p2', 1, '0');

% OUTLET   - SINGLE PHASE FLOW
singlephase.feature.create('out1', 'Outlet', 1);
singlephase.feature('out1').selection.set( ph.fd.Neumann );

% INLET    - SINGLE PHASE FLOW
singlephase.feature.create('inl1', 'Inlet', 1);
singlephase.feature('inl1').selection.set( ph.fd.Diriclet );
singlephase.feature('inl1').set('U0in', 1, ...
     ph.fd.DIR.iwrite );

% DIRICLET - HEAT TRANSFER
heatransfer.feature.create('Diricl', 'TemperatureBoundary', 1);
heatransfer.feature('Diricl').selection.set([ ph.ht.Diriclet ]);  
heatransfer.feature('Diricl').name('Temperature');
heatransfer.feature('Diricl').set(...
    'T0', 1, ph.ht.DIR.iwrite );

% ROBIN    - HEAT TRANSFER
heatransfer.feature.create('Robin', 'ConvectiveCooling', 1);
heatransfer.feature('Robin').selection.set( ph.ht.Robin );  
heatransfer.feature('Robin').name('Convection');
heatransfer.feature('Robin').set('h', 1, 'HT_ROB_h');
heatransfer.feature('Robin').set('Text', 1, ...
    ph.ht.ROB.iwrite  );

% NEUMANN  - HEAT TRANSFER
heatransfer.feature.create('Neumann', 'HeatFluxBoundary', 1);
heatransfer.feature('Neumann').selection.set([ ph.ht.Neumann ]);  
heatransfer.feature('Neumann').name('HeatFlux');
heatransfer.feature('Neumann').set('q0', 1, ...
	 ph.ht.NEU.iwrite  );


% %% OPTIONS
% % On the physics
% 
% % Equation form
% switch Solver
%     case 'TranSol'
% singlephase.prop('EquationForm').set('form', 1, 'Transient');
%     case 'StatSol'
% singlephase.prop('EquationForm').set('form', 1, 'Stationary');
% end
% 
% % On the solver
% % TranSol
% model.study('std1').feature('time').set('useadvanceddisable', 'on');
%  
% % StatSol
% StatSol = model.study('std1').feature('stat');
% % Set advanced physics
% StatSol.set('useadvanceddisable', 'on');
% % Disable physics
% StatSol.set('disabledphysics', {'ht_AUX' 'spf_AUX'});
% StatSol.activate('spf_AUX', false);
% StatSol.activate('ht_AUX', false);
% % Eneble physics
% model.study('std1').feature('time').set('disabledphysics', {});
% StatSol.activate('ht_AUX', true);







end
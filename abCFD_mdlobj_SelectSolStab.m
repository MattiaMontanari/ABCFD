% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	CHANGING THE MODEL OBJECT SELECTS WHEATER:     %
%  /----\ |  \|    |--  |   |    - USE THE STABILIZATION OR NOT, AND           %
% /      \|__/ \__ |    |__/     - USE STEADY OR TRANSIENT SOLVER              %
%                                                                              %
% * * REMIND * * IT DOES NOT RUN THE SOLVER, IT JUST CHANGES THE MODEL OBJ     %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%abCFD_ChangeSolPhy( model, SolTag, SPF, HT, Stability )
% INPUT:
%       - model     : Comsol's structure
%       - SolTag    : Sovler tag to activated: 
%                       if == 1 use Transient Solver   
%                       if == 2 use Stationary solver
%       - SPF       : SinglePhaseFlow physics' tag
%       - HT        : Heat Transfer physics' tag
%       - Stability : Specify if use or not stabilization techniques
%                       'Stable' to include default stabilization techniques
%                       'UnStable' to avoid stabilization techniques

function abCFD_mdlobj_SelectSolStab( model, SolTag, SPF, HT, Stability )

%% WALKIN COMSOL WITH NEW MODEL

% Initialize  
SolvTag = cell( model.sol.tags );
SolvTag = SolvTag{ SolTag };

%% MODIFY STABILIZATION TECHNIQUES

if     strcmp( Stability, 'UnStable' )
        % Avoid stabilization techniques
        model.physics( SPF ).prop('Stabilization').set('CrosswindDiffusion', 1, '0');
        model.physics( SPF ).prop('Stabilization').set('StreamlineDiffusion', 1, '0');
        model.physics( HT ).prop('HeatConsistentStabilization').set('heatStreamlineDiffusion', 1, '0');
elseif strcmp( Stability, 'Stable' )
        % Consider stabilization techniques
        model.physics( SPF ).prop('Stabilization').set('CrosswindDiffusion', 1, '1');
        model.physics( SPF ).prop('Stabilization').set('StreamlineDiffusion', 1, '1');
        model.physics( HT ).prop('HeatConsistentStabilization').set('heatStreamlineDiffusion', 1, '1');
else
        error('Uncorrect stabilization')
end

%% MODIFY SOLVER SETTING

switch SolTag
    case 1      % TRANSIENT SOLVER
        
        % Check if 'FakeTime' exists it yes deactivate it         
        AllVar = cell(model.variable.tags);
        if sum( cellfun(@(x) strcmp(x,'FakeTime'),AllVar) ) > 0
            model.variable('FakeTime').active(false);
        end
        
        % Activate transient
        model.study('std1').feature('time').active(true);
        % Deactivate strady
        model.study('std1').feature('stat').active(false);
        % Activate solver
        model.sol( SolvTag ).attach('std1');
     
    case 2      % STEADY SOLVER
        % Use steady solver
        % Check if 'FakeTime' already exists
        if isempty(find(ismember( cell(model.variable.tags), 'FakeTime')))
            model.variable.create('FakeTime');
        end
        model.variable('FakeTime').model('mod1');
        model.variable('FakeTime').set('t', 'timestep');
        % Deactivate transient
        model.study('std1').feature('time').active(false);
        % Activate strady
        model.study('std1').feature('stat').active(true);
        % Ensure Convergency: Number of iterations to converge
        model.sol( SolvTag ).feature('s1').feature('fc1').set('ntermauto', 'itertol');
        model.sol( SolvTag ).feature('s1').feature('fc1').set('niter', '1');
        % Avoid scaling
        model.sol( SolvTag ).feature('v1').set('scalemethod', 'none');
        model.sol( SolvTag ).feature('s1').feature('aDef').set('autorescale', 'off');
%         % For STOKES FLOW use linear solver
%         model.sol( SolvTag ).feature('s1').set('nonlin', 'off');
%         % Store linearization points
%         model.sol( SolvTag ).feature('s1').set('storelinpoint', 'on');
        % Detailed log file
        model.sol( SolvTag ).feature('s1').feature('aDef').set('convinfo', 'detailed');
        % Activate solver
        model.sol( SolvTag ).attach('std1'); 
end

end
% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.4                                 date:  APRIL 2013             % 
% ---                                                                      --- %
%   0.4 - Stability and Solver type are selected only within this   05/04/2013 %
%       rountine by the user. 
%   0.3 - Scaling, linearziation points for steady solver           03/04/2013 %
%   0.2 - Convergence criteria for steady solver                    02/04/2013 %
%   0.1 - kick-off                                                  21/03/2013 %
% ---------------------------------------------------------------------------- %
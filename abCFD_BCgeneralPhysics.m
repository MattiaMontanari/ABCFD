% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \  	ASSEMBLES STRUCTURES TO BE USED IN THE         %
%  /----\ |  \|    |--  |   |   MODELOBJECT. IT DEALS WITH BOUNDARY CONDITIONS %
% /      \|__/ \__ |    |__/    OF A GENERAL PHYSICS.                          %
%                               ALSO BODY FORCES ARE DEFINED                   %
%   * * * CALLS * * *                                                          %
%           i. abCFD_SlopeRamp                                                 % 
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ phase ] = abCFD_BCgeneralPhysics( phase, dt )
% INPUTS
%   - phase : usual structure with phase's info
%   - dt    : time step taken 
% OUTPUT
%   - phase : usual structure with phase's info


function [ phase ] = abCFD_BCgeneralPhysics( phase, dt )
% Initialize

gDIR = phase.DIR.iwrite;
gNEU = phase.NEU.iwrite;
gROB = phase.NEU.iwrite;
gF = phase.Q;

%%  DEFINE BOUNDARY CONDITIONS ----------------------------------------------- %

% DIRICLET BC's
    %  abCFD:   T(x,t) = g(x,t)               COMSOL:  T = T0(x,t)
%Define g(x,t)
phase.DIR.iwrite = gDIR ;
%Define starting function
phase.DIR.start.cutoff = 1;
phase.DIR.start.t_initial = 0;
phase.DIR.start.t_final = dt;
phase.DIR.start.slope = abCFD_SlopeRamp( phase.DIR.start );

% ROBIN BC's
    %  abCFD:            pi * dT = h*(g - T)      COMSOL:  k*dT = h*(Text - T)
%Define g(x,t), k and h
phase.ROB.iwrite = gROB ;
phase.ROB.pi     = 1 ;
phase.ROB.kappa  = 1;
%Define starting function
phase.ROB.start.cutoff = 1;
phase.ROB.start.t_initial = 0;
phase.ROB.start.t_final = dt;
phase.ROB.start.slope = abCFD_SlopeRamp( phase.ROB.start );
% NEUMANN BC's 
    %  abCFD:            pi * dT = g(x,t)      COMSOL:  k*dT  = dT0
%Define g(x,t) and k
phase.NEU.iwrite = gNEU ;
phase.NEU.pi  = 1;
%Define starting function
phase.NEU.start.cutoff = 1;
phase.NEU.start.t_initial = 0;
phase.NEU.start.t_final = dt;
phase.NEU.start.slope = abCFD_SlopeRamp( phase.NEU.start );

% INSULATION (NEUMANN)
    %  abCFD:   dT(x,t) = 0                    COMSOL:  k*dT  = 0
    %  for:     dT0 = 0                     these formulation are equal 
  
%% DEFINE SPECIFIC FORCING TERM ---------------------------------------------- %
    %  abCFD:   Q                              COMSOL:  Q
phase.Q = gF ;        % Constant in specific unit: [ [Q] /m^3 ]

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.3                                 date:   MAY  2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.3 - Inputs modified                                           02/05/2013 %
%   0.2 - General revision, some comment added                      27/03/2013 %
%   0.1 - kick-off                                                  27/02/2013 %
% ---------------------------------------------------------------------------- %
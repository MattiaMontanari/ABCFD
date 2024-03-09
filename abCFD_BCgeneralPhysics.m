% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \  	ASSEMBLES STRUCTURES TO BE USED IN THE         %
%  /----\ |  \|    |--  |   |   MODELOBJECT. IT DEALS WITH BOUNDARY CONDITIONS %
% /      \|__/ \__ |    |__/    OF A GENERAL PHYSICS.                          %
%                               ALSO BODY FORCES ARE DEFINED                   %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ physics ] = abCFD_BCgeneralPhysics( physics, gDIR, gNEU, gROB, gF, dt )
% INPUTS
%   See parent function: abCFD_BCtoComsol
% OUTPUT
%   - physics: structure with fields detailing Diriclet, Robin and Neumann 
%               boundary conditions for any general physics.


function [ physics ] = abCFD_BCgeneralPhysics( physics, gDIR, gNEU, gROB, gF, dt )

%%  DEFINE BOUNDARY CONDITIONS ----------------------------------------------- %

% DIRICLET BC's
    %  abCFD:   T(x,t) = g(x,t)               COMSOL:  T = T0(x,t)
%Define g(x,t)
physics.DIR.iwrite = gDIR ;
%Define starting function
physics.DIR.start.cutoff = 1;
physics.DIR.start.t_initial = 0;
physics.DIR.start.t_final = dt;
physics.DIR.start.slope = abCFD_SlopeRamp( physics.DIR.start );

% ROBIN BC's
    %  abCFD:            pi * dT = h*(g - T)      COMSOL:  k*dT = h*(Text - T)
%Define g(x,t), k and h
physics.ROB.iwrite = gROB ;
physics.ROB.pi     = 1 ;
physics.ROB.kappa  = 1;
%Define starting function
physics.ROB.start.cutoff = 1;
physics.ROB.start.t_initial = 0;
physics.ROB.start.t_final = dt;
physics.ROB.start.slope = abCFD_SlopeRamp( physics.ROB.start );
% NEUMANN BC's 
    %  abCFD:            pi * dT = g(x,t)      COMSOL:  k*dT  = dT0
%Define g(x,t) and k
physics.NEU.iwrite = gNEU ;
physics.NEU.pi  = 1;
%Define starting function
physics.NEU.start.cutoff = 1;
physics.NEU.start.t_initial = 0;
physics.NEU.start.t_final = dt;
physics.NEU.start.slope = abCFD_SlopeRamp( physics.NEU.start );

% INSULATION (NEUMANN)
    %  abCFD:   dT(x,t) = 0                    COMSOL:  k*dT  = 0
    %  for:     dT0 = 0                     these formulation are equal 
  
%% DEFINE SPECIFIC FORCING TERM ---------------------------------------------- %
    %  abCFD:   Q                              COMSOL:  Q
physics.Q = gF ;        % Constant in specific unit: [ [Q] /m^3 ]

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.2                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.2 - General revision, some comment added                      27/03/2013 %
%   0.1 - kick-off                                                  27/02/2013 %
% ---------------------------------------------------------------------------- %
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
function [ ph, el_type , el_order ] = abCFD_BCtoComsol( ph, elements )
% The Structure ph includes the following fields:
%       ph.step_t = time step 
%       ph.t_final = phase's end time
%   BOUNDARY CONDITIONS ------------------------------------------------------ %
% DIRICLET BC's on \Gamma_1
    %  abCFD:   T(x,t) = g1(x,t)               COMSOL:  T = T0
    %  for:     T0 = g2(x,t)             	these formulation are equal 
ph.g1.cutoff = 1;          ph.g1.t_initial = 0;       ph.g1.t_final = ph.step_t;
ph.g1.iwrite = '1/4*sin(3*2*pi*s+2*pi*t[1/s])*sin(4*pi*t[1/s])*';  % 5*s*(1-s)*sin(3/2*pi*t[1/s]*t[1/s])
% ROBIN BC's on \Gamma_2 
    %  abCFD:            pi * dT = g2 - T      COMSOL:  k*dT = h*(Text - T)
    %  for:     h = 1
    %           Text(x,t) = g2(x,t)
    %           k = pi                    	these formulation are equal 
ph.g2.cutoff = 1;          ph.g2.t_initial = 0;       ph.g2.t_final = ph.step_t;
ph.g2.pi = 1;              ph.g2.iwrite = '1*g2(t[1/s])*';
% INSULATION (NEUMANN) BC's on \Gamma_3 
    %  abCFD:   dT(x,t) = 0                    COMSOL:  k*dT  = 0
    %  for:     dT0 = 0                     these formulation are equal 
% NEUMANN BC's on \Gamma_4
    %  abCFD:            pi * dT = g4           COMSOL:  k*dT  = dT0
    %  for:     k = pi
    %           dT = g4(x,t)                these formulation are equal 
ph.g4.cutoff = 1;          ph.g4.t_initial = 0;       ph.g4.t_final = ph.step_t; 
ph.g4.pi = 1;              ph.g4.iwrite = '1*g4(t[1/s])*';   

%   HEAT SOURCE  ------------------------------------------------------------- %
% Heat source in the cylinder. Constant in time and space.
    %  abCFD:   ???                          COMSOL:  ??
ph.Q = 0;        % Constant heat source                             [W/m^3]


%%  PRE-COMPUTATIONS
% Mesh definition
[el_type , el_order ] = abCFD_elements( elements );
% Function definition
ph.g1.slope = abCFD_SlopeRamp(ph.g1.cutoff , ph.g1.t_final , ph.g1.t_initial );
ph.g2.slope = abCFD_SlopeRamp(ph.g2.cutoff , ph.g2.t_final , ph.g2.t_initial );
ph.g4.slope = abCFD_SlopeRamp(ph.g4.cutoff , ph.g4.t_final , ph.g4.t_initial );
%   PHASE 1 - SET UP --------------------------------------------------------- %
% Include the four boundary conditions

% ph = struct('g1',g1 ,'g2', g2 ,'g4', g4 );

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.2                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.2 - Better definition of the boundaries' plots. Each plot is  11/02/2013 %     
%         showing the relative prescribed quantity. Moreover the variable
%         'iwrite' has been introduced to have a better control of the BC's
%         directly from the MAIN.m. An array 'matrices' is added as output
%   0.1 - kick-off                                                  04/02/2013 %     
%  
% ---------------------------------------------------------------------------- %
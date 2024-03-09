% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \  	COMPUTES AND PRINTS DIMENSIOLESS PARAMETERS    %
%  /----\ |  \|    |--  |   |   FOR FORCED CONVECTION AROUND A CYLINDER        %
% /      \|__/ \__ |    |__/                                                   %
%                                                                              %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ nuv] = abCFD_pi( L,  ph )
% INPUTS:
%   - L     : Characteristic lenght scale
%   - ph    : phase structure 

function [nuv] = abCFD_pi( L, ph)

%% INITIALIZE

rho = ph.mat.rho;
mud = ph.mat.mud;
k = ph.mat.k;
Cp = ph.mat.Cp;
thex = ph.mat.thex;
Tsur = ph.Tsur;
Tinit = ph.Tinit;
Umax = ph.Umax;
g = ph.g;

nuv = mud / rho;    % Kinematic viscosity   [m^2/s]
U_mean = 2/3 * Umax;% U mean at inlet due to the velocity profile

%% FROM LITERATURE
% REcritic = 2 * 10^5;
% REdrag = 10;  For RE < 10 the friction drag dominate
% REpres = 5000;  For RE > 5000  the pressure drag dominate

%% COMPUTE DIMENSIONLESS PARAMETERS
% Reynolds number: Internal forces Vs viscous forces
RE = round( U_mean * L / nuv );
% Grashof number: Bouyancy Vs viscous forces
GR = round( ( thex * L^3 * g * (Tsur - Tinit) ) / nuv^2 );
% Prandtl: Hydraulic and thermal property of the fluid
PR = mud * Cp / k;
% P?clet: advection effect oveer diffusion of same quantity
PE = RE * PR;
% Rayleigh number: indicates the way the heat is transfered
RA = round(GR*PR);
% Nusselt averaged around the cylinder
NU = 0.3 + ((0.62*RE^(1/2)*PR^(1/3))/(1+(0.4/PR)^(2/3) ))*(1+(RE/282.000)^(5/8))^(4/5);
%               ref. Churchill and Bernstein: Average Heat Transfer  ?30%
%               Coefficient actually this changes along the cylinder! 

%% PRINT TO SCREEN THE COEFF.
fprintf('\n><  DIMENSIONLESS PARAMETERS \n\t<> Reynolds:\t %g \t||\t%d\n', RE, RE )
fprintf('\t<> Grashof:\t\t %d \t||\t%e \n\t<> Prandtl:\t\t %f \t||\t%d\n',  GR , GR, PR, PR)
fprintf('\t<> Peclet:\t\t %f \t||\t%d \n', PE,PE )
fprintf('\t<> Nusselt:\t\t %f \t||\t%d\n', NU,NU)
fprintf('\t<> Rayleigh:\t  %d \t||\t%d\n', RA,RA)

%% OTHER DIMENSIONELSS NUMBERS NOT USEFUL IN FORCED CONVECTION
% Froude: Only if there's free surface/static pressure
% FR = U_mean^2 / ( G0 * L );
% Strouhal number: 
% ST = L / ( U_mean * t0 );



end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.3                                 date:  APRIL 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.4 - Introduced in the inputs the sturcture 'ph' which has     05/04/2013 %
%       all the material properties required.                                  %
%   0.3 - General revision. Few comments added                      27/03/2013 %
%   0.2 - General imprements                                        08/03/2013 %
%   0.1 - kick off. the Nusselt number is obtained from the web:    04/03/2013 %
% http://wwwme.nchu.edu.tw/Enter/html/lab/lab516/Heat%20Transfer/chapter_7.pdf %
% ---------------------------------------------------------------------------- %
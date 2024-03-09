% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \  	USER INTERFACE TO SET UP BOUNDARY CONDITIONS   %
%  /----\ |  \|    |--  |   |   FOR ALL THE PHYSICS. EACH PHYSICS ALLOWS ANY   %
% /      \|__/ \__ |    |__/    BOUNDARY CONDITIONS TYPE: DIRICLET, NEUMANN,   %
%                               AND ROBING. SOURCE TERMS AND BODY FORCES ARE   %
% TO BE DEFINED HERE AS A SIMPLE STRING.                                       %
% * * * CALLS * * *                                                            %    
%        i. abCFD_BCgeneralPhysics                                             %    
%       ii. abCFD_BCgeneralPhysics                                             %    
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ ph ] = abCFD_BCtoComsol( ph )
% INPUTS:
%   - ph    : Common phase structure
% OUTPUT:
%   - ph    : Common phase structure

function [ ph ] = abCFD_BCtoComsol( ph )

%%  HEAT TRANSFER PROBLEM ---------------------------------------------------- %

% Boundary conditions and heat source
htDIR = [ num2str( ph.Tmax-ph.Tbul ) , '*' ];
htNEU = '100*';
htROB = '0*';
htSRC = '0';
[ ph.ht ] = abCFD_BCgeneralPhysics( ph.ht, htDIR, htNEU, htROB, htSRC, ph.step_t );
% Initial conditions
ph.ht.init = ph.Tbul;         % homogeneous initial temperature               [K]

%%  MASS TRANSFER ------------------------------------------------------------ %

% Boundary conditions and heat source
fdDIR = [ num2str( ph.Umax ) , '*4*s*(1-s)*' ];
fdNEU = '1*';
fdROB = '1*';
fdSRC = [ num2str(-ph.g*ph.mat.thex*ph.mat.rho),'*0*( 0 -',num2str( ph.Tbul ),')'];
[ ph.fd ] = abCFD_BCgeneralPhysics( ph.fd, fdDIR, fdNEU, fdROB, fdSRC, ph.step_t );
% Initial conditions
ph.fd.init = 0;         % homogeneous initial Velcity


%%      CLOSURE
fprintf('Boundary conditions'' setting completed in %f sec. \n',toc)

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
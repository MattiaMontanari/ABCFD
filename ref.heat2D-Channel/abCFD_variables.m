% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SUB-ROUTINE TO SET COMSOL's VARIABLES FOR A    %
%  /----\ |  \|    |--  |   |   NAME AND VALUES OF EACH PARAMTER               %
% /      \|__/ \__ |    |__/    AD DEFINED AD-HOC                              %
%                                                                              %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%abCFD_variables( model ) adds to the current model a set of variables

function abCFD_variables( model )

%% VARIABLES FOR BOUNDARY FLUX EVALUATIONS

% On \Gamma1
gamma = '1';
IntVar = model.variable.create( ['var',gamma] );
IntVar.model( 'mod1' );
IntVar.set( [ 'ImprtTgamma',gamma] , ['InterpGamma',gamma,'(x,y)[1/K]'] );

IntVar.active(false);

% On \Gamma2
gamma = '2';
IntVar = model.variable.create( ['var',gamma] );
IntVar.model( 'mod1' );
IntVar.set( [ 'ImprtTgamma',gamma] , ['InterpGamma',gamma,'(x,y)[1/K]'] );

IntVar.active(false);

% On \Gamma3
gamma = '3';
IntVar = model.variable.create( ['var',gamma] );
IntVar.model( 'mod1' );
IntVar.set( [ 'ImprtTgamma',gamma] , ['InterpGamma',gamma,'(x,y)[1/K]'] );

IntVar.active(false);

% On \Gamma4
gamma = '4';
IntVar = model.variable.create( ['var',gamma] );
IntVar.model( 'mod1' );
IntVar.set( [ 'ImprtTgamma',gamma] , ['InterpGamma',gamma,'(x,y)[1/K]'] );

IntVar.active(false);

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 1.0                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick off                                                  25/02/2013 %     
% ---------------------------------------------------------------------------- %
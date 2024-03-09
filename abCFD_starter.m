% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \  	STARTS THE TIME COUNTER FOR THE WHOLE          %
%  /----\ |  \|    |--  |   |   SIMULATION. IT ALSO SHOWS THE LOGO AND CLEANS  %
% /      \|__/ \__ |    |__/    UP MONITOR AND WORKSPACE                       %
%                                                                              %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%abCFD_starter Starts() 
%       - Not inputs nor outputs.

function abCFD_starter 
 
% close all; 
colordef black; 
clc; 
addpath('Imported');    
set(0,'DefaultFigureWindowStyle','docked') 

clc;
disp(' >      >     >    >   >  > > >> x << < <  <   <    <     <      <')
disp('   /\   |     __  ___  __      > x <       ')
disp('  /  \  |__  /   |    |  \     > x <   COMPUTE HEAT MASS STRANSFER   ')
disp(' /----\ |  \|    |--  |   |    > x <   PROBLEM WITH COMSOL AND ROM    '  )
disp('/      \|__/ \__ |    |__/     > x <   ')
disp(' >      >     >    >   >  > > >> x << < <  <   <    <     <      <')

fprintf('\n Computation starts at %s  \n',  datestr(now,15) )

%% Time monitor
tic

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.3                                 date:  AUGUST 2013            % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.3 - Directory name change to 'Imported'                       16/08/2013 %
%   0.2 - Added general initialization lines and general resion     27/03/2013 %
%   0.1 - kick off                                                  01/03/2013 %
% ---------------------------------------------------------------------------- %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \    SUB-ROUTINE TO CALCULATE THE SLOPE OF A LINEAR %
%  /----\ |  \|    |--  |   |   FUNCTION PROVIDED THE CUT-OFF VALUE, THE       %
% /      \|__/ \__ |    |__/    INITIAL AND FINAL VALUES ON THE X-AXIS         %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[slope] = abCFD_SlopeRamp( inputstructure ) returns the slope of linear 
% function which passes by the two points. The points are defined with the 
% following coordinates:
%   1st point: [ 't_initial'    ,     0     ]
%   2nd point: [ 't_final'      , 'cutOff'  ]
% 
%   The slope is calculated as: slope = cutOff / ( t_final - t_initial);

function [slope] = abCFD_SlopeRamp( inputstructure )

cutOff = inputstructure.cutoff;
t_initial = inputstructure.t_initial;
t_final = inputstructure.t_final;
 
%% CALCULATE SLOPE
slope = cutOff / ( t_final - t_initial);

if slope == Inf
    slope = 0;
end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 1.0                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   1.0 - kick off.                                                 27/02/2013 %
% ---------------------------------------------------------------------------- %
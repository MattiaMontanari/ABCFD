% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \    SUB-ROUTINE TO CALCULATE THE SLOPE OF A LINEAR %
%  /----\ |  \|    |--  |   |   FUNCTION PROVIDED THE CUT-OFF VALUE, THE       %
% /      \|__/ \__ |    |__/    INITIAL AND FINAL VALUES ON THE X-AXIS         %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%abCFD_SlopeRamp returns the slope of linear function which passes by the
%   two points. The points are defined with the following coordinates:
%   1st point: [ 't_initial'    ,     0     ]
%   2nd point: [ 't_final'      , 'cutOff'  ]
% 
%   The slope is calculated as: slope = cutOff / ( t_final - t_initial);
% 
%   If t_initial is not specified, by default is equal to zero

function [slope] = abCFD_SlopeRamp( cutOff, t_final, varargin )


%% INPUT CHECK
% If t_initial is not speficied this is assumed to be zero
nVarargs = nargin;
t_initial = 0;
if nVarargs > 2;
    t_initial = varargin{1};
end
 
%% CALCULATE SLOPE
slope = cutOff / ( t_final - t_initial);

if slope == Inf
    slope = 0;
end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 1.0                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   1.0 - IMPLEMENTATION FOR A SIMPLE COMSOL RAMP FUNCTION.         01/02/2013 %
%         NO SMOOTH ZONE IS CONSIDERED                                         %
% ---------------------------------------------------------------------------- %
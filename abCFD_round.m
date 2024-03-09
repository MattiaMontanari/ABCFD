% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	ROUNDS ALL THE ENTRIES OF A MATRIX UP TO A     %
%  /----\ |  \|    |--  |   |   DEFINED PRECISION                              %
% /      \|__/ \__ |    |__/                                                   %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%A = abCFD_round( A, {tol})
% INPUT:
%       - tol is a small value, eg 1e-7, which defines the rounding precision
%       - A is the matrix to be rounded

function A = abCFD_round( A, varargin )

if numel( varargin ) == 0
    tol = 1e-8;
else
    tol = varargin{ 1 };
end
A = A./tol;
A = round(A);
A = A.*tol;
end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.2                                 date:  April 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.2 - Made 'tol' optional input                                 09/04/2013 %     
%   0.1 - kick-off.                                                 21/03/2013 %     
% ---------------------------------------------------------------------------- %
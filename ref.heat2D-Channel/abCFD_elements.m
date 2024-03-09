% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \   	SUB-ROUTINE WHICH TRANSLATES FOR COMSOL THE    %
%  /----\ |  \|    |--  |   |   USER'S REQUIREMENT OF USING A PARTICULAR KIND  %
% /      \|__/ \__ |    |__/    OF GRID USING A SPECIFIC DEGREE FOR THE SHAPE  %
%                               FUNCTIONS OF THE ELEMENTS. THE GRID HAS ONLY   %
% ONE TYPE OF ELEMENT.                                                         %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[el_type , el_order ] = abCFD_elements( elements ) possible entries are
%   strings with 'T'(triangular) or 'Q'(quadratic) or 'M'(mapped) as first
%   entry and the second entry is an integer within (1,5) to specify the
%   shape functions' order.

function [el_type , el_order ] = abCFD_elements( elements )

switch elements(1)
    case 'T'
        el_order = str2double( elements(2) );
        el_type = 'Free_Tria';
    case 'Q'
        el_order = str2double( elements(2) );
        el_type = 'Free_Quad';        
    case 'M'
        el_order = str2double( elements(2) );
        el_type = 'Map_Quad';    
        
end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick-off                                                  07/02/2013 %     
% ---------------------------------------------------------------------------- %
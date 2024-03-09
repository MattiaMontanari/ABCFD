% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \   	SUB-ROUTINE WHICH TRANSLATES FOR COMSOL THE    %
%  /----\ |  \|    |--  |   |   USER'S REQUIREMENT OF USING A PARTICULAR KIND  %
% /      \|__/ \__ |    |__/    OF GRID USING A SPECIFIC DEGREE FOR THE SHAPE  %
%                               FUNCTIONS OF THE ELEMENTS. THE GRID HAS ONLY   %
% ONE TYPE OF ELEMENT.                                                         %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[el_type , T_order,V_order, P_order ] = abCFD_elements( elements )
% INPUTS: 
%   - elements: string which defines mesh elements and shape functions
%               IF ==  'T'(triangular) 
%               IF ==  'Q'(quadratic) 
%               IF ==  'M'(mapped)
%             as first entry. the second entry is an integer within (1,5) 
%             to specify the shape functions' order.

function [el_type , T_order,V_order, P_order ] = abCFD_elements( elements )

    switch elements(1)
        case 'T'
            el_type = 'Free_Tria';
        case 'Q'
            el_type = 'Free_Quad';        
        case 'M'
            el_type = 'Map_Quad';    
    end

    T_order = str2double( elements(2) );
    V_order = str2double( elements(4) );

    if V_order == 1
        P_order = V_order;
    else
        P_order = V_order - 1;
    end

end
% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.3                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.3 - General revision, few comments added                      27/03/2013 %     
%   0.2 - P_order added as output                                   28/02/2013 %     
%   0.1 - kick-off                                                  07/02/2013 %     
% ---------------------------------------------------------------------------- %
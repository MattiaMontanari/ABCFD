% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \   	GIVEN A PARTICULAR ELEMENT THIS ROUTINE GIVES  %
%  /----\ |  \|    |--  |   |   THE TAGS AND THE GLOBAL COORDINATES OF THE     %
% /      \|__/ \__ |    |__/    VERTECES OF THE ELEMENT                        %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%vtx = abCFD_elem_FindVtx ( DofTagGlobal, DofCordGloba, ElType, FieldOrd)
% INPUT:
%   - DofTagGlobal: Ordered dofs' of an element
%   - DofCordGloba: Global coordinate of each ordered dof
%   - ElType:       Element type, e.g. 'tri','quad'
%   - FieldOrd:     Element interp. order
 
function vtx = abCFD_elem_FindVtx ( DofTagGlobal, DofCordGloba, ElType,FieldOrd)

switch ElType 
    case 'tri'
        % Store Global coordinates
        vtx.GlobalCoord( : , 1 ) = DofCordGloba( : , 1);
        vtx.GlobalCoord( : , 2 ) = DofCordGloba( : , 1 + FieldOrd );
        vtx.GlobalCoord( : , 3 ) = DofCordGloba( : , end );
        % Store Global tag
        vtx.tag( 1 ) = DofTagGlobal( 1 );
        vtx.tag( 2 ) = DofTagGlobal( 1 + FieldOrd );
        vtx.tag( 3 ) = DofTagGlobal( end);
    case 'quad'
        vtx.GlobalCoord( : , 1 ) = DofCordGloba( : , 1);
        vtx.GlobalCoord( : , 2 ) = DofCordGloba( : , 1   + FieldOrd );
        vtx.GlobalCoord( : , 3 ) = DofCordGloba( : , end - FieldOrd );
        vtx.GlobalCoord( : , 4 ) = DofCordGloba( : , end );
        % Store Global tag
        vtx.tag( 1 ) = DofTagGlobal( 1 );
        vtx.tag( 2 ) = DofTagGlobal( 1   + FieldOrd );
        vtx.tag( 3 ) = DofTagGlobal( end - FieldOrd );
        vtx.tag( 4 ) = DofTagGlobal( end);
end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   0.1 - kick-off                                                  22/03/2013 %     
%  
% ---------------------------------------------------------------------------- %
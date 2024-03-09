% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SUB-ROUTINE TO ASSEMBLE A SPARCE MATRIX IN     %
%  /----\ |  \|    |--  |   |         %
% /      \|__/ \__ |    |__/      %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ RobinContribute ] = abCFD_AssembleEdgCntrb 

function [ EleCntrb ] = abCFD_AssembleEleCntrb( dom ,dof, ele, elements )

dom = 2; 
if elements(1) == 'M'
    disp('Exit. You are dealing with structured grid which has NOT heat source')
    return
end

% Obtain dof's included in domain 'dom'
DofDom = ele.DOFToNode( : , ele.GeoToNode == dom )';

%% Initialize
% Obtain number of elements and connectivity for elements on \Gamma
[NumEle , NumLocalDof] = size( DofDom ) ;

%% Loop over each element
% Initialized output
EleCntrb = sparse( dof.number_V1,  dof.number_V1);
for id_ele = 1 : NumEle
    % Map row, i.e. test functions
    map_i = DofDom( id_ele, : );
    % Map columns, i.e. shape functions
    map_j = DofDom( id_ele, : );
    % Coordinates of nodes of id_ele-th element on \Gamma
    XYcoord =  dof.dofCoord( :, DofDom( id_ele, :)  );
    % Define integral boundaries
    xa = min( XYcoord(1, : ));
    xb = max( XYcoord(1, : ));
    ya = min( XYcoord(2, : ));
    yb = max( XYcoord(2, : ));
    % Compute local coeff. matrix
    [ Q ] = abCFD_GPquad2D( NumLocalDof , xa, xb,ya,yb, ele,XYcoord, 'Y');
    % Map local coeff. matrix
    EleCntrb( map_i , map_j ) = EleCntrb( map_i , map_j ) + Q;
end

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick-off. Routine NOT WORKING AND NOT COMPLETED.          19/02/2013 %     
%         Mistakes in the formulationa and boubts about why should I do
%         this myself. the numerical integration might be hard for a
%         general purpose of different elements.
%  
% ---------------------------------------------------------------------------- %
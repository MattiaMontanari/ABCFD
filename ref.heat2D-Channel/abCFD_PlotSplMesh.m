% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \    SUB-ROUTINE TO PLOT THE SIMPLE MESH.           %
%  /----\ |  \|    |--  |   |   THE SUPPORTED ELEMENTS ARE: TRIANGLES and      %
% /      \|__/ \__ |    |__/    QUADRILATERAL ELEMENTS (FREE or MAPPED)        %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%abCFD_PlotSplMesh( Answer , coord , connectivity, elements, varargin ) 
%   plots in new window a 2D SIMPLE MESH with N vertexes defining M polygons
%   of type S = 'T' or'M' or 'Q'
%       where: coord        [ N x 2 ]   array
%              elements     [ 1 x 2 ]   string
%                   S = elements(1)
%              connectivity [ M x S ]
% abCFD_PlotSplMesh( Answer , coord , connectivity, elements, varargin ) 
%   plots in the FCth windows
 
% EXAMPLE:
%   abCFD_PlotSplMesh( 'Y' , mesh_data.vertex' , Tag.ELM', elements, 3 )

function abCFD_PlotSplMesh( Answer , coord , connectivity, elements, varargin )
%% PLOT SIMPLE MESH

if Answer == 'Y'
    % Select which windows should be used
    nVarargs = nargin - 4;
    if nVarargs == 0
        figure();
        FC_tag = gcf;
    elseif nVarargs == 1
        FC_tag = varargin{1};
    end
    
%   Reorder connectivity in case of QUADRILATERAL elements

    switch elements(1)
    % The Simple mesh is affected only by the TYPE of element
        case 'T'    
            % Nothing to do
        case 'Q'    % Consider Free Quadrilateral element
            commutation = [1,2,4,3];
            connectivity = connectivity( : , commutation );
        otherwise  % Consider MAPPED Quadrilateral element
            commutation = [1,2,4,3];
            connectivity = connectivity( : , commutation );
    end
    
    % Plotting    
    figure( FC_tag )
    title('Simple Mesh')    
    axis equal
    patch('faces', connectivity  ,'Vertices',coord);
    
    
end



end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick-off                                                  07/02/2013 %     
%  
% ---------------------------------------------------------------------------- %
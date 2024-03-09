% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \    SUB-ROUTINE TO PLOT THE EXTENDED MESH GIVING   %
%  /----\ |  \|    |--  |   |   A 3D TIME-VARYING SURFACE OF THE SOLUTION      %
% /      \|__/ \__ |    |__/    VECTOR. WORKS FINE ONLY FOR LINEAR AND T2      %                     
%                               ELEMENTS                                       %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%abCFD_PlotExtMesh( model, A , coord , connectivity, Z , elements  )
%   plots in new window the 3D EXTENDED MESH in n time-frames with N vertexes 
%   defining M polygons of type S = 'T' or'M' or 'Q'
%       where: coord        [ N x 2 ]   array
%              elements     [ 1 x 2 ]   string
%                   S = elements(1)
%              connectivity [ M x S ]
%              Z            [ N x t ]
%abCFD_PlotExtMesh( model, A , coord , connectivity, Z , elements, FC )
%   plots in the FCth windows
 
% EXAMPLE:
% abCFD_PlotExtMesh( model, 'Y', dof.Coord_V1', elem.DOFToNode', U , elements)

function abCFD_PlotExtMesh( model, A , coord , connectivity, Z , elements, varargin )
%% PRE-COMPUTATIONS
BoundingBox = model.geom('geom1').getBoundingBox;

%% PLOT EXTENDED MESH

if A == 'Y'
    % Select which windows should be used
    nVarargs = nargin - 6;
    if nVarargs == 0
        figure();
        FC_tag = gcf;
    elseif nVarargs > 0
        FC_tag = varargin{1};
    end
    
%   Reorder connectivity 

    switch elements
        case 'T1'
            % Nothing to do
        case 'T2'
            commutation = [ 1,2,3,5,6,4 ];
            connectivity = connectivity( : , commutation );
        case 'T3'
            commutation = [ 1,2,3,4,7,9,10,8,5,6 ];
            connectivity = connectivity( : , commutation );
        case 'Q1'
            commutation = [1,2,4,3];
            connectivity = connectivity( : , commutation );
        case 'Q2'
            commutation = [1,2,3,6,9,8,7,4,5];
            connectivity = connectivity( : , commutation );
        case 'Q3'
            commutation = [1,2,3,4,8,12,16,15,14,13,9,10,11,7,6,5];
            connectivity = connectivity( : , commutation );    
        case 'M1'
            commutation = [1,2,4,3];
            connectivity = connectivity( : , commutation );
        case 'M2'
            commutation = [1,2,3,6,9,8,7,4,5];
            connectivity = connectivity( : , commutation );
    end

    % Setup plot  
    figure( FC_tag )
    hold on
    title('Extended mesh')
    view( 49 , 29 )
    axis equal
    axis([ BoundingBox(1) BoundingBox(2) BoundingBox(3) BoundingBox(4) ...
                min(min(Z)) max(max(Z))])
    axis tight
    xlabel('X-coord. [m]')
    ylabel('Y-coord. [m]')
%     set(gcf,'position',get(0,'screensize'))

	% Plot in time
%     patch('faces', connectivity  ,'Vertices',coord);
    
    ExtendedMesh = patch('faces', connectivity  ,'Vertices',[coord , Z(:,1)],...
        'FaceColor','interp','LineWidth',1,'EdgeColor','interp',...
        'FaceVertexCData', Z(:,1) );
    for id_t = 1 : size( Z , 2 )
        set( ExtendedMesh , 'Vertices', [ coord ,Z(:,id_t)] );
        pause( 0.1 )
    end
    
end



end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.2                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.2 - COMMENT MARK ON: 'patch' (TO AVOID THE PLOT OF THE        08/02/2013 %
%   	SIMPLE MESH), AND 'set' TO AVOID THE FULL SCREEN VIEW.                 %
%   0.1 - kick-off                                                  07/02/2013 %     
%  
% ---------------------------------------------------------------------------- %
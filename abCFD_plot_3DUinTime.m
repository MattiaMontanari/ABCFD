% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \  	3D PLOT OF A SOLUTION FIELD, WHICH IS  A TIME  %
%  /----\ |  \|    |--  |   |   SCALAR OR VECTORIAL FIELD.                     %
% /      \|__/ \__ |    |__/    IN CASE OF VECTORIAL FILED THE COMPONENTS ARE  %
%                               COMBINED AS FOLLOS:  sqrt( var1.^2 + var2.^2 ) %
%                                                                              %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%abCFD_plot_3DUinTime(model,mesh,var,i_mesh,i_ele,i_dof,data,TimeSteps,{FC_tag})
% INPUTS:
%   - model     : common structure
%   - mesh      : common structure
%   - var       : common structure
%   - i_mesh    : identifies the mesh tag into the structyure 'mesh'
%   - i_ele     : identifies the element tag into the structyure 'mesh'
%   - i_dof     : identifies the DOF(s) tag into the structyure 'mesh'
%                 IF == c1 (scalar) plot the scalar field with tag c1
%                 IF == [c1,c2] (vector) plot the vector field with tag c1 & c2
%   - data      : can be a structure with field data.d1 or an array with
%       solution in time. 
%                 IF == U OR U.d1 plot scalar solution U in time
%                 IF == cat(3,U.d1,U.d1) plots vector solution U in time
%   - TimeSteps : time values taken by the solver
%   - {FC_tag}  : eventually specify figure's tag

function abCFD_plot_3DUinTime( model, mesh, var, i_mesh, i_ele, i_dof, data, TimeSteps, varargin)

%% INITIALIZE
% X-Y bounding for plot
BoundingBox = model.geom('geom1').getBoundingBox;
% Select which windows should be used
if isempty( varargin )
    figure();
    FC_tag = gcf; close ( FC_tag ) 
else
    FC_tag = varargin{1};
end

% Check data input
if isfield( data,'d1') ; data = data.d1; end

% Initialize new plot
figure( FC_tag )
hold on
view( 2)
axis equal
% Coordinates
X = var.dofCoord( 1 , :)';
Y = var.dofCoord( 2 , :)';
    
%% PLOT SCALAR FIELD
if ndims ( data ) == 2                  %#ok<*ISMAT>
    
    % extract data from 'mesh'
    Tags1   = mesh{ i_mesh }.ele{ i_ele }.dof{ i_dof }.uniqTAG;
    Name1   = mesh{ i_mesh }.ele{ i_ele }.dof{ i_dof }.name;
    Type1   = mesh{ i_mesh }.ele{ i_ele }.dof{ i_dof }.types;
    NumDof1 = mesh{ i_mesh }.ele{ i_ele }.dof{ i_dof }.numb;
    
    % Data to be plotted, based 'i_var'
    Z = data( Tags1, : );
    coord = [ X( Tags1 , 1) , Y( Tags1 , 1)];
    
    % Triangulation
    tri = delaunay( coord );
    
    % Axes and labels of the plot
    axis equal
    axis( [BoundingBox', min(min( Z )) , max(max( Z )) ])
        % axis tight
    axis tight
    xlabel('X-coord. [m]')
    ylabel('Y-coord. [m]')
    
    % Setup view
    view(3); view( -120 , 16)
    %   P L O T   I N   T I M E
    if ~isscalar( TimeSteps )
        % Print info
        fprintf( 'SCALAR PLOT INFO: \n Mesh type %s with %d DOFS \n', Type1,NumDof1)
        % Plot at initial time
        ExtendedMesh = patch('faces', tri  ,'Vertices',[coord , Z(:,1)],...
            'FaceColor','interp','LineWidth',1,'EdgeColor','interp',...
            'FaceVertexCData', Z(:,1),'CData', Z(:,1));
        
        % Loop plot on time steps
        for id_t = 1 : numel( TimeSteps ) 
            
            set( ExtendedMesh , 'Vertices', [ coord ,Z(:,id_t)],'CData', Z(:,id_t) );
            pause( 0.001 ) % REQUIRED!
            title( ['Variable ',Name1,' Solution at time = ' , num2str( TimeSteps(id_t))])
            
        end
        
    %   S I N G L E    P L O T    
    elseif isscalar( TimeSteps )
        % Plot at single time
        patch('faces', tri  ,'Vertices',[coord , Z(:,TimeSteps)],...
            'FaceColor','interp','LineWidth',1,'EdgeColor','interp',...
            'FaceVertexCData', Z(:,TimeSteps),'CData', Z(:,TimeSteps));
        title( ['Solution''s tag = ' , num2str( TimeSteps )])
    end
    
%% PLOT VECTOR FIELD
elseif ndims ( data ) == 3
    
    % Variables tags
    if numel( i_dof ) == 2; dof1 = i_dof(1); dof2 = i_dof(2); end
    % extract data from 'mesh'
    Tags1   = mesh{ i_mesh }.ele{ i_ele }.dof{ dof1 }.uniqTAG;
    Tags2   = mesh{ i_mesh }.ele{ i_ele }.dof{ dof2 }.uniqTAG;
    Name1   = mesh{ i_mesh }.ele{ i_ele }.dof{ dof1 }.name;
    Name2   = mesh{ i_mesh }.ele{ i_ele }.dof{ dof2 }.name;
    Type1   = mesh{ i_mesh }.ele{ i_ele }.dof{ dof1 }.types;
    Type2   = mesh{ i_mesh }.ele{ i_ele }.dof{ dof2 }.types;
    NumDof1 = mesh{ i_mesh }.ele{ i_ele }.dof{ dof1 }.numb;
    NumDof2 = mesh{ i_mesh }.ele{ i_ele }.dof{ dof2 }.numb;
    
    % Data to be plotted, based 'i_var'
    Z1 = data( Tags1, : , 1);
    Z2 = data( Tags2, : , 2);
    coord = [ X( Tags1 , 1) , Y( Tags1 , 1)];
    
    Z = sqrt( Z1.^2 + Z2.^2);
    
    % Triangulation
    tri = delaunay( coord );
    
    % Axes and labels of the plot
%     axis equal
    axis( [BoundingBox', min(min( Z )) , max(max( Z )) ])
        
    axis tight
    xlabel('X-coord. [m]')
    ylabel('Y-coord. [m]')
    
    % Setup view
    view(3); view( -120 , 16)
    %   P L O T   I N   T I M E
    if ~isscalar( TimeSteps )
        % Print info
        fprintf( 'VECTORIAL PLOT INFO: \n Mesh type %s with %d DOFS \n', Type1,NumDof1)
        % Plot at initial time
        ExtendedMesh = patch('faces', tri  ,'Vertices',[coord , Z(:,1)],...
            'FaceColor','interp','LineWidth',1,'EdgeColor','interp',...
            'FaceVertexCData', Z(:,1),'CData', Z(:,1));
        
        % Loop plot on time steps
        for id_t = 1 : numel( TimeSteps ) 
            
            set( ExtendedMesh , 'Vertices', [ coord ,Z(:,id_t)],'CData', Z(:,id_t) );
            pause( 0.001 )% REQUIRED!
title( ['Variables ',Name1,' && ',Name2,' Solution at time = ' , num2str( TimeSteps(id_t))])
            
        end
        
    %   S I N G L E    P L O T    
    elseif isscalar( TimeSteps )
        % Print to screen info
        fprintf( 'VECTORIAL PLOT INFO: \n Mesh type %s with %d DOFS \n', Type1,NumDof1)
                % Plot at single time
        patch('faces', tri  ,'Vertices',[coord , Z(:,TimeSteps)],...
            'FaceColor','interp','LineWidth',1,'EdgeColor','interp',...
            'FaceVertexCData', Z(:,TimeSteps),'CData', Z(:,TimeSteps));
        title( ['Variables ',Name1,' && ',Name2])
        
    end
end


end
% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  April 2013             %
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick-off                                                  17/04/2013 %
% ---------------------------------------------------------------------------- %
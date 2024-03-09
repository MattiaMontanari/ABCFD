% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	ROUTINE TO PLOT THE MESH FROM SmeshInfo DATA   %
%  /----\ |  \|    |--  |   |   AND VARIOUS SCALAR FIELDS ON TOP OF THE MESH   %
% /      \|__/ \__ |    |__/                                                   %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[FC]=abCFD_plot_MeshWithTags( 
% Texts on figure 'FC' the 'TAGS' with coordinates 'coord' using 'color' to
% plot the tags on the quadrant 'quad' = [1:4]

function [ FC ] = abCFD_plot_MeshWithTags( SmeshInfo, Fields, mesh, var )

%% PLOT MESH
% open new window
figure(); FC = gcf;
% Create data structure
toto = struct( 't',SmeshInfo.data.elem{2},'p',SmeshInfo.data.vertex,...
                        'd1',zeros(1, size( SmeshInfo.data.vertex,2) ));
% Plot mesh
pat = abCFD_patch( toto, 1);
set(pat,'EdgeColor', 'black','FaceColor',[.5 .5 .5]);

%% PLOT SOLUTION FIELDS

% Define color map
coloring = [ 1 0 0; 0 1 0; 0 0 1; 1 1 0; rand(  numel( Fields),3) ];
% Initialize legend
legendary = blanks(6);
legendary(1,1:5) = [ SmeshInfo.stat.meshtag];
% Loop over the scalar field specied in 'Fields'
for i_f = 1 : numel( Fields)
    % Get dofs' tags
    TAGS = unique( mesh{2}.ele{1}.dof{ Fields(i_f) }.TAG, 'stable');
    % Get dofs' coords
    coord = var.dofCoord( :,TAGS );
    % Plot tags on the mesh
    abCFD_plot_SingleTagOnMesh( coord, FC, coloring(i_f, :), i_f, TAGS );
    % Useless plot, just to give proper color to the legend
    plot( coord(1,1), coord(2,1),'color', coloring(i_f, :));
    legendary = [legendary; mesh{2}.ele{1}.dof{  i_f  }.name{1} ]; %#ok<AGROW>
end
% Show legend
legend(legendary);
end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick-off                                                  06/02/2013 %     
%  
% ---------------------------------------------------------------------------- %
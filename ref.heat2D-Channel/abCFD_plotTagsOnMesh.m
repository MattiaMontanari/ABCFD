% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SUB-ROUTINE TO PLOT ON THE MESH PLOT ANY KIND  %
%  /----\ |  \|    |--  |   |   OF ORDERED TAG PROVIDED AN ORDERED COUPLE OF   %
% /      \|__/ \__ |    |__/    VECTORS COINTAINIG THE TAGs' COORDINATES       %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ FC_mesh ] = abCFD_plotTagsOnMesh( model,X , Y , 'title' ) plots into a new 
%   window the mesh and and the i-th tag of the point with 
%   coordinates [X(i) , Y(i)]. X and Y are two row vectors of equal size. 
%[ FC_mesh ] = abCFD_plotTagsOnMesh( model,X , Y , FC_mesh, 'color') does the 
%   plot in the window figure numbered 'FC_mesh' plotting the tags using 'color'

% Example: abCFD_plotTagsOnMesh( model,mesh_data.vertex(1,:)' ,mesh_data.vertex(2,:)' , 'mesh' )
function [ FC_mesh ] = abCFD_plotTagsOnMesh( model, X , Y , varargin )

%% INPUT CHECK
% Put X and Y in a row vector if not correctly inseret
if size(X,1) > size(X,2)
    X = X';
    Y = Y';
end

% If FC_mesh is not speficied this a new figure will be created and
% color as well as allignment parameter are defined by default.
nVarargs = nargin;
if nVarargs == 5;
    FC_mesh = varargin{1};
    color = varargin{2};
    VerticalAlignment = 'top';
    HorizontalAlignment = 'left';
else
    figure()
    FC_mesh = gcf;
    title(varargin{1});
    color = 'w';
    VerticalAlignment = 'bottom';
    HorizontalAlignment = 'right';
end

figure( FC_mesh );
mphmesh(model, 'mesh1','Edgemode','off','edgecolor', 'b','facemode','on' );   
hold on ;      axis off

for id_n = 1 : size( X , 2 )
        text( X(id_n) , Y(id_n)     ,   num2str(id_n)     	, ...
            'VerticalAlignment'     ,   VerticalAlignment   ,...
            'HorizontalAlignment'   ,   HorizontalAlignment ,...
            'color'     , color   )
end


% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick-off                                                  06/02/2013 %     
%  
% ---------------------------------------------------------------------------- %
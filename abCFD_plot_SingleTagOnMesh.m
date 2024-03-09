% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SUB-ROUTINE TO PLOT ON THE MESH PLOT ANY KIND  %
%  /----\ |  \|    |--  |   |   OF ORDERED TAG PROVIDED COORDINATED AND TAGS   %
% /      \|__/ \__ |    |__/                                                   %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[FC]=abCFD_plot_SingleTagOnMesh(model,coord,FC,color,quad,TAGS)
% Texts on figure 'FC' the 'TAGS' with coordinates 'coord' using 'color' to
% plot the tags on the quadrant 'quad' = [1:4]

function [ FC ] = abCFD_plot_SingleTagOnMesh( coord, FC, color, quad, TAGS )

%% INPUT CHECK
% Turn 'coord' in a row vector if not correctly inseret
if size(coord,1) > size(coord,2)
    coord = coord';
end
X = coord( 1 , :);
Y = coord( 2 , :);
% Check if is specified the figure handle number
if FC == 0
    figure()
    FC = gcf;
    title('Mesh');
end
  
% Define the text position based on 'quad' - quadrant

switch quad
    case 1
    	VerticalAlignment = 'top';
        HorizontalAlignment = 'left';
    case 2
        VerticalAlignment = 'top';
        HorizontalAlignment = 'right';
    case 3
        VerticalAlignment = 'bottom';
        HorizontalAlignment = 'left';
    case 4
        VerticalAlignment = 'bottom';
        HorizontalAlignment = 'right';     
end

figure( FC );
% mphmesh(model, 'mesh1','Edgemode','off','edgecolor', 'b','facemode','on' );   
hold on ;      axis off; axis tight;

for id_n = 1 : size( X , 2 )
        text( X(id_n) , Y(id_n)     ,   num2str( TAGS( id_n,1 ) ) , ...
            'VerticalAlignment'     ,   VerticalAlignment   ,...
            'HorizontalAlignment'   ,   HorizontalAlignment ,...
            'color'     , color   );
end


% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.2                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.2 - Changed to sub-routine in order to support another        09/03/2013 %
%         routine to plot multiple tags. Inputs have been changed
%   0.1 - kick-off                                                  06/02/2013 %     
%  
% ---------------------------------------------------------------------------- %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SUB-ROUTINE TO PLOT ISOCONTOUR SURFACES        %
%  /----\ |  \|    |--  |   |           %
% /      \|__/ \__ |    |__/                                          %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%FC = abCFD_isocontour( sol, NumSol)

function abCFD_isocontour( data, time , varargin )

%% Initialize
% Interpolation resulution - point per direction
Npoints = 50;
Nisolines = 25;
% Problem size
[ GeoDim,  NumNode] = size( data.p );
% x-coordinate
x = data.p(1,:)';
% Y-coordinate
y = data.p(2,:)';
% Define z-coordinate
switch nargin 
    case 2  % Within the strct 'data'
        % z-coordinate
        z = data.d1(time ,:)' ;
        % Unit
        unit = data.unit;
    case 3
        z = varargin{1};
        z = z( : , time );
        unit = '-';
    case 4
        z = varargin{1};
        z = z( : , time );
        unit = '-';
        
        if isnumeric( varargin{2} )
                figure( varargin{2} )
        end         
end


% figure;

% Perform reguar mapping using delunay triang
DT = DelaunayTri(x,y);
    % F1 = TriScatteredInterp(x,y,z); 
% Interpolate z-coordinates on Delaunary Triangulation
F2 = TriScatteredInterp(DT, z);
% New map and meshgrid for plot
xlin = linspace( min(x), max(x), Npoints);
ylin = linspace( min(y), max(y), Npoints);
[X,Y] = meshgrid(xlin,ylin);
% Calculate interupolated values
    % Z1 = F1(X,Y);
Z2 = F2(X,Y);
% Plot
surf(X,Y,Z2,'EdgeColor', 'none' )
    % mesh(X,Y,Z)
hold on
contour(X,Y,Z2,Nisolines, 'LineColor', 'black','LineWidth',1.5)
plot3( x,y,z,'.','MarkerSize',0.5,'Color','blue' )
axis equal
axis tight


end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 1.0                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   1.0 - kick-off                                                  28/02/2013 %
% ---------------------------------------------------------------------------- %
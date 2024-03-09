% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SUB-ROUTINE TO PLOT ISOCONTOUR SURFACES        %
%  /----\ |  \|    |--  |   |           %
% /      \|__/ \__ |    |__/                                          %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%FC = abCFD_streamline( sol, NumSol)

function abCFD_streamline( x,y,u,v )

% sol.d1 = sqrt( U.u.d1( time,:).^2' + U.v.d1( time,:).^2' );

% x = sol.p(1,:)';
% y = sol.p(2,:)';
% u = u(time, :)';
% v = v(time, :)';
% z = sol.d1(time,:)';

DT = DelaunayTri(x,y);

% F1 = TriScatteredInterp(x,y,z); 
Fu = TriScatteredInterp(DT, u);
Fv = TriScatteredInterp(DT, v);


xlin = linspace( min(x), max(x), 50);
ylin = linspace( min(y), max(y), 50);

[X,Y] = meshgrid(xlin,ylin);

U = Fu(X,Y);
V = Fv(X,Y);
  
% hold on
% plot3( x,y,sqrt( u.^2 + v.^2 ),'.','MarkerSize', 1 , 'Color','black');
% hold on

h = streamline( X, Y, U,V, linspace(min(x),max(x),20), linspace(min(y),min(y),20)  ) ;
% 'LineColor', 'black','LineWidth',1.5)
set(h,'Color','white','LineWidth', 1)
% axis equal
% axis tight

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 1.0                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   1.0 - kick-off                                                  28/02/2013 %
% ---------------------------------------------------------------------------- %
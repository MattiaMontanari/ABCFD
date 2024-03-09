function [phi,phi_x,phi_y] = eval_N_TRI_2D(XY,xi)
%Evaluates the shape functions and their derivatives on an element for some
%points specified in the parent element (usually hammer points)
%
%
%xi is a 2xN matrix of points coordinates in the parent element
%
%Author: Adrien Leygue (adrien.leygue@ec-nantes.fr)
%Last modification: 04/03/2013

%Check the inputs 
assert( all(size(XY)==[3 2]),'coordinates of the element as 3x2 matrix');

NP = size(xi,2);
phi = [1-xi(1,:)-xi(2,:); xi(1,:); xi(2,:)];

    dphidxi = [-1 1 0; -1 0 1]';
    dxdxi = XY'*dphidxi;
    dphidx = dphidxi/dxdxi;
    phi_x = dphidx(:,1);
    phi_y = dphidx(:,2);
end
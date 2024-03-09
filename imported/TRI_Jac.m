function J = TRI_Jac(XY)
%Evaluates the jacobian of a triangular element at some
%points in the parent element (usually hammer points)
%
%xi is a 2xN matrix of points coordinates in the parent element
%
%Author: Adrien Leygue (adrien.leygue@ec-nantes.fr)
%Last modification: 16/10/2012

%Check the inputs 
assert(all(size(XY)==[3 2]),'coordinates of the element as 2x3 matrix');

dphidxi = [-1 1 0; -1 0 1]';
dxdksi = XY'*dphidxi;
J = (det(dxdksi)); 
end


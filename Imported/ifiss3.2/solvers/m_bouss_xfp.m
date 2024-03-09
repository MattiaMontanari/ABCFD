function y = m_bouss_xfp(x_it,aparams,mparams);
%M_BOUSS_XFP modified ideal PCD preconditioner
%   y = m_bouss_xfp(x_it,aparams,mparams);
%   input
%          x_it         operand for preconditioning operator
%          aparams      structure defining coefficient matrix
%          mparams      structure defining preconditioning matrix
%   output
%          y            result of preconditioning operation
%
%   IFISS function: DJS; 27 April 2012.
% Copyright (c) 2012 D.J. Silvester, M.L. Mihajlovic.

nv = length(aparams.F);
np = size(aparams.B,1); 
nt = size(aparams.Mt,2);
Gdiag=spdiags(diag(mparams.G),0,nv,nv);
rv=x_it(1:nv); rp=x_it(nv+1:nv+np); rt=x_it(nv+np+1:nv+np+nt);

%% temperature solve
zt = aparams.Ft \ rt;

%% pressure Poisson setup
xB = (Gdiag\aparams.B')';
BBt = aparams.B*xB'; 

%% pressure
zp = -(BBt)\((mparams.Fp)*((mparams.Mp)\rp));

%% velocity solve
rv = rv +(aparams.Mt)*zt - (aparams.B')*zp;
zv = aparams.F \ rv;
y = [zv;zp;zt];
return

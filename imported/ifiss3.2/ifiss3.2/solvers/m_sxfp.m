function y = m_sxfp(x_it,aparams,mparams)
%M_SXFP modified stabilized ideal preconditioner
%   y = m_sxfp(x_it,aparams,mparams);
%   input
%          x_it         operand for preconditioning operator
%          aparams      structure defining coefficient matrix
%          mparams      structure defining preconditioning matrix
%   output
%          y            result of preconditioning operation
%
%   IFISS function: HCE; 29 December 2009.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage

nv = length(aparams.F);
%nu = nv/2;
np = size(aparams.B,1);


rv=x_it(1:nv); rp=x_it(nv+1:nv+np);

%% pressure Poisson setup
Gdiag=spdiags(diag(mparams.G),0,nv,nv);
xB = (Gdiag\aparams.B')';
BBt = aparams.B*xB' + (.5/mparams.viscosity)*mparams.Cp1;

%% pressure
if mparams.domain==1,
   n_null = mparams.n_null;
   minor = [1:n_null-1,n_null+1:np]';
   yp = (mparams.Fp)*((mparams.Mp)\rp);
   zp = zeros(np,1);
   zp(minor) = - BBt(minor,minor) \ yp(minor);
else
   zp = -(BBt)\((mparams.Fp)*((mparams.Mp)\rp));
end
%% velocity solve
rv = rv-(aparams.B')*zp;
zv = aparams.F \ rv;
y = [zv;zp];
return

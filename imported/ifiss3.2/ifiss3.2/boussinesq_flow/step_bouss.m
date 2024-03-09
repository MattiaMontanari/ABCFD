%STEP_BOUSS set up Boussinesq flow on backward step domain 
%  step_bouss;
%
%  Grid defining data is saved to two files: step_grid1h, step_bouss_nobc
%   IFISS scriptfile: DJS; 29 April 2012.
% Copyright (c) 2012 D.J. Silvester, M.D. Mihajlovic.

%% Define domain geometry 
%
fprintf('\n--------------------------------------------------------');
fprintf('\nBoussinesq flow over a backward facing step');
fprintf('\n--------------------------------------------------------\n');
%
%% Generate the grid hierarchy
%
bouss_stepdomain 
load step_grid1h.mat 
pde=15; domain=3;   % the Boussinesq equation in step domain
%
%% Set up the problem parameters
%
fprintf('\nSPACE DISCRETISATION OF THE BOUSSINESQ PROBLEM (u,p,T): \n');
qmethod=12;  %% Q2-Q1-Q2 
fprintf('Setting up Q2-Q1-Q2 matrices:  \n');
nngpt=9;
%
%  Copy the geometry data from the grid data structure
%
lh=1;
   x=grid(lh).x;
   y=grid(lh).y;
   xy2=grid(lh).xy2;
   xy1=grid(lh).xy1;
   mv2=grid(lh).mv2;
   mp1=grid(lh).mp1;
   bnd_d   =grid(lh).bnd_d;
   bnd_dv2 =grid(lh).bnd_dv2;
   mt2     =grid(lh).mt2;
   bnd_dn2 =grid(lh).bnd_dn2;
   bnd_dnt2=grid(lh).bnd_dnt2;
         [x,y,xyv,xyp,xyt]=q2q1q2grid(x,y,xy2,xy1,mv2,mp1,mt2,bnd_d,bnd_dn2,...
                                      bnd_dv2,bnd_dnt2,H,L);
         [Av,Qv,B,Bx,By,Qp,M,At,Qt,BBx,BBy,f,g,h]= ...
                            bouss_q2q1q2(xyv,xyp,xyt,mv2,mp1,mt2,nngpt,lh);
%
%% Store the matrices into the data structure
%
   spmat(lh).Av=Av;
   spmat(lh).Qv=Qv;
   spmat(lh).At=At;
   spmat(lh).Qt=Qt;
   spmat(lh).B=B;
   spmat(lh).Bx=Bx;
   spmat(lh).By=By;
   spmat(lh).Qp=Qp;
   spmat(lh).M=M;
   spmat(lh).f=f;
   spmat(lh).g=g;
   spmat(lh).h=h;
   spmat(lh).BBx=BBx;
   spmat(lh).BBy=BBy;
%
%% Store grid data into the data structure
   grid(lh).xyv=xyv;
   grid(lh).xyp=xyp;
   grid(lh).xyt=xyt;
%
gohome; cd datafiles;
save step_grid1h.mat grid str H L gh nh pde domain
save step_bouss_nobc.mat spmat qmethod x y pde domain
fprintf('Grid data re-saved in step_grid1h.mat.\n')
fprintf('System matrices saved in step_bouss_nobc.mat.\n')

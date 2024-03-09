function box_unsteadyflowref(qmethod,ev,sol,tt,A,By,Bx,G,xy,xyp,x,y,bound,hty,snaptime)
%BOX_UNSTEADYFLOWREF plots cavity flow data at snapshot times 
%   box_unsteadyflowref(qmethod,mv,sol,tt,A,By,Bx,G,xy,xyp,x,y,bound,hty,snaptime);
%   input
%          qmethod    mixed method 
%          ev         mv/ev  Q2/Q1 element mapping matrix
%          sol        flow solution vector
%          tt         snapshot time vector
%          A          vector diffusion matrix
%          By         velocity  y-derivative matrix    
%          Bx         velocity x-derivative matrix    
%          G          veclocity mass matrix
%          xy         velocity nodal coordinate vector  
%          xyp        pressure nodal coordinate vector  
%          x          vector of x-axis interpolation points
%          y          vector of y-axis interpolation points
%          bound      boundary vertex vector
%          hty        horizontal/vertical (1/2) hot wall switch 
%          snaptime   vector of snapshot time step levels
%
% calls function xxstreambc.m to set boundary values
% calls function boxx.m to color code imposed temperature
%   IFISS function: DJS; 3 May 2012.
% Copyright (c) 2011 D.J. Silvester, H.C. Elman, A. Ramage 
fprintf('\n   Plotting flow field snapshots ... ')
L=max(x); H=max(y); nstep=length(snaptime);
if nstep>9, error('Too many snapshots!'), end
nvtx=length(xy); nu=2*nvtx; np=length(xyp);
[LG,UG]=lu(G(1:nvtx,1:nvtx)); 
Asv=A(1:nvtx,1:nvtx); fzero=zeros(nvtx,1);
[Abc,fzero]=streambc(Asv,fzero,xy,bound);
[LA,UA]=lu(Abc); 
%
fprintf('\n   step   time    mean_vorticity    min_phi  max_phi\n')
%
% ------------------ loop over snapshots
for k=1:nstep
kk=snaptime(k); ttk=tt(kk);
% compute derived quantites
u=sol(:,kk);
ux=u(1:nvtx); uy=u(nvtx+1:nu);  utotal=sqrt(ux.*ux+ uy.*uy);
fsv=-[By,-Bx]*u;
omega=UG\(LG\fsv);
f=[By,-Bx]*u;
[fsv]=xxstreambc(Asv,f,xy,bound,ttk);   
phi=UA\(LA\fsv);  
if qmethod > 1, wev = vorticity_q2(xy,ev,omega,0);
else, wev = vorticity_q1(xy,ev,omega,0); end
%
% interpolate to a cartesian product mesh
[X,Y]=meshgrid(x,y);
%
% plot stream function
figure(101)
indx=100*nstep +10 +k;
subplot(indx)
xysol = griddata(xy(:,1),xy(:,2),phi,X,Y);
maxphi=max(max(xysol)); minphi=min(min(xysol));
fprintf('  %4i  %7.3f  %11.3e   %12.5f  %9.3e\n', ...
            kk, ttk,  sum(wev), min(phi), max(phi));
vneg=[minphi:-minphi/24:0];
vpos=[maxphi/6:maxphi/6:maxphi];
vpospos=[0: maxphi/48:maxphi/12];
contour(X,Y,xysol,[vneg,vpos,vpospos])
title(['Streamlines: time = ',num2str(ttk,'%5.2f')],'FontSize',12), 
axis('equal'), boxx, axis('off')
%
% plot vorticity
figure(102)
subplot(indx)
xysol = griddata(xy(:,1),xy(:,2),omega,X,Y);
solheight = max(max(xysol))-min(min(xysol));
contour(X,Y,xysol,20)
axis('equal'), boxx, axis('off')
title(['Vorticity: time = ',num2str(ttk,'%5.2f')],'FontSize',12),  
end
% ------------------ end loop over snapshots
%
fprintf('   All done\n')
return
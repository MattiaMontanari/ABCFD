%IFISS returns IFISS version number
%IFISS scriptfile: DJS; 24 June 2012. 
%Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage (see readme.m)
%
%For help, type 'helpme'.
fprintf('This is IFISS version 3.2: released on 24 June 2012\n') 
fprintf('For help, type "helpme".\n') 
% display IFISS logo
gohome, cd plotfiles
ldata = imread('ifisslogo.jpg');
image(ldata),axis off, axis tight, 
title('IFISS Toolbox', 'FontSize', 12), 
set(gcf,'Position',[0,800,300,250]);
gohome

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \    SUB-ROUTINE PLOT A SINGLE ELEMENT AND ITS      %
%  /----\ |  \|    |--  |   |   DOFs IN LOCAL COORDINATES                      %
% /      \|__/ \__ |    |__/                                                   %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%abCFD_PlotElement( 'Y' , coord ) plots in new window a single element 
%   with local 2D coordinates defined into coord
%abCFD_PlotElement( 'Y' , coord , FC) plots into the FCth windows
 
% EXAMPLE:
%           abCFD_PlotElement(  'Y' , elem.localCoor , 1 )

function abCFD_PlotElement( Anserw , coord ,varargin )
%% PLOT ELEMENT

if Anserw == 'Y'
    % Select which windows should be used
    nVarargs = nargin - 2;
    if nVarargs == 0
        figure();
        FC_tag = gcf;
    elseif nVarargs == 1
        FC_tag = varargin{1};
    end
    
    % Plotting
	figure( FC_tag )
	set(gca,'LineStyleOrder', '-*|:|o')
	plot(coord(1,:),coord(2,:),'*r')
	title(['ELEMENT with ',num2str(size(coord,2)),' local dofs'])
	axis equal
 
end



end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  FEBRUARY 2013          % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick-off                                                  07/02/2013 %     
%  
% ---------------------------------------------------------------------------- %
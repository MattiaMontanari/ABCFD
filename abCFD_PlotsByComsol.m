% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	IMPORT COMSOL's PLOTS INTO MATLAB's FIGURES    %
%  /----\ |  \|    |--  |   |   THE USER WILL BE ASKED TO CHOOSE AMONG THREE   %
% /      \|__/ \__ |    |__/    TYPES OF PLOT                                  %    
%                               * * * CALLS * * *                              %
%                                   i. mphgeom                                 %
%                                  ii. mphmesh                                 %
%                                 iii. mphplot                                 %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%abCFD_PlotsByComsol( model , A )
% INPUTS
%   - model : common COMSOL structure
%   - A     : Trigger to access plot options
%           IF == 'Y' goes to the plot session
%           ELSE -> Skip any plot

function abCFD_PlotsByComsol( model , A )
%% Choose what to plot
if A == 'Y'
    plot_answer = input( ['\nWhat to plot? [ {0} / 1 / 2 / 3 ] \n',...
        '0)Nothing  \n1)Geometry  \n2)Solution  \n3)BCs \n'],'s' );

    if strcmp( plot_answer , '1') == 1 
        figure()
        mphgeom(model,'geom1','Edgelabels','on','Edgelabelscolor','y')
    elseif strcmp( plot_answer , '2') == 1
        % Heat flux plots
        figure()
        subplot(2 , 3, [1,4] ); mphmesh(model, 'mesh1','Edgemode','off','edgecolor',...
                                  'b','facemode','on','facealpha',0.5,'meshcolor','r');
                                hold on;
                                mphplot(model,'temperature','rangenum',1);
        subplot(2 , 3, [2,5] ); mphplot(model,'Velocity','rangenum',1);
                                legend('off')
        subplot(2 , 3, 3); mphplot(model,'Pressure','rangenum',1);
                                axis tight
        title('Heat Flux normal to \Gamma_2')
        legend('a','b','c','d','a','b','c','t = 0','t = tfinal')
        subplot(2 , 3, 6); %mphplot(model,'pg4','rangenum',1);
                                axis tight
        title('Temperature profile on \Gamma_1')
        legend('a','b','c','d','a','b','c','t = 0','t = tfinal')
    elseif strcmp( plot_answer , '3') == 1
        %  Tempearture plot                        
        figure(); colordef white
        subplot( 2,2,1)
        mphplot(model,'pg4');      
        grid on; axis tight;
        legend('a','b','c','d','a','b','c','t = 0','t = tfinal')
        title('Temperature profile on \Gamma_1')
                
        subplot( 2,2,2)
        mphplot(model,'pg3');    
        grid on; axis tight;
        legend('a','b','c','d','a','b','c','t = 0','t = tfinal')
        title('Heat Flux normal to \Gamma_2')
        
      	subplot( 2,2,3)
        mphplot(model,'pg6');    
        grid on; axis tight;
        legend('a','b','c','d','a','b','c','t = 0','t = tfinal')
        title('Heat Flux normal to \Gamma_3')
        
        subplot( 2,2,4)
        mphplot(model,'pg5');    
        grid on; axis tight;
        legend('a','b','c','d','a','b','c','t = 0','t = tfinal')
        title('Heat Flux normal to \Gamma_4')
    end
end
 
%% CLOSURE
fprintf('Plots from Comsol imported at  %s \n',  datestr(now,15) )

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.5                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.5 - Few comments added. Should be adpted for new plots        27/02/2013 %
%   0.4 - Updated to the new bonudaries' numbering. Plots with      13/02/2013 %
%       choice 3 and 4 are validated, however: in 3 the time is missing in     %
%       the first two plots, and the legens (in 3 and 4) are not completed.    %
%   0.3 - CHOICE OF NOT TO PLOT OR PLOT ONLY A SINGLE FIGURE        06/02/2013 %
%   0.2 - ADDED A TEMPERATURE PLOT ALONG SIGMA 1 AND SIGMA 4        05/02/2013 %
%         THIS INCLUDES TWO LINES TO INSPECT THE TEMPERATURE AND TO            %
%         ENSURE THE DIRICLET BC'S IS APPLIED CORRECTLY                        %
%   0.1 - kick-off                                                  05/02/2013 %     
%  
% ---------------------------------------------------------------------------- %
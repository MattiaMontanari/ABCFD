% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \  	SUB-ROUTINE TO COMPUTE AND ORDER THE           %
%  /----\ |  \|    |--  |   |   EMPIRICAL MODES OF THE REDUCED BASIS.          %
% /      \|__/ \__ |    |__/    PROPER ORGHOGONAL DECOMPOSITION METHOD         %                     
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[U]=abCFD_Modes( U , ee , A) 

function [ U ] = abCFD_Modes( U , ee , fun )

%% INITIALIZE
minE = 3;
if ee < minE
    ee = minE;
end

%% COMPUTE TEMPERATURE MODES
% Compute modes and store in the U.T structure
tiin = toc;
fprintf('Start eigenvalue problem at %s  \n', datestr(now,15) )

[ U.T.eigval , U.T.B, U.T.energy ] = abCFD_EigProblem( U.T , ee, fun );
fprintf('Temperature modes computed in %f minutes \n', (toc-tiin)/60 )

%% COMPUTE PRESSURE MODES
% Not required!
% [ U.p.eigval , U.p.B ] = abCFD_EigProblem( U.p , ee, fun );

%% COMPUTE X and Y VELOCITY COMPONENTS MODES
% Compute and store modes
tiin = toc;
[ U.u.eigval , U.u.B, U.u.energy] = abCFD_EigProblem( U.u , ee, fun );
fprintf('U-component modes computed in %f minutes \n', (toc-tiin)/60 )

tiin = toc;
[ U.v.eigval , U.v.B, U.u.energy ] = abCFD_EigProblem( U.v , ee, fun );
fprintf('V-component modes computed in %f minutes \n\n', (toc-tiin)/60 )

% %% Plots  * ** *  * * * * * * * N O T    W O R K I N G * * * * * * * * * * *
% %                   issue due to the new numbering. the actual plot
% %                   sequence uses abCFD_patch, thur requires faces, points
% %                   etc...
% %                   In the new numbering this is not provided yet
% % 
% % TEMPERATURE
% figure;
% if strcmp( A(1) , 'T') == 1
%     % Plot Energy per modes rate
%     subplot( 2, 2 , 1);
%     bar( T_energy )
%     title('Energy per modes rate')
%     
%     % Plot each single mode
%     for id_e = 1 : minE - 1
%         subplot( 2, 2 , 1 + id_e );
%         abCFD_streamline( U, U.T.B' , U.T.B' , id_e );
%         view( 2 )
%         abCFD_patch( U.T , id_e, U.T.B)
%         title( ['\fontsize{14}Eigenmode #', num2str(id_e) ,...
%                             ' for quantity: ', char( U.T.expr )])
%     end
% end
% 
% % % VELOCITY
% % figure;
% % if strcmp( A(3:4) , 'UV') == 1
% %     % Plot Energy U 
% %     subplot( 3, 2 , 1);
% %     bar( U_energy )
% %     title('Energy per U_modes rate')
% %     % Plot Energy V 
% %     subplot( 3, 2 , 2);
% %     bar( V_energy )
% %     title('Energy per V_modes rate')
% %     
% %     % Plot each single mode
% %     for id_e = 1 : minE
% %         subplot( 3, 2 , 2 + id_e );
% %         abCFD_streamline( U.u, U.u.B' , U.v.B' , id_e );
% %         
% %         view( 2 )
% %         title( ['\fontsize{14}Eigenmode #', num2str(id_e) ,' for quantity: ', char(U.v.expr)])
% %     end
% % end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.3                                 date:  April 2013             % 
% ---                                                                      --- %
%   0.3 - Removed plot sections                                     15/04/2013 %  
%   0.2 - New Concept. This routine collects subroutines to compute 01/03/2013 %  
%         and plot modes separately. The computed modes is minium minE         %
%         but can be specified byy ee if > 3. However the plots concern        %
%         only first minE modes.                                               %
%   1.0 - kick-off                                                  28/02/2013 %      
% ---------------------------------------------------------------------------- %
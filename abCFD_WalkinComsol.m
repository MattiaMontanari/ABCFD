% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \  	ADRESSEES A NEW MODEL TO BE SOLVED IN COMSOL   %
%  /----\ |  \|    |--  |   |   OR LOADS A MODEL, AND ASSEMBLE FILE, STORED IN %
% /      \|__/ \__ |    |__/    CURRENT MATLAB's FOLDER                        %
%                                                                              %
% * * * CALLS * * *                                                            %
%               i.abCFD_mdlobj_FEM                                            % 
%              ii.mphload                                                      %      
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ model, assemb ] = abCFD_WalkinComsol( ph, A_SaveLoad, ModelVersion )
% INPUTS:
%   - A_SaveLoad:   Defines if to load or to start a model
%                   IF == 'S' saves&solves a new model
%                   IF == 'L' loads an old model, 
%                   ELSE  -> solves the model without saving it
%   - ph        :   common phase structure
%   - ModelVersion: string for model identification
% OUTPUT:
%   - models    :   common COMSOL's structure
%   - assemb    :   common assembly structure


function [ model, assemb ] = abCFD_WalkinComsol( ph, A_SaveLoad, ModelVersion )

    if strcmp( A_SaveLoad, 'S' ) == 1
        % RUN AND SAVE A NEW MODEL
        [ model, assemb ] = abCFD_mdlobj_FEM( ph,ModelVersion  );
        save('input_assemb.mat', 'assemb');
        mphsave(model,['heat_mode14_',ModelVersion,'.mph'])
        disp(['NEW MODEL SAVED IN: ',cd,'heat_mode14_',ModelVersion,'.mph']);
        
    elseif strcmp( A_SaveLoad, 'L' ) == 1
        
        % LOAD A OLD MODEL
        disp(['LOAD MODEL: ',cd,'heat_mode14_',ModelVersion,'.mph']);
        [ model ] = mphload( ['heat_mode14_',ModelVersion,'.mph'] );    
        load('input_assemb.mat', 'assemb');
        fprintf('Model loaded after %f sec. \n',toc)
    else
        
        % RUN A NEW MODEL
        [ model, assemb ] = abCFD_mdlobj_FEM( ph,ModelVersion  );    
        
    end

end
% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.2                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.2 - General revision. Few comments added                      27/03/2013 %
%   0.1 - kick-off                                                  27/02/2013 %
% ---------------------------------------------------------------------------- %
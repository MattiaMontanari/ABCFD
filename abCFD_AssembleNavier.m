% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \   	ASSEMBLE x x x x x  MATRIX GIVEN A SOLUTION    %
%  /----\ |  \|    |--  |   |   VECTOR. LINEAR TRIANGULAR ELEMENTS ONLY        %
% /      \|__/ \__ |    |__/                                                   %
%                                                                              %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%[ K, M] = abCFD_AssembleNavier( var, TRI, SOL )
% INPUTS:
%   - var   : Usual structure
%   - TRI   : Connectivity matrix for linear triangular elements
%   - SOL   : Solution vector 'included' into the matrix
% OUTPUT:

function [ K11, K12, K21, K22, M] = abCFD_AssembleNavier( var, TRI  )

%% I N I T I A L I Z E
% check coordinates
% if size( XY, 2 ) > size( XY, 1 )
%     XY = XY';
% end
X = var.dofCoord( 1 , :)';
Y = var.dofCoord( 2 , :)';

IntPrecision = 2;                % Numerical integration order
NumDof = sum(var.fieldNDofs); 
[ NumEle, ShpOrd ] = size( TRI );    % Number elements , number of dofs per element

if ShpOrd == 3
%% VELOCITY MATRICES
    % Initialize unkwon matrix
    K11 = sparse( NumDof, NumDof);
    K12 = sparse( NumDof, NumDof);
    K21 = sparse( NumDof, NumDof);
    K22 = sparse( NumDof, NumDof);
    M = sparse( NumDof, NumDof);
    % Get Gauss poits and weights
    [ sGP , tGP, WGP ] = abCFD_GaussQuad2D( 'tri' , IntPrecision );

    % inner loop over elements  
    Xloc = zeros( NumEle, ShpOrd);
    Yloc = zeros( NumEle, ShpOrd);
%     Uloc = zeros( NumEle, ShpOrd);
%     Vloc = zeros( NumEle, ShpOrd);
    
    for i_el = 1 : ShpOrd
       Xloc( : , i_el ) = X( TRI( : ,  i_el) );
       Yloc( : , i_el ) = Y( TRI( : ,  i_el) );
%        Uloc( : , i_el ) = SOL( TRI( : ,  i_el) , 1);
%        Vloc( : , i_el ) = SOL( TRI( : ,  i_el) , 1);
    end
    
    % initialize local matrix
    Kele11 = zeros( NumEle, 3,3);
    Kele12 = zeros( NumEle, 3,3);
    Kele21 = zeros( NumEle, 3,3);
    Kele22 = zeros( NumEle, 3,3);
    Mele = zeros( NumEle, 3,3);
    % Loop over Gauss points
    for i_GP = 1 : numel( WGP )
       
        % Shape functions are theirs derivatives
        phi_e(1) =   1    - sGP( i_GP )   - tGP( i_GP );
        phi_e(2) =        + sGP( i_GP )                ;
        phi_e(3) =                        + tGP( i_GP );

        dphids(1) = - 1 ;
        dphids(2) = + 1 ;
        dphids(3) =   0 ;
      
        dphidt(1) = - 1 ;
        dphidt(2) =  0  ;
        dphidt(3) = + 1 ;
        
        % Evaluate Jabocian - - - - - - - - - - - - - - - - - - - - - - -  
        % Initialize vectors
       	dxds = zeros( NumEle , 1);
        dxdt = zeros( NumEle , 1);
        dyds = zeros( NumEle , 1);
        dydt = zeros( NumEle , 1);
        jac = zeros( NumEle , 1);
        invjac = zeros( NumEle , 1); 
        
        for i_el = 1 : ShpOrd
            dxds(:) = dxds(:) + Xloc(:,i_el) .* ones( NumEle, 1 ) * dphids(i_el);
            dxdt(:) = dxdt(:) + Xloc(:,i_el) .* ones( NumEle, 1 ) * dphidt(i_el);
            dyds(:) = dyds(:) + Yloc(:,i_el) .* ones( NumEle, 1 ) * dphids(i_el);
            dydt(:) = dydt(:) + Yloc(:,i_el) .* ones( NumEle, 1 ) * dphidt(i_el);
        end
        % compute Jacobian
        jac(:) = dxds(:).*dydt(:) - dxdt(:).*dyds(:);
        % Compute Jacobian inverse
        invjac(:) = ones( NumEle, 1 ) ./ jac(:);
        % compute x derivative of phi and y derivative of phi 
        for i_el = 1 : ShpOrd
            phi(:,i_el) = phi_e(i_el)*ones( NumEle, 1 );
            dphidx(:,i_el) = ( dphids(:,i_el).*dydt(:) - dphidt(:,i_el).*dyds(:));
            dphidy(:,i_el) = (-dphids(:,i_el).*dxdt(:) + dphidt(:,i_el).*dxds(:));        
        end
        % Compute velocity derivatives
%         u_x = zeros( NumEle ,1); 
%         u_y = zeros( NumEle ,1);
%         for i_el = 1 : ShpOrd
% 				u_x(:) = u_x(:) + Uloc(:, i_el ) .* phi(:, i_el );
% 				u_y(:) = u_y(:) + Yloc(:, i_el ) .* phi(:, i_el );	 
%         end
        
        % Assemble matrix
        for j = 1 : ShpOrd
            for i = 1 : ShpOrd              
                                            
                    Kele11(:,i,j)  = Kele11(:,i,j)  + WGP(i_GP) * dphidx(:,i).*dphidx(:,j).*invjac(:) ;
                    Kele12(:,i,j)  = Kele12(:,i,j)  + WGP(i_GP) * dphidy(:,i).*dphidx(:,j).*invjac(:) ;
                    Kele21(:,i,j)  = Kele21(:,i,j)  + WGP(i_GP) * dphidx(:,i).*dphidy(:,j).*invjac(:) ;
                    Kele22(:,i,j)  = Kele22(:,i,j)  + WGP(i_GP) * dphidy(:,i).*dphidy(:,j).*invjac(:) ;
                    
                    Mele(:,i,j)  = Mele(:,i,j)  + WGP(i_GP) * phi(:,i).*phi(:,j).*jac(:) ;
            end
        end        
    end
    
    %  element assembly into global matrix
    for krow = 1 : 3
        nrow = TRI(:,krow);	 
          for kcol = 1 : 3
              ncol = TRI(:,kcol);	  
              K11 = K11 + sparse( nrow , ncol , Kele11(:,krow,kcol) , NumDof , NumDof);
              K12 = K12 + sparse( nrow , ncol , Kele12(:,krow,kcol) , NumDof , NumDof);
              K21 = K21 + sparse( nrow , ncol , Kele21(:,krow,kcol) , NumDof , NumDof);
              K22 = K22 + sparse( nrow , ncol , Kele22(:,krow,kcol) , NumDof , NumDof);
              M = M + sparse( nrow , ncol , Mele(:,krow,kcol) , NumDof , NumDof);
          end
    end
    
else
    % Other convective matrices
end
end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  April 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick off                                                  15/04/2013 %
% ---------------------------------------------------------------------------- %
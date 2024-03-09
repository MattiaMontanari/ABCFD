% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %
%   /  \  |__  /   |    |  \   	COMPUTES THE PARTIAL DERIVATES OF A FUNCTION   %
%  /----\ |  \|    |--  |   |   EXPRESSED AS POLYNOMIAL WITH N-INDEPENDENT     %
% /      \|__/ \__ |    |__/    VARIABLES.                                     %
%                                                                              %
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%PartDerivat = abCFD_PolynDer ( p ) 


function PartDerivat = abCFD_PolynDer ( p )
% Number of independent varibles
IndVar_Num = length( p.VarNames);
% Number of terms/coeff of the polynomial
Terms_Num = length( p.Coefficients );


% Initialize space for partial derivatives
PartDerivat = cell( 1, IndVar_Num );
for i_var = 1 : IndVar_Num
    
    Coef = zeros( 1, Terms_Num );
    Expo = zeros( Terms_Num , IndVar_Num );
    % Loop over polynomial's term
    for i_term = 1 : Terms_Num
       
         
        % Check exponent
        if p.ModelTerms( i_term , i_var) == 0
            % Derivative is null
            Coef( i_term ) = 0;
        elseif p.ModelTerms( i_term , i_var) > 0
            % Derivative is not null and the coeffcient is:
            Coef( i_term ) = p.ModelTerms( i_term , i_var)* p.Coefficients( i_term );
            % While the new exponent is
            Expo(  i_term ,  : ) = p.ModelTerms( i_term , :) ;
            Expo(  i_term , i_var ) = p.ModelTerms( i_term , i_var ) - 1;
        end
        
    end
    
    % Store exponents and coefficients
    PartDerivat{ 1, i_var } = struct( 'Coefficients', Coef,...
                                      'ModelTerms', Expo  );
%                                       'VarNames', p.VarNames );
    
end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   revision hystory:                                               -  date -  %
%   0.1 - kick off                                                  22/03/2013 %
% ---------------------------------------------------------------------------- %
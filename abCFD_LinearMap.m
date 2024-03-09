% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %%
%    /\   |     __  ___  __                                                    %            
%   /  \  |__  /   |    |  \    APPLY A FUNCTION f(x1,..xN) TO THE SUB-ARRAY   %
%  /----\ |  \|    |--  |   |   OF A BIGGER MATRIX, MAPPING THE MODIFED        %
% /      \|__/ \__ |    |__/    ENTRIES IN PLACE OF THE ORGINAL ENTRIES.       %
%                               THE SUB-MATRIX CAN BE RECTANGULAR              %
% f(x1,..xN) IS FUNCTION OF N-VARIABLES
%  <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>    ><    <>  %
%BIG = abCFD_LinearMap( BIG, RowsTag, ColsTag, Fun, )
% Returns the array 'BIG' with the entries BIG(RowsTag, ColsTag) modified
% INPUT:
%       - BIG       : Original matrix to be modified and returned
%       - RowsTag   : Rows' tags intended to be modifed
%       - ColsTag   : Columns' tags intended to be modifed
%       - Fun       : function handle of N-variables
%       - {FunVar}  : if N variables must be specified
%                   if {FunVar} == [] the unique {FunVar} is W = BIG(RowsTag, ColsTag)

function BIG = abCFD_LinearMap( BIG, RowsTag, ColsTag, Fun, varargin)

%% Initialize

% working matrix
W = BIG( RowsTag, ColsTag );

%% OPERATION
switch numel( varargin )
    case 0
        W = Fun( W );
    case 1
        W = Fun( W , varargin{1} );
    case 2
        W = Fun( W , varargin{1} ,varargin{2} );
    case 3
        W = Fun( W , varargin{1} ,varargin{2} ,varargin{3} );
end

%% MAPPING TO 'BIG'

BIG( RowsTag, ColsTag ) = W;

end

% ---------------------------------------------------------------------------- %
%   Author: MATTIA MONTANARI         mattia.montanari@eleves.ec-nantes.fr      % 
% ---                                                                      --- %
%   Version: 0.1                                 date:  MARCH 2013             % 
% ---                                                                      --- %
%   0.1 - kick-off                                                  21/03/2013 %
% ---------------------------------------------------------------------------- %
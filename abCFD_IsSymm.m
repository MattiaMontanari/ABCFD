function [ issym ] = abCFD_IsSymm( A )

if isequal(A,A.') == 1
    
    issym = 'Symmetric';
    
elseif isequal(A,A.') == 0
    
    issym = 'Non symm.';
    
end
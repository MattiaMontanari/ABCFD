function abCFD_RelativeDiff( A, B )

if issparse( A )
    A = full( A );
end


if issparse( B )
    B = full( B );
end

dif = arrayfun( @(d,a) 100*d/a, A-B, B);

MIN = -1* full(max( max( -dif ) ));
MAX = full(max( max( dif ) ));
sprintf('%% Error: MIN = %d MAX %d\n', round(abs(MIN)), round(abs(MAX)) ) 


%% CHECK PATTER
Ar = abCFD_round(  1e-9, A);
Br = abCFD_round(  1e-9, B);
% A(A == 0)=[];
A(A ~= 0)=1;
% B(B == 0)=[];
B(B ~=0 )=1;

if full( sum(sum( A-B ))) == 0
    disp( '  * * * same patter * * * ')
else
    figure
    spy(Ar)
    hold on
    spy(Br,'ro')
    spy(Ar-Br,'gx')
end

end
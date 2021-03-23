function A = mlm_testmats_mct(k, n)
% Set of 32 matrices that are used in testing the algorithms for computing
% the matrix Mittag-Leffler functions (with circul frank magic removed 
% since they cause overflow problems).

matrices = [       
    'chebspec'; 'chebvand'; 'chow    ';
    'clement '; 'cycol   '; 'dramadah'; 'forsythe'; 
                'gearmat '; 'grcar   '; 'hanowa  '; 
    'invhess '; 'jordbloc'; 'kahan   '; 'leslie  '; 
    'lesp    '; 'lotkin  '; 'parter  '; 'randcolu'; 
    'rando   '; 'randsvd '; 'redheff '; 'riemann '; 
    'sampling'; 'smoke   '; 'triw    '; 'orthog  ']; 
n_gall = size(matrices,1);

% Matrices from Matrix Computation Toolbox:
matrices = [matrices;
    'gfpp    '; 'makejcf '; 'rschur  '; 'vand    ';];

% Other MATLAB matrices:
matrices = [matrices;
                'rand    '; 'randn   ';];

if nargin == 1
    if k == 0
        A = size(matrices,1);
    elseif k > 0
        A = deblank(matrices(k,:));
    else
        A = 'Modified on Jul 31 2020';
    end
else
    if k < n_gall
        A = eval( ['gallery(''' deblank(matrices(k,:)) ''',n)'] );
    elseif k == n_gall % gallery('orthog',n,-2);
        A = eval( ['gallery(''' deblank(matrices(k,:)) ''',n, -2)'] );
    else
        A = eval( [deblank(matrices(k,:)) '(n)'] );
    end
end
end
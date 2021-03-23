function A = testmats_mct(k, n)
% Set of 35 matrices that are used in testing the algorithms for computing
% matrix functions in section 5.

matrices = [       
    'chebspec'; 'chebvand'; 'chow    '; 'circul  ';
    'clement '; 'cycol   '; 'dramadah'; 'forsythe'; 
    'frank   '; 'gearmat '; 'grcar   '; 'hanowa  '; 
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
    'magic   '; 'rand    '; 'randn   ';];

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
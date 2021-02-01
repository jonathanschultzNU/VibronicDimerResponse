function output = vibron_double_exc(v_max)

%	Author: Jonathan D. Schultz
%	Email: jonathanschultz2022@u.northwestern.edu
%	Last revision date: February 1st, 2021
%
%	Copyright: Jonathan D. Schultz, 2021

%   Please see readme file for information about this package

length = (v_max^2 + 3*v_max+2)/2;
basis = zeros(length,4);

k = 1;
for v_1 = 0:v_max
    for v_2 = 0:v_max-v_1
        basis(k,:) = [1, v_1, 1, v_2];
        k = k+1;
    end
end

output = basis;


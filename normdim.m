function out = normdim(data)

%	Author: Jonathan D. Schultz
%	Email: jonathanschultz2022@u.northwestern.edu
%	Last revision date: February 1st, 2021
%
%	Copyright: Jonathan D. Schultz, 2021

%   Please see readme file for information about this package

d = ndims(data);

if d == 1
    out = data./max(abs(data));
end

if d == 2
    out = data./max(max((abs(data))));
end

end
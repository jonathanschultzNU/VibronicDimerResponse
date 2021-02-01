function c = elcre(basis,n_states,molnum)

%	Author: Jonathan D. Schultz
%	Email: jonathanschultz2022@u.northwestern.edu
%	Last revision date: February 1st, 2021
%
%	Copyright: Jonathan D. Schultz, 2021

%   Please see readme file for information about this package

molnum = 2*(molnum-1)+1;
c = zeros(n_states);

for i = 1:n_states
    for j = 1:n_states
        bra = basis(i,:);
        ket = basis(j,:);
        
        if ket(molnum) > 0
           ket(molnum) = ket(molnum)-1; 
           if bra == ket
              c(i,j) = 1;
           end  
        end
        
    end
end 

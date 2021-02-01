function b = vibcre(basis,n_states,molnum)

%	Author: Jonathan D. Schultz
%	Email: jonathanschultz2022@u.northwestern.edu
%	Last revision date: February 1st, 2021
%
%	Copyright: Jonathan D. Schultz, 2021

%   Please see readme file for information about this package

b = zeros(n_states);

for i = 1:n_states
    for j = 1:n_states
        bra = basis(i,:);
        ket = basis(j,:);
        
        if ket(molnum*2) > 0
           ket(molnum*2) = ket(molnum*2)-1; 

           if bra == ket
              b(i,j) = sqrt(ket(molnum*2)+1);      
           end  
        end    
    end
end
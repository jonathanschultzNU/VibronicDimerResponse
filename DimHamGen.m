function Hstruct = DimHamGen(Par)

%	Author: Jonathan D. Schultz
%	Email: jonathanschultz2022@u.northwestern.edu
%	Last revision date: February 1st, 2021
%
%	Copyright: Jonathan D. Schultz, 2021

%   Please see readme file for information about this package

n_states = 2*(Par.v_max^2 + 3*Par.v_max+2);         % number of states
corr = -Par.e1*eye(n_states*0.25);                  % Removal of fast-oscillating frequencies from diagonal elements of Hamiltonian - added back to Fourier Transform axis before plotting
fock = eye(n_states);                               % define Fock space (site basis vectors)
lambda_sqr = (Par.lambda^2)*eye(n_states);          % Matrix of Huang-Rhys factor squared for Hamiltonian term

basis_g = vibron_zero_exc(Par.v_max);               % generate ground-state basis set
basis_s = vibron_single_exc(Par.v_max);             % generate singly-excited basis set
basis_d = vibron_double_exc(Par.v_max);             % generate doubly-excited basis set
basis = [basis_g; basis_s; basis_d];                % combine into total basis set

%construct vibrational annhiliation operators
b1 = vibcre(basis,n_states,1);
b2 = vibcre(basis,n_states,2);

% construct electronic operators   
c1 = elcre(basis,n_states,1);
c2 = elcre(basis,n_states,2);

%Build the Hamiltonian
H1 = Par.e1*(c1')*c1;
H2 = Par.e2*(c2')*c2;
H3 = Par.J.*((c1')*c2+(c2')*c1);
H4 = Par.w0*(b1')*b1;
H5 = Par.w0*(b2')*b2;
H6 = Par.w0*(c1')*c1*(Par.lambda.*(b1'+b1)+lambda_sqr);
H7 = Par.w0*(c2')*c2*(Par.lambda.*(b2'+b2)+lambda_sqr);
H = H1+H2+H3+H4+H5+H6+H7; 

%construct transition dipole operators
mu1E = 1;   
mu2E = 1;   %set aribrary scaling factor for influence of electric field 
mu = mu1E*(c1'+c1) + mu2E*(c2'+c2);
mu = mu./2;
mu = mu*Par.c;  %unit conversion

%Rotate all matrices and vectors to eigenbasis of Hamiltonian
[V,D] = eig(H);       
A = zeros(4); A(2,2)=1; A(3,3)=1; A(4,4)=1; % state which blocks need to be anti-aliased 
corr = kron(A,corr);                        % create correction matrix to avoid aliasing signals
Hwvn = H;
H = (H+corr).*Par.c;                        % anti-aliasing and unit conversion

BC = @(A) V'*A*V;       % Define change of basis function
fock_e = BC(fock);      % Rotate basis
mueig = BC(mu);

n = length(H)/4;        % Define blocks that will be separated into individual matrices

% Separate the Hamiltonian and transition dipole operator into various state
% components. This is done solely to decrease computational cost, since we
% are avoiding propagating the wavefunction under the entire system
% Hamilonian at each step. Rather, we are propagating the system under
% whichever block of the Hamiltonian is relevant at a given time.
Hgg = H(1:n,1:n); Hggwvn = Hwvn(1:n,1:n);
Hse = H(n+1:3*n,n+1:3*n); Hsewvn = Hwvn(n+1:3*n,n+1:3*n);
Hde = H(3*n+1:end,3*n+1:end); Hdewvn = Hwvn(3*n+1:end,3*n+1:end);
fockgg = fock(1:n,1:n);
fockse = fock(n+1:3*n,n+1:3*n);
fockde = fock(3*n+1:end,3*n+1:end);
MUg_se = mu(1:n,n+1:3*n);
MUse_g = mu(n+1:3*n,1:n);
MUse_de = mu(n+1:3*n,3*n+1:end);
MUde_se = mu(3*n+1:end,n+1:3*n);

%Pack the output structure
Hstruct.Hgg = Hgg; Hstruct.Hggwvn = Hggwvn;
Hstruct.Hse = Hse; Hstruct.Hsewvn = Hsewvn;
Hstruct.Hde = Hde; Hstruct.Hdewvn = Hdewvn;
Hstruct.Hdiag = D;
Hstruct.Hwvn = Hwvn;
Hstruct.eigvec = V;
Hstruct.fockgg = fockgg;
Hstruct.fockse = fockse;
Hstruct.fockde = fockde;
Hstruct.fock_e = fock_e;
Hstruct.fock = fock;
Hstruct.MUg_se = MUg_se;
Hstruct.MUse_g = MUse_g;
Hstruct.MUse_de = MUse_de;
Hstruct.MUde_se = MUde_se;
Hstruct.MU = mu;
Hstruct.mueig = mueig;
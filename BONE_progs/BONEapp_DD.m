function mu = BONEapp_DD(a,A,B);
%% Nemat-Nasser (1993) result for a straight crack;
%% gives effective shear modulus divided by \mu_0,
%% effective shear modulus for host;

denom = 1-2*A./(pi*B).*log(cos(pi*a./(2*A)));
mu    = 1./denom;

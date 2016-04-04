function mu = BONEapp_AchLi(a,A,B);
%% Nemat-Nasser (1993) result for a straight crack;
%% gives effective shear modulus divided by \mu_0,
%% effective shear modulus for host;

beta  = 2*A./(pi*B).*log(cos(pi*a./(2*A)));
%mu    = sqrt(1-beta)./(1-2*beta);

mu    = 1./(1-beta);

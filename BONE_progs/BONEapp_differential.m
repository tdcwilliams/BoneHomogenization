function mu = BONEapp_differential(a,A,B);
%% Deng & Nemat-Nasser (1992) result for a straight crack;
%% gives effective shear modulus divided by \mu_0,
%% effective shear modulus for host;

%exponent = (8*A^2/pi./a.^2).*log(cos(pi*a/2/A));
exponent = -(8*A^2/pi./a.^2).*log(cos(pi*a/2/A));%%opposite to D&N(1992) but must be a typo;
%%
f  = a.^2/4/A/B;
mu = NaN*ones(size(f));
%%
jj       = find(f<1);
mu(jj)   = abs(1-f(jj)).^exponent(jj);

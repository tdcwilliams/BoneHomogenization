function mu = BONEapp_NemNas(a,A,B);
%% Nemat-Nasser (1993) result for a straight crack;
%% gives effective shear modulus divided by \mu_0,
%% effective shear modulus for host;

N     = 5000;
s10   = 0;
%%
for n=1:N
   v1    = n*pi*a./A;
   v2    = n*pi*B./A;
   t1    = 2*besselj(1,v1)./v1;
   t2    = v2./tanh(v2);

   %% sum defined as over \Zfield-{0}, but
   %% I sum over \Nfield (since it is even in n),
   %% so need to multiply by 2;
   s10   = s10+2*t2.*t1.^2;
end
%%
mu = 1./(1+1./s10);

function [soL,ip_stuff,Mlog,w_qd]=...
   BONE_ipmatrices_cos(Nint,Nterms)

Del=2/Nint;
soL=-1+( (0:Nint-1)'+.5 )*Del;
w_qd=Del+0*soL;
%%
kn=pi*(1:Nterms);
DD=diag(1./kn);
COS=cos(soL*kn);
SIN=sin(soL*kn);
%%
d1chi_vals=[.5+0*soL,COS,SIN];
chi_vals=[.5+0*soL,SIN*DD,-COS*DD];

ip_chi=Del*chi_vals';
ip_d1chi=Del*d1chi_vals';
hn=[2,1+0*kn,1+0*kn]';
%%
ip_stuff={{ip_d1chi,ip_chi},hn,{d1chi_vals,chi_vals}};
Mlog=-diag([0,1/2./[kn,kn]]);

if 0%% do test
  ff=cos(3*pi*soL-.57);
  fn=ip_d1chi*ff;
  ff_ap=fn(1)+d1chi_vals(:,2:end)*fn(2:end);
  %%
  plot(soL,ff), hold on;
  plot(soL,ff_ap,'--g'), hold off;
end
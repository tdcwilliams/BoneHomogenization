function [soL,ip_stuff,Mlog,w_qd]=...
   BONE_ipmatrices_exp(Nint,Nterms)

Del=2/Nint;
soL=-1+( (0:Nint-1)'+.5 )*Del;
w_qd=Del+0*soL;
%%
kn=pi*(-Nterms:Nterms);
DD=diag(kn);
Exp=exp(1i*soL*kn);
ip_chi=Del/2*Exp';
ip_d1chi=1i*DD*ip_chi;
hn=1/2+0*kn';
%%
ip_stuff={{ip_chi,ip_d1chi},hn,{Exp/2,1i*Exp*DD/2}};
kn(Nterms+1)=Inf;
Mlog=-diag(1/4./kn);

if 0%% do test
  ff=cos(3*pi*soL-.57);
  fn=ip_chi*ff;
  ff_ap=Exp*fn;
  %%
  plot(soL,ff), hold on;
  plot(soL,ff_ap,'--g'), hold off;
end
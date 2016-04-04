function FF=BONE_forcing_crk_Tn(GEOM_stuff,IP_stuff,GFxn,GF_args);

%% stuff for s integral - Chebyshev quadrature points:
geom_stuff=GEOM_stuff{1};
soL=geom_stuff{1};
xx=geom_stuff{2};
yy=geom_stuff{3};
LC=geom_stuff{4};
Xtra=geom_stuff{5};
th_vec=Xtra{1};

%% stuff for s_0 integral - Legendre quadrature points:
geom0_stuff=GEOM_stuff{2};
s0oL=geom0_stuff{1};
xx0=geom0_stuff{2};
yy0=geom0_stuff{3};
LC0=geom0_stuff{4};
Xtra0=geom0_stuff{5};
%%
d2xy_ds2=Xtra0{5};
th_vec0=Xtra0{1};
ff=-exp(1i*th_vec0);
%%
[XX0,XX]=meshgrid(xx,xx0);
[YY0,YY]=meshgrid(yy,yy0);
[S0oL,SoL]=meshgrid(s0oL,soL);
%%
dX=XX-XX0;
dY=YY-YY0;
dSoL=SoL-S0oL;
[RR,varTH]=GEN_polar_coords(dX,dY);
%%
ip_stuff=IP_stuff{1};
ip_chi=ip_stuff{1}{2}(2:end,:);
w_Pn=IP_stuff{2};
ip_right=w_Pn.*ff;
%%
Gn=feval(GFxn,dX,dY,GF_args{:});
FF=ip_chi*Gn*ip_right;
%%
KEEPSING=GF_args{2};
%%
%% LIMIT OF sin(\Theta-\theta_\Delta)/r_\Delta
%% (COMMON TO INNER AND OUTER ROUTINES):



if ~KEEPSING
  hard_part=0*RR;
  [jz,rz]=find(RR==0);
  LIM=-.5*( d2xy_ds2(jz,2).*cos(th_vec(jz)) - ...
           d2xy_ds2(jz,1).*sin(th_vec(jz)) );
  hard_part(jz,rz)=LIM;
  %%
  jnz=find(RR);
  hard_part(jnz)=sin(TH(jnz)-varTH(jnz))./RR(jnz);
  FF=FF+ip_chi*(hard_part/2/pi)*ip_right;
end
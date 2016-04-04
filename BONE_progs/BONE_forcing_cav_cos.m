function FF=BONE_forcing_cav_cos(GEOM_stuff,...
              IP_stuff,GFxn,GF_args);

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
Xtra0=geom0_stuff{5};%={th_vec,dtheta_ds,...
%   ds_dt,d2s_dt2,d2theta_ds2,d2xy_ds2}
%%
d2xy_ds2=Xtra0{6};
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
w_qd=IP_stuff{2};
ip_right=w_qd.*ff;
%%
if 0
  Gn=feval(GFxn,dX,dY,GF_args{:});
else
  GF_args=[{dX,dY},GF_args];
  Gn=feval(GFxn,GF_args);
end
FF=ip_chi*Gn*ip_right;
%%
KEEPSING=GF_args{4};
%%
%% LIMIT OF sin(\Theta-\theta_\Delta)/r_\Delta
%% (COMMON TO INNER AND OUTER ROUTINES):

if ~KEEPSING
  hard_part=0*RR;
  [jz,rz]=find(RR==0);%[jz,rz]
  LIM_thingee=.5*( d2xy_ds2(:,2).*cos(th_vec0) - ...
           d2xy_ds2(:,1).*sin(th_vec0) );
  [LIM,TH]=meshgrid(LIM_thingee,th_vec);
  hard_part(jz,rz)=LIM(jz,rz);
  %%
  jnz=find(RR);
  hard_part(jnz)=sin(TH(jnz)-varTH(jnz))./RR(jnz);
  FF=FF+ip_chi*(hard_part/2/pi)*ip_right;
end

if 0%%plot Gn:
  for j=1:length(soL)
    y=Gn(j,:)+hard_part(j,:)/2/pi;
%      y=hard_part(j,:);
    plot(s0oL,y'), hold on;
%      plot(s0oL(j),LIM_thingee(j),'ok');
    plot(s0oL(j),y(j),'.r'), hold off;
    pause;
  end
end

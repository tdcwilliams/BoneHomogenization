function MK=BONE_kernelmat_dnlog(GEOM_stuff,IP_stuff,whichdiff);

if nargin==2
  whichdiff=1;
end

%% stuff for s integral - Fourier quadrature points:
geom_stuff=GEOM_stuff{1};
soL=geom_stuff{1};
xx=geom_stuff{2};
yy=geom_stuff{3};
LC=geom_stuff{4};
theta=geom_stuff{5};

%% stuff for s_0 integral - Fourier quadrature points:
geom0_stuff=GEOM_stuff{2};
s0oL=geom0_stuff{1};
xx0=geom0_stuff{2};
yy0=geom0_stuff{3};
LC0=geom0_stuff{4};
theta0=geom0_stuff{5};
d2xy_ds2=geom0_stuff{6};
%={th_vec,dtheta_ds,...
%   ds_dt,d2s_dt2,d2theta_ds2,d2xy_ds2}
%%
%%
[XX0,XX]=meshgrid(xx0,xx);
[YY0,YY]=meshgrid(yy0,yy);
[S0oL,SoL]=meshgrid(s0oL,soL);
%%
dX=XX-XX0;
dY=YY-YY0;
dSoL=SoL-S0oL;
[RR,varTH]=GEN_polar_coords(dX,dY);
%%
ip_left=IP_stuff{1};
ip_right=IP_stuff{2};%(2:end,:);
%%
Gn=0*XX;
[jz,rz]=find(RR==0);%[jz,rz]
%%
jnz=find(RR);
LIM_thingee=.5*( d2xy_ds2(:,2).*cos(theta0) + ...
                   - d2xy_ds2(:,1).*sin(theta0) );
%  Gn(jnz)=sin(TH(jnz)-varTH(jnz))./RR(jnz);
if whichdiff==1
  [LIM,TH]=meshgrid(LIM_thingee,theta);
  Gn(jz,rz)=LIM(jz,rz);
  %%
  Gn(jnz)=( sin(TH(jnz)).*dX(jnz)+...
           -cos(TH(jnz)).*dY(jnz) )./RR(jnz).^2;
else
  [TH0,LIM]=meshgrid(theta,LIM_thingee);
  Gn(jz,rz)=LIM(jz,rz);
  %%
  Gn(jnz)=( sin(TH0(jnz)).*dX(jnz)+...
            - cos(TH0(jnz)).*dY(jnz) )./RR(jnz).^2;
end
MK=ip_left*(Gn/2/pi)*ip_right';

if 0%nargin==3%%plot Gn:
  z01=xx0(1)+1i*yy0(1)
%    subplot(1,3,1)
  for j=10%1:length(soL)
    y=Gn(:,j)/2/pi;
%      y=hard_part(j,:);
    plot(soL,y), hold on;
%      plot(s0oL(j),LIM_thingee(j),'ok');
%      plot(s0oL(j),y(j),'.r'), hold off;
%      pause;
  end
end

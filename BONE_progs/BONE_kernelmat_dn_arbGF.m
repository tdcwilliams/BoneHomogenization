function MK=BONE_kernelmat_dn_arbGF(GEOM_stuff,...
              IP_stuff,GF,whichdiff);

if nargin<=3
  whichdiff=1;
end

if nargin==2
  MK=BONE_kernelmat_dnlog(GEOM_stuff,IP_stuff,whichdiff);
  return;
elseif isempty(GF)
  MK=BONE_kernelmat_dnlog(GEOM_stuff,IP_stuff,whichdiff);
  return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%GET GREEN'S FXN:
GFxn=GF{1};
GF_args=GF(2:end);
if whichdiff==1
  normal_deriv={1,theta};
else
  normal_deriv={-1,theta0};
end
USE_SING=0;
GF_args=[{dX,dY},GF_args,{normal_deriv,USE_SING}];
Gn0=feval(GFxn,GF_args);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ip_left=IP_stuff{1};
ip_right=IP_stuff{2};%(2:end,:);
%%
Gn=0*XX;
[jz,rz]=find(RR==0);%[jz,rz]
LIM_thingee=.5*( d2xy_ds2(:,2).*cos(theta0) + ...
                   - d2xy_ds2(:,1).*sin(theta0) );
%%
jnz=find(RR);
%  Gn(jnz)=sin(TH(jnz)-varTH(jnz))./RR(jnz);
if whichdiff==1%% calc \pa_n:
  [LIM,TH]=meshgrid(LIM_thingee,theta);
  Gn(jz,rz)=LIM(jz,rz);
  %%
  Gn(jnz)=( sin(TH(jnz)).*dX(jnz)+...
            - cos(TH(jnz)).*dY(jnz) )./RR(jnz).^2;
else%% calc \pa_n':
  [TH0,LIM]=meshgrid(theta,LIM_thingee);
  Gn(jz,rz)=LIM(jz,rz);
  %%
  Gn(jnz)=( sin(TH0(jnz)).*dX(jnz)+...
            - cos(TH0(jnz)).*dY(jnz) )./RR(jnz).^2;
end
MK=ip_left*(Gn0+Gn/2/pi)*ip_right';

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

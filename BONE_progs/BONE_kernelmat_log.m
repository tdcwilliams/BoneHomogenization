function MK=BONE_kernelmat_log(...
              GEOM_stuff,IP_stuff,CLOSED_CURVE)

geom_stuff=GEOM_stuff{1};
soL=geom_stuff{1};
xx=geom_stuff{2};
yy=geom_stuff{3};
LC=geom_stuff{4};
theta=geom_stuff{5};
%%
geom0_stuff=GEOM_stuff{2};
s0oL=geom0_stuff{1};
xx0=geom0_stuff{2};
yy0=geom0_stuff{3};
LC0=geom0_stuff{4};
theta0=geom0_stuff{5};
%%
[XX0,XX]=meshgrid(xx0,xx);
[YY0,YY]=meshgrid(yy0,yy);
[S0oL,SoL]=meshgrid(s0oL,soL);
%%
dX=XX-XX0;
dY=YY-YY0;
dSoL=SoL-S0oL;
RR=abs(dX+1i*dY);
%%
ip_left=IP_stuff{1};
ip_right=IP_stuff{2};
Mlog=IP_stuff{3};
%%
if isempty(Mlog)
  GG=log(RR)/2/pi;
  MK=ip_left*GG*ip_right';
elseif CLOSED_CURVE
  log_RoET=0*RR+log(LC/pi);
  jnz=find(RR);
  %%
  Exp_term=abs(1-exp(1i*pi*dSoL));
  log_RoET(jnz)=log(RR(jnz)./Exp_term(jnz));
  MK=Mlog+ip_left*(log_RoET/2/pi)*ip_right';
else
  log_RoDS=0*RR;
  jnz=find(RR);
  %%
  DS=LC*abs(dSoL);
  log_RoDS(jnz)=log(RR(jnz)./DS(jnz));
  MK=Mlog+ip_left*(log_RoDS/2/pi)*ip_right';
end

if 0%%plot G:
  disp('plotting G:');
  for j=1:length(soL)
    y=GG(j,:)+log_RoET(j,:)/2/pi;
%      y=hard_part(j,:);
    plot(s0oL,y'), hold on;
%      plot(s0oL(j),LIM_thingee(j),'ok');
    plot(s0oL(j),y(j),'.r'), hold off;
    pause;
  end
end
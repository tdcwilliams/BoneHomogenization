function MK=BONE_kernelmat_G(...
              GEOM_stuff,IP_stuff,GFxn,GF_args)

geom_stuff=GEOM_stuff{1};
soL=geom_stuff{1};
xx=geom_stuff{2};
yy=geom_stuff{3};
LC=geom_stuff{4};
Xtra=geom_stuff{5};
%%
geom0_stuff=GEOM_stuff{2};
s0oL=geom0_stuff{1};
xx0=geom0_stuff{2};
yy0=geom0_stuff{3};
LC0=geom0_stuff{4};
Xtra0=geom0_stuff{5};
%%
[XX0,XX]=meshgrid(xx,xx0);
[YY0,YY]=meshgrid(yy,yy0);
[S0oL,SoL]=meshgrid(s0oL,soL);
%%
dX=XX-XX0;
dY=YY-YY0;
dSoL=SoL-S0oL;
RR=abs(dX+1i*dY);
%%
ip_d1chi=IP_stuff{1}{1};%(2:end,:);
GG=feval(GFxn,dX,dY,GF_args{:});
MK=ip_d1chi*GG*ip_d1chi';
%%
KEEPSING=GF_args{4};
if ~KEEPSING
  Mlog=IP_stuff{end}(2:end,2:end);
  log_RoET=0*RR+log(LC/pi);
  jnz=find(RR);
  %%
  Exp_term=abs(1-exp(1i*pi*dSoL));
  log_RoET(jnz)=log(RR(jnz)./Exp_term(jnz));
  MK=MK+Mlog+ip_d1chi*(log_RoET/2/pi)*ip_d1chi';
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
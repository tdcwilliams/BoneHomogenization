function MK=BONE_kernelmat_crk_Tn(...
              GEOM_stuff,ip_stuff,GFxn,GF_args)

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
ipT=ip_stuff{1}{1}(2:end,:);
GG=feval(GFxn,dX,dY,GF_args{:});
MK=ipT*GG*ipT';
%%
KEEPSING=GF_args{2};
if ~KEEPSING
  Mlog=ip_stuff{end}(2:end,2:end);
  log_RoDSoL=0*RR+log(LC);
  jnz=find(RR);
  log_RoDSoL(jnz)=log(RR(jnz)./abs(dSoL(jnz)));
  MK=MK+Mlog+ipT*(log_RoDSoL/2/pi)*ipT';
end
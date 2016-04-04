function MK=BONE_kernelmat_arbGF(...
              GEOM_stuff,IP_stuff,CLOSED_CURVE,GF)

if nargin==3
  MK  = BONE_kernelmat_log(GEOM_stuff,IP_stuff,CLOSED_CURVE);
  return;
elseif isempty(GF)
  MK  = BONE_kernelmat_log(GEOM_stuff,IP_stuff,CLOSED_CURVE);
  return;
end
%  CLOSED_CURVE,pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET GEOMETRIC STUFF FOR (x,y) CURVE:
geom_stuff  = GEOM_stuff{1};
soL         = geom_stuff{1};
xx          = geom_stuff{2};
yy          = geom_stuff{3};
LC          = geom_stuff{4};
theta       = geom_stuff{5};
if 0%%plot curve
   %figure;
   plot(xx,yy);
   hold on;
%  plot(0*yy,yy,'k');
%  plot(xx,0*yy,'k');
   %hold off;
   axis([-1 1 -1 1]);
   %pause;
end

%% GET GEOMETRIC STUFF FOR (x',y') CURVE:
geom0_stuff = GEOM_stuff{2};
s0oL        = geom0_stuff{1};
xx0         = geom0_stuff{2};
yy0         = geom0_stuff{3};
LC0         = geom0_stuff{4};
theta0      = geom0_stuff{5};
%%
[XX0,XX]    = meshgrid(xx0,xx);
[YY0,YY]    = meshgrid(yy0,yy);
[S0oL,SoL]  = meshgrid(s0oL,soL);
%%
%  [XX0;YY0],[xx,yy],
dX    = XX-XX0;
dY    = YY-YY0;
dSoL  = SoL-S0oL;
RR    = abs(dX+1i*dY);

if 0%%plot curve
   %figure;
   plot(xx0,yy0);
   hold on;
%  plot(0*yy,yy,'k');
%  plot(xx,0*yy,'k');
   %hold off;
   axis([-1 1 -1 1]);
   disp('paused'),pause;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET NON-SINGULAR BIT OF GREEN'S FXN:
GFxn           = GF{1};
GF_args        = GF(2:end);
normal_deriv   = [];
USE_SING       = 0;
GF_args        = [{dX,dY},GF_args,{normal_deriv,USE_SING}];
GG0            = feval(GFxn,GF_args);%GG0(1:10,1:10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ip_left  = IP_stuff{1};
ip_right = IP_stuff{2};
Mlog     = IP_stuff{3};
%%
if isempty(Mlog)
  GG  = GG0+log(RR)/2/pi;
  MK  = ip_left*GG*ip_right';
elseif CLOSED_CURVE
  log_RoET  = 0*RR+log(LC/pi);
  jnz       = find(RR);
  %%
  Exp_term        = abs(1-exp(1i*pi*dSoL));
  log_RoET(jnz)   = log(RR(jnz)./Exp_term(jnz));
  GG              = GG0+log_RoET/2/pi;
  MK              = Mlog+ip_left*GG*ip_right';
else
  log_RoDS  = 0*RR;
  jnz       = find(RR);
  %%
  DS              = LC*abs(dSoL);
  log_RoDS(jnz)   = log(RR(jnz)./DS(jnz));
  GG              = GG0+log_RoDS/2/pi;
  MK              = Mlog+ip_left*GG*ip_right';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0%%plot G:
  disp('plotting G:');
  for j=1:length(soL)
    y = GG(j,:)+log_RoET(j,:)/2/pi;
%      y=hard_part(j,:);
    plot(s0oL,y'), hold on;
%      plot(s0oL(j),LIM_thingee(j),'ok');
    plot(s0oL(j),y(j),'.r'), hold off;
    pause;
  end
end

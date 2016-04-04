%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alp_c,fvalue,Mhat] = BONE_alpcrit_strt_rot(aoA,rot,guess)

do_test  = 0;
if nargin==0
   aoA      = .3;
   rot      = 90;
   do_test  = 1;
end

if nargin==3
   [alp_c,fvalue,exit_flag]   = fminsearch(@(alp)fn_alp_crit_strt(alp,aoA,rot),guess);
   if do_test
      int                     = alp_c*[.95 1.05];
      [fmin_value90,Mhat90]   = fn_alp_crit_strt(1,aoA,90);
      mx90                    = Mhat90(1,1);
      alp_c90                 = 1/sqrt(mx90)
   end
else
   %%put rot=90 & use that as a guess;
   [fmin_value90,Mhat90] = fn_alp_crit_strt(1,aoA,90);
   %%
   mx90     = Mhat90(1,1);
   alp_c90  = 1/sqrt(mx90);
   alp_c0   = 1;
%  if 0
%     guess = alp_c0+rot/90*(alp_c90-alp_c0);
%     alp_c = fminsearch(@(alp)fn_alp_crit_strt(alp,aoA,rot),guess);
%  else
      int      = [alp_c0,alp_c90];
      options  = optimset('TolX',1e-8);
      [alp_c,fvalue,exitflag] = fminbnd(@(alp)fn_alp_crit_strt(alp,aoA,rot),int(1),int(2),options);
%  end
end

if nargout==3
   [fmin_value,Mhat] = fn_alp_crit_strt(alp_c,aoA,rot);
   tst               = eig(Mhat);
end

if do_test
   [fmin_value,Mhat] = fn_alp_crit_strt(alp_c,aoA,rot)
   tst               = eig(Mhat)
   %%
   if 1
      [fmin_value90,Mhat90] = fn_alp_crit_strt(alp_c90,aoA,90)
   end
   %%
   if 1
      di = int(2)-int(1);
      xx = linspace(int(1),int(2)+.2*di,25)';
      yy = fn_alp_crit_strt(xx,aoA,rot);
      pause(0.1);
      plot(xx,yy,'k',alp_c,fvalue,'or',alp_c90,fmin_value90,'xg');
   end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fmin,Mhat] = fn_alp_crit_strt(alp_vec,aoA,rot)

crk_solver  = @BONE_GF_cracks_rect_cell;
%%
A     = .5;
Bhat  = .5;
a     = aoA*A;
%%
Npolys   = 50;
Na       = length(alp_vec);
for j = 1:Na
   alp      = alp_vec(j);
   Dalp     = diag([1 1/alp]);
   B        = Bhat*alp;
   AB       = [A,B];
   %%
   z_rot    = exp(1i*pi*rot/180);
   Irr_vars = CURVEget_strtline(-z_rot*a,z_rot*a);
   %%
   C_eff    = feval(crk_solver,AB,Irr_vars,Npolys);
   Mhat     = Dalp*C_eff{1}*Dalp;
   %%
   %fzero = Mhat(1,1)-Mhat(2,2);
   m1       = Mhat(1,1);
   m2       = Mhat(2,2);
   m12      = Mhat(1,2);
   fmin(j)  = (m1-m2)^2+4*m12^2;
   if Na>10
      disp([j,Na])
   end
end

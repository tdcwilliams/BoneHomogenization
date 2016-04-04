%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alp_c,fvalue,Mhat] = BONE_alpcrit_cos_phase(aoA,boB,phi,guess)

%[fzero,Mhat]   = fn_alp_crit_cos(1,aoA,boB);
%alp0           = sqrt(Mhat(2,2));
if nargin==0
   aoA   = .5;
   boB   = .5;
   phi   = 0.05;
end


if ~exist('guess')
   guess = BONE_alpcrit_cos(aoA,boB);
end

if phi==0
   alp_c          = guess;
   [fvalue,Mhat]  = fn_alp_crit_cos(alp_c,aoA,boB,phi);
   return;
end

[alp_c,fvalue,exit_flag]   = fminsearch(@(alp)fn_alp_crit_cos(alp,aoA,boB,phi),guess);
if nargout==3
   [fvalue,Mhat]  = fn_alp_crit_cos(alp_c,aoA,boB,phi);
end

if 0
   int                     = alp_c*[.95 1.05];
   [fmin_value90,Mhat90]   = fn_alp_crit_cos(1,aoA,boB,0);
   mx90                    = Mhat90(1,1);
   alp_c90                 = 1/sqrt(mx90)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fzero,Mhat] = fn_alp_crit_cos(alp,aoA,boB,phi)

crk_solver  = @BONE_GF_cracks_rect_cell;
crk_fxn     = @CURVEprof_cos;
%%
A     = .5;
Bhat  = .5;
a     = aoA*A;
%%
Npolys   = 100;
Dalp     = diag([1 1/alp]);
B        = Bhat*alp;
b        = boB*B;
AB       = [A,B];
%%
phase    = phi;
crk_vars = {1,phase};
srt      = {[a,b],0,[0 0]};
Irr_vars = {crk_fxn,crk_vars,srt};
C_eff    = feval(crk_solver,AB,Irr_vars,Npolys);
Mhat        = Dalp*C_eff{1}*Dalp;
%%
fzero = (Mhat(1,1)-Mhat(2,2))^2+4*Mhat(1,2)^2;

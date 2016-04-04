%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function boB_c = BONE_alpcrit_cos_VboB(aoA,alp)

if nargin==0
   alp   = 1;
   aoA   = .5;
end

disp('finding alp_c for b=0');
alp_c0   = BONE_alpcrit_strt(aoA)
if alp<alp_c0
   boB_c = NaN;
   return;
elseif alp==alp_c0
   boB_c = 0;
   return;
end

disp('finding alp_c for b=B');
alp_c1   = BONE_alpcrit_cos(aoA,1)
if alp>alp_c1
   boB_c = NaN;
   return;
elseif alp==alp_c1
   boB_c = 1;
   return;
end

db    = 1e-5;
int   = [db,1-db];

%[fzero,Mhat]   = fn_alp_crit_cos(boB,aoA,alp);
%alp0           = sqrt(Mhat(2,2));


disp('finding b_c for given alp');
boB_c =...
   GEN_findroots_bisection(@fn_alp_crit_cos,int,aoA,alp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fzero,Mhat] = fn_alp_crit_cos(boB,aoA,alp)

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
phase    = 0;
crk_vars = {1,phase};
srt      = {[a,b],0,[0 0]};
Irr_vars = {crk_fxn,crk_vars,srt};
C_eff    = feval(crk_solver,AB,Irr_vars,Npolys);
Mhat        = Dalp*C_eff{1}*Dalp;
%%
fzero = Mhat(1,1)-Mhat(2,2);

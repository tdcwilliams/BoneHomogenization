%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alp_c = BONE_get_alp_crit_strt(aoA)

[fzero,Mhat]   = fn_alp_crit_strt(1,aoA);
alp0           = sqrt(Mhat(2,2));

%% find lower limit
al = .95*alp0;
ah = max(1,1.001*alp0);
fl = fn_alp_crit_strt(al,aoA);
its   = 1;
while fl>0
   its   = its+1;
   ah    = al;
   al    = al-.05*alp0;
   fl    = fn_alp_crit_strt(al,aoA);
end

ints  = [al,ah];
alp_c = GEN_findroots_bisection(@fn_alp_crit_strt,ints,aoA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fzero,Mhat] = fn_alp_crit_strt(alp,aoA)

crk_solver  = @BONE_GF_cracks_rect_cell;
%%
A     = .5;
Bhat  = .5;
a     = aoA*A;
%%
Npolys   = 50;
Dalp     = diag([1 1/alp]);
B        = Bhat*alp;
AB       = [A,B];
%%
Irr_vars = CURVEget_strtline(-a,a);
C_eff    = feval(crk_solver,AB,Irr_vars,Npolys);
Mhat     = Dalp*C_eff{1}*Dalp;
%%
fzero = Mhat(1,1)-Mhat(2,2);

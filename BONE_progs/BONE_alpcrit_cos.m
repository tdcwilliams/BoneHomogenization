%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alp_c = BONE_alpcrit_cos(aoA,boB)

[fzero,Mhat]   = fn_alp_crit_cos(1,aoA,boB);
alp0           = sqrt(Mhat(2,2));

if fzero==0%% find lower limit
   alp_c = 1;
   return;
elseif fzero>0
   al    = .95*alp0;
   ah    = 1.001*alp0;
   fl    = fn_alp_crit_cos(al,aoA,boB);
   fh    = fn_alp_crit_cos(ah,aoA,boB);
   if fh==0
      alp_c = ah;
      return;
   elseif fh<0
      ints  = [ah,1];
   else
      its   = 1;
      while fl>0
         its   = its+1;
         ah = al;
         al = al-.05*alp0;
         fl = fn_alp_crit_cos(al,aoA,boB);
      end
      ints  = [al,ah];
   end
else
   al    = .999*alp0;
   ah    = 1.05*alp0;
   fh    = fn_alp_crit_cos(ah,aoA,boB);
   fl    = fn_alp_crit_cos(al,aoA,boB);
   its   = 1;
   if fl==0
      alp_c = al;
      return;
   elseif fl>0
      ints  = [1,al];
   else
      while fh<0
         its   = its+1;
         al    = ah;
         ah    = ah+.05*alp0;
         fh = fn_alp_crit_cos(ah,aoA,boB);
      end
      ints  = [al,ah];
   end
end

alp_c =...
   GEN_findroots_bisection(@fn_alp_crit_cos,ints,aoA,boB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fzero,Mhat] = fn_alp_crit_cos(alp,aoA,boB)

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

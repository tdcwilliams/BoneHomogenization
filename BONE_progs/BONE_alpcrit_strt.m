%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alp_c,Mhat] = BONE_alpcrit_strt(aoA,is_vert)

do_test  = 0;
if nargin==0
   aoA      = .5;
   is_vert  = 1;
   do_test  = 1;
end

if ~exist('is_vert')
   is_vert  = 0;
end

[fzero,Mhat]   = fn_alp_crit_strt(1,aoA,is_vert);
if is_vert==0
   alp0           = sqrt(Mhat(2,2));
   
   %% find lower limit
   al = .95*alp0;
   ah = max(1,1.001*alp0);
   fl = fn_alp_crit_strt(al,aoA,is_vert);
   its   = 1;
   while fl>0
      its   = its+1;
      ah    = al;
      al    = al-.05*alp0;
      fl    = fn_alp_crit_strt(al,aoA,is_vert);
   end
   
   ints  = [al,ah];
else
   alp0  = 1/sqrt(Mhat(1,1));
   ints  = [.95 1.05]*alp0;
end

alp_c = GEN_findroots_bisection(@fn_alp_crit_strt,ints,aoA,is_vert);

if nargout==2
   [fzero_val,Mhat]  = fn_alp_crit_strt(alp_c,aoA,is_vert);
end

if do_test
   [fzero_val,Mhat]  = fn_alp_crit_strt(alp_c,aoA,is_vert)
   if 1
      di = ints(2)-ints(1);
      xx = linspace(ints(1),ints(2)+.2*di,100)';
      yy = fn_alp_crit_strt(xx,aoA,is_vert);
      plot(xx,yy,'k',xx,0*yy,'b',alp_c,fzero_val,'or');
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fzero,Mhat] = fn_alp_crit_strt(alp_vec,aoA,is_vert)

if ~exist('is_vert')
   is_vert  = 0;
end
Na    = length(alp_vec);
fzero = zeros(Na,1);

crk_solver  = @BONE_GF_cracks_rect_cell;
%%
A     = .5;
Bhat  = .5;
a     = aoA*A;
%%
Npolys   = 50;

for j = 1:Na
   alp      = alp_vec(j);
   %%
   Dalp     = diag([1 1/alp]);
   B        = Bhat*alp;
   AB       = [A,B];
   %%
   if is_vert==0
      Irr_vars = CURVEget_strtline(-a,a);
   else
      Irr_vars = CURVEget_strtline(-a*1i,a*1i);
   end
   C_eff    = feval(crk_solver,AB,Irr_vars,Npolys);
   Mhat     = Dalp*C_eff{1}*Dalp;
   %%
   fzero(j) = Mhat(1,1)-Mhat(2,2);
   if Na>10
      disp([j Na]);
   end
end

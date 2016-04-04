function fig_alp_crit(j0)

na       = 40;
aoA_vec  = linspace(.05,.95,na)';
DO_MORE  = 0;
if DO_MORE
   aoA_vec_ = [aoA_vec;(.955:.005:.995)'];
   na       = length(aoA_vec_);
end
if 1
   aoA_vec  = [aoA_vec;(.955:.005:.995)'];
   na       = length(aoA_vec);
end
%%
if 0
   boB_vec  = [.1 .2 .3];
   filename = 'out/alp_crit_strt.mat'
elseif 0
   boB_vec  = [.2 .4 .6 .8];
   filename = 'out/alp_crit_all.mat'
elseif 0
   boB_vec  = .99;
   filename = 'out/alp_crit_lims.mat'
   figname  = 'out/alp_crit_lims.eps'
elseif 0
   boB_vec  = 1;
   filename = 'out/alp_crit_lims2.mat'
elseif 0
   boB_vec  = [.85 .9 .95]
   filename = 'out/alp_crit_fill.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else%%combine all
   tmp1  = [];%%boB
   tmp2  = [];%%alp_crit
   if 0
      fname = {'out/alp_crit_strt.mat',...
               'out/alp_crit_all.mat',...
               'out/alp_crit_lims.mat',...
               'out/alp_crit_lims2.mat',...
               'out/alp_crit_combo.mat'};
      for r=1:4
         load(fname{r});
         tmp1  = [tmp1,boB_vec];
         tmp2  = [tmp2,alp_crit(:,2:end)];
      end
      [boB_vec,JJ]   = sort(tmp1);
      alp_crit       = [alp_crit(:,1),tmp2(:,JJ)];
      %%
      save(fname{5},'boB_vec','aoA_vec','alp_crit');
      filename = fname{5}
      boB_vec
   elseif 0
      fname = {'out/alp_crit_combo.mat',...
               'out/alp_crit_fill.mat',...
               'out/alp_crit_combo2.mat'};
      for r=1:2
         load(fname{r});
         tmp1  = [tmp1,boB_vec];
         tmp2  = [tmp2,alp_crit(:,2:end)];
      end
      [boB_vec,JJ]   = sort(tmp1);
      alp_crit       = [alp_crit(:,1),tmp2(:,JJ)];
      %%
      save(fname{3},'boB_vec','aoA_vec','alp_crit');
      filename = fname{3}
      boB_vec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   else
      if 0
         filename = 'out/alp_crit_combo.mat'
         boB_vec  = [0.1000,0.2000,0.3000,...
                     0.4000,0.6000,0.8000,.99,1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      else
         filename = 'out/alp_crit_combo2.mat'
         boB_vec  = [0.1000,0.2000,0.3000,...
                     0.4000,0.6000,0.8000,...
                     .85,.9,.95,.99,1]
      end
   end
   %figname  = 'out/alp_crit_combo.eps'
   figname  = 'out/alp_crit_combo2.eps'
   %return
end
nb       = length(boB_vec);
%%
alp_crit = zeros(na,nb+1);

if nargin==0
   j0       = 1;
   DO_RUN   = (~exist(filename));
else
   DO_RUN   = 1;
   if exist(filename)
      load(filename);
      if DO_MORE
         aoA_vec  = aoA_vec_;
      end
      %na = 10
   end
end

if DO_RUN
   %smh   = 0;%% time in seconds;
   smh   = 1;%% time in mins;
   t0    = GEN_time(smh);
   for j=j0:na
      aoA            = aoA_vec(j);
      alp_crit(j,1)  = get_alp_crit_strt(aoA);
      for r=1:nb
         disp([r nb]);
         boB               = boB_vec(r);
         alp_crit(j,r+1)   = get_alp_crit_cos(aoA,boB);
      end
      GEN_progrep([j na],t0,smh);
      save(filename,'aoA_vec','alp_crit','boB_vec');
   end
else
   load(filename);
end
%%
GEN_set_plotorder;
%alp_crit = fliplr(alp_crit);
if 0
   plot(aoA_vec,alp_crit);
else%%alp_crit_combo2.mat;
   X     = [.075 .99];
   Y     = [0.25 3.5];
   ax    = .7;
   ay    = .75;
   x     = X(1)+ax*(X(2)-X(1));
   y     = Y(1)+ay*(Y(2)-Y(1));
   %%
   JJ = [1,7:8,10:13];
   plot(aoA_vec,alp_crit(:,JJ(1:4)));
   hold on;
   plot(aoA_vec,alp_crit(:,JJ(5:end)),'r');
   txt   = text(x,y,'(a)');
   GEN_font(txt,17);
   boB_vec(JJ(2:end)-1)
   %%
   hold off;
end
GEN_proc_fig('a/A','ac');
xlim(X);
ylim(Y);

if exist('figname')
   %saveas(gcf,figname);
   saveas(gcf,figname,'epsc');
end
%aoA   = .5;
%boB   = .3;
%alp_c = get_alp_crit_cos(aoA,boB)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alp_c = get_alp_crit_strt(aoA)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function alp_c = get_alp_crit_cos(aoA,boB)

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

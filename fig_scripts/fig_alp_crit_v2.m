function fig_alp_crit_v2(j0)

%%
filename = 'out/alp_crit_combo2.mat'
boB_vec  = [0.1000,.2,0.2000,0.3000,...
            0.4000,0.6000,0.8000,...
            .85,.9,.95,.99,1];
%figname  = 'out/alp_crit_combo.eps'
JJ       = [1,7:8,10:13];
%%
load(filename);
mx_crit  = 0*alp_crit;
na       = length(aoA_vec);
alp_crit = alp_crit(:,JJ);
%boB_vec
boB_vec(JJ(2:end)-1)
boB_vec([2 4:6]-1)
nb       = length(boB_vec);
%%
figname  = 'out/alp_crit.eps'
fname    = 'out/alp_crit_final.mat'
   %return
if 1
   subplot(1,2,1);
   GEN_set_plotorder;
   X     = [.075 .99];
   Y     = [0.25 3.5];
   ax    = .8;
   ay    = .925;
   x     = X(1)+ax*(X(2)-X(1));
   y     = Y(1)+ay*(Y(2)-Y(1));
   %%
   JJ = [1,7:8,10:13];
   plot(aoA_vec,alp_crit(:,1:4));
   hold on;
   plot(aoA_vec,alp_crit(:,5:end),'r');
   txt   = text(x,y,'(a)');
   GEN_font(txt,17);
   %%
   GEN_proc_fig('a/A','ac');
   xlim(X);
   ylim(Y);
   hold off;
end

if ~exist(fname)
   %smh      = 0;%% time in seconds;
   smh     = 1;%% time in mins;
   t0       = GEN_time(smh);
   for j=1:na
      aoA   = aoA_vec(j);
      alp   = alp_crit(j,1);
      %%
      [fzero,Mhat]   = fn_alp_crit_strt(alp,aoA);
      mx_crit(j,1)   = mean(diag(Mhat));
      for r=1:nb
         boB            = boB_vec(r);
         [fzero,Mhat]   = fn_alp_crit_cos(alp,aoA,boB);
         mx_crit(j,r+1) = mean(diag(Mhat));
      end
      GEN_progrep([j na],t0,smh);
      save(fname,'mx_crit','aoA_vec','alp_crit','boB_vec');
   end
else
   load(fname);
end

subplot(1,2,2);
GEN_set_plotorder;
Y     = [0.56 1.03];
ax    = .8;%755;
ay    = .075;
x     = X(1)+ax*(X(2)-X(1));
y     = Y(1)+ay*(Y(2)-Y(1));
   %%
if 1
   plot(aoA_vec,mx_crit(:,JJ(1:4)));
   hold on;
   plot(aoA_vec,mx_crit(:,JJ(5:end)),'r');
   jp    = [2,4:6];
   ls    = {'-or','-vr','-xr','-^r'};
   Np    = length(aoA_vec);
   jp2   = find(aoA_vec<.90);
   jp2   = (1:2:max(jp2))';
   jp2   = [jp2;( max(jp2)+2:3:Np )'];
   %jp2   = [(1:2:Np )'];
   for j=1:4
      if j==3
         plot(aoA_vec(jp2),mx_crit(jp2,jp(j)),ls{j},'markersize',9.5);%mx_crit(:,2:6),boB_vec(2:6)
      else
         plot(aoA_vec(jp2),mx_crit(jp2,jp(j)),ls{j});%mx_crit(:,2:6),boB_vec(2:6)
      end
   end
else
   plot(aoA_vec,mx_crit);
   hold on;
end
txt   = text(x,y,'(b)');
GEN_font(txt,17);
   %%
GEN_proc_fig('a/A','mx');
xlim(X);
ylim(Y);
hold off;
%%
saveas(gcf,figname,'epsc');
!gv out/alp_crit.eps &
return;

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
   ah = al;
   al = al-.05*alp0;
   fl = fn_alp_crit_strt(al,aoA);
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
Mhat     = Dalp*C_eff{1}*Dalp;
%%
fzero = Mhat(1,1)-Mhat(2,2);

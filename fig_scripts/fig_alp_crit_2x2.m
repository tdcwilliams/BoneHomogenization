%% fig_alp_crit_2x2.m

%%
nrows = 1;
figname  = 'out/alp_crit_2x2.eps'

if nrows==2
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
   fname    = 'out/alp_crit_final.mat'
      %return

   if 1
      subplot(nrows,2,1);
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

   subplot(nrows,2,2);
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
   subplot(nrows,2,2*nrows-1);
   fig_alp_crit_rot_strt;
   %%
   subplot(nrows,2,2*nrows);
   fig_alp_crit_cos_phase;
else
   mar   = .15;
   %gap   = .18;
   gap   = .13;
   wid   = (1-2*mar-gap)/2;
   hei   = 1-2*mar;
   xpos  = [mar mar wid hei;mar+wid+gap mar wid hei];

   subplot('position',xpos(1,:));
   fig_alp_crit_rot_strt;
   set(xxl1,'position',xlp1-[0.245 0 0]);
   %%
   subplot('position',xpos(2,:));
   fig_alp_crit_cos_phase;
   set(xxl2,'position',xlp2-[0 .7e-4 0]);
   set(yyl2,'position',ylp2+[.75e-1 0 0]);
end
%%
saveas(gcf,figname,'epsc');
!gv out/alp_crit_2x2.eps &
return;

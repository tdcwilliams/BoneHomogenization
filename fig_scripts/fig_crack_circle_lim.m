%% fig_crack_circle_lim.m
%% plots (a) m_x and (b) m_y for elliptical arc cracks
%%  as the fraction of the ellipse covered tends to 1;
%% we check that m_x & m_y tend to those of an elliptical cavity
%%  (m_xy=0 here as axes are aligned with sides of unit cell);
%% eccentricity of ellipse ('ecc') is 0.25;

filename = 'out/crack_circle_lim.mat'
if exist(filename)
   load(filename)
else

   A  = .5;
   B  = .5;
   AB = [A,B];
   np = 50;
   %%
   aoA_vec  = (1:np)'/np*.5;
      %% a = major axis of ellipse (parallel to x axis);
   ecc      = 0.25;
      %% b/a = eccentricity of ellipse;

   %%
   Npolys      = 10;
   fraxn_vec   = [.9 .7 .5];
   Nterms      = 15;
   %  m_y=0*Ach_Li;
   for j=1:np
     a   = aoA_vec(j)*A;[j np]
     %%
     crk_fxn   = @CURVEprof_circarc;
     crk_prams = {1};%% fraction of circle
     radius    = a;
     srt       = {radius*[1 ecc],0,[0 0]};
     Irr_vars  = {crk_fxn,crk_prams,srt};
     %%
     C_eff     = BONE_GF_cavs_rect_cell(AB,Irr_vars,Nterms);
     %m_y(j,1)  = C_eff{1}(4);%[m_y(j),Ach_Li(j)]
     %m_x(j,1)  = C_eff{1}(1);
     %m_xy(j,1) = C_eff{1}(2);
     Mhat      = C_eff{1};
     %%
     [U,Lam]   = eig(Mhat);
     r         = 0;
     m1(j,r+1) = min(diag(Lam));
     m2(j,r+1) = max(diag(Lam));
     %%
     jj          = find(diag(Lam)==m1(j,r+1));
     pdir(j,r+1) = angle(U(1,jj)+1i*U(2,jj));
     %%
     if j>1
        if pdir(j,r+1)-pdir(j-1,r+1)>.5*pi
           pdir(j,r+1)   = pdir(j,r+1)-pi;
        elseif pdir(j,r+1)-pdir(j-1,r+1)<-.5*pi
           pdir(j,r+1)   = pdir(j,r+1)+pi;
        end
     end
     %%
     for r=1:length(fraxn_vec)
       crk_fxn    = @CURVEprof_circarc;
       crk_prams  = {fraxn_vec(r)};%% fraction of circle
       radius     = a;
       srt        = {radius*[1 ecc],0,[0 0]};
       Irr_vars   = {crk_fxn,crk_prams,srt};
       %%
       C_eff         = BONE_GF_cracks_rect_cell(AB,Irr_vars,Npolys);
       %pause
       %m_y(j,r+1)    = C_eff{1}(4);
       %m_x(j,r+1)    = C_eff{1}(1);
       %m_xy(j,r+1)   = C_eff{1}(2);
       Mhat = C_eff{1};
       %%
       [U,Lam]    = eig(Mhat);
       m1(j,r+1)  = min(diag(Lam));
       m2(j,r+1)  = max(diag(Lam));
       %%
       jj            = find(diag(Lam)==m1(j,r+1));
       pdir(j,r+1)   = angle(U(1,jj)+1i*U(2,jj));
       %%
       if j>1
          if pdir(j,r+1)-pdir(j-1,r+1)>.5*pi
             pdir(j,r+1)   = pdir(j,r+1)-pi;
          elseif pdir(j,r+1)-pdir(j-1,r+1)<-.5*pi
             pdir(j,r+1)   = pdir(j,r+1)+pi;
          end
       end
     end
   end

   %save(filename,'aoA_vec','m_x','m_y','m_xy');
   save(filename,'aoA_vec','m1','m2','pdir');
end

if 0
   subplot(1,3,1);
   GEN_set_plotorder;
   %plot(aoA_vec,m_y,'k'), hold on;
   plot(aoA_vec,m1,'k'), hold on;
   %%
   subplot(1,3,2);
   GEN_set_plotorder;
   %plot(aoA_vec,m_x,'k')
   plot(aoA_vec,m2,'k')
   %%
   subplot(1,3,3);
   GEN_set_plotorder;
   %plot(aoA_vec,m_x,'k')
   plot(aoA_vec,pdir,'k')
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   subplot(1,3,1);
   xlim([0 .5]);
   box on;
   %GEN_proc_fig('\ita/\itA','\itm^*_{\ity}');
   %GEN_proc_fig('a/A','my');
   GEN_proc_fig('a/A','m1');
   %%
   subplot(1,3,2);
   xlim([0 .5]);
   box on;
   %GEN_proc_fig('\ita/\itA','\itm^*_{\itx}');
   %GEN_proc_fig('a/A','mx');
   GEN_proc_fig('a/A','m2');
   %%
   subplot(1,3,3);
   xlim([0 .5]);
   box on;
   %GEN_proc_fig('\ita/\itA','\itm^*_{\itx}');
   %GEN_proc_fig('a/A','mx');
   GEN_proc_fig('a/A','chi/pi');
else
   subplot(1,2,1);
   GEN_set_plotorder;
   %plot(aoA_vec,m_y,'k'), hold on;
   plot(aoA_vec,m1,'k'), hold on;
   %%
   subplot(1,2,2);
   GEN_set_plotorder;
   %plot(aoA_vec,m_x,'k')
   plot(aoA_vec,m2,'k')

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   subplot(1,2,1);
   xlim([0 .5]);
   box on;
   %GEN_proc_fig('\ita/\itA','\itm^*_{\ity}');
   GEN_proc_fig('a/A','my');
   %GEN_proc_fig('a/A','m1');
   %%
   subplot(1,2,2);
   xlim([0 .5]);
   box on;
   %GEN_proc_fig('\ita/\itA','\itm^*_{\itx}');
   GEN_proc_fig('a/A','mx');
   %GEN_proc_fig('a/A','m2');
end

%saveas(gcf,'out/crack_circle_lim.fig')
%saveas(gcf,'out/crack_circle_lim.jpg')
saveas(gcf,'out/crack_circle_lim.eps')

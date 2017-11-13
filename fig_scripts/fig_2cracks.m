%% fig_crack_2ringres.m
%% plots (a) m_x and (b) m_y for 2 concentric
%%  circular arc cracks. The openings point in different
%%  directions - like in the 2 ring resonator problem. Because
%%  there is a resonance for low freq's our results become invalid
%%  as the rings close up, but including this problem shows the
%%  numerical method can handle it;
%% parameters of the problem looked at are
%%  - fraction of circles closed up;
%%  - ratio of 2 radii;
%% ?? expect m_y to be less reduced ??

clear;
outdir   = 'out'
if ~exist(outdir,'dir')
   mkdir(outdir);
end
filename = [outdir,'/crack_2cracks'];

DO2      = 1;

if exist([filename,'.mat'])
   load(filename)
   if 0
      fraxn_vec   = [.9 .8 .7 .6];
      %b_vec       = [.4 .3 .2 .1 [.457 .4565]/2];
      b_vec       = [.4 .3 .2 .1];
      save(filename,'fraxn_vec','b_vec','-append');
   end
else
   if DO2==0
      load(filename);
   end

   A  = .5;
   B  = .5;
   AB = [A,B];
   a  = A/2;
   %%
   np       = 50;
   Npolys   = 10;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%prams for 2-ring res;
   theta_vec   = (0:np-1)'/(np-1)*90;
      %% rotation of 2-ring resonator;
   fraxn_vec   = [.9 .8 .7 .6];
      %% fraction of circle;
   ecc      = 1.5;
         %% b/a = eccentricity of ellipse => circular arcs;
   radratio = .5;
      %% (radius of inner ring)/(radius of outer ring);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   AB2   = [.5 1];
   %%prams for sine wave;
   phi_vec  = (0:np-1)'/(np-1);
      %% phase of sine wave;
   %b_vec    = [.4 .3 .2 .1  .20875 .2075];
   b_vec    = [.4 .3 .2 .1];
      %% amplitude of sine wave;

   %%
   %  m_y=0*Ach_Li;
   for j=1:np
      disp([num2str(j),' runs done, out of ',num2str(np)]);
      th  = theta_vec(j);
      phi = phi_vec(j);
      %%
      for r=1:length(fraxn_vec)
         crk_fxn     = @CURVEprof_circarc;
         crk_prams   = {fraxn_vec(r)};%% fraction of circle
         radius      = a;%% outer radius;
         rot         = 90+th;%% rotation of outer circ (opening points left if rot=90);
         srt         = {radius*[1 ecc],rot,[0 0]};
         srt2        = {radratio*radius*[1 ecc],rot+180,[0 0]};
%        srt2        = {radratio*radius*[1 ecc],270,[0 0]};
         %%
         Irr_vars = {crk_fxn,crk_prams,srt;
                      crk_fxn,crk_prams,srt2};
         %%
         C_eff = BONE_GF_cracks_rect_cell(AB,Irr_vars,Npolys);
         Mhat  = C_eff{1};
         %%
         [U,Lam]     = eig(Mhat);
         m1(j,r,1)   = min(diag(Lam));
         m2(j,r,1)   = max(diag(Lam));
         %%
         jj          = find(diag(Lam)==m1(j,r,1));
         pdir(j,r,1) = angle(U(1,jj)+1i*U(2,jj));
         if j>1
            if pdir(j,r,1)-pdir(j-1,r,1)>.5*pi
               pdir(j,r,1)   = pdir(j,r,1)-pi;
            elseif pdir(j,r,1)-pdir(j-1,r,1)<-.5*pi
               pdir(j,r,1)   = pdir(j,r,1)+pi;
            end
         else
            if pdir(j,r,1)<0
               pdir(j,r,1) = pdir(j,r,1)+pi;
            end
         end
      end%%r loop
      for r=1:length(b_vec)
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if DO2
            %%
            crk_fxn     = @CURVEprof_cos;
            crk_prams   = {1,0};%% fraction of cos wave & phase;
            crk_prams2  = {1,phi};%% fraction of cos wave & phase;
            b0          = .25;
            b           = b_vec(r);
            %%
            srt  = {[a b0],0,[0 .5]};
            srt2 = {[a b],0,[0 -.5]};
            %%
            Irr_vars  = {crk_fxn,crk_prams,srt;
                         crk_fxn,crk_prams2,srt2};
            %%
            C_eff = BONE_GF_cracks_rect_cell(AB2,Irr_vars,Npolys);
            Mhat  = C_eff{1};
            %%
            [U,Lam]     = eig(Mhat);
            m1(j,r,2)   = min(diag(Lam));
            m2(j,r,2)   = max(diag(Lam));
            %%
            jj          = find(diag(Lam)==m1(j,r,2));
            pdir(j,r,2) = angle(U(1,jj)+1i*U(2,jj));
            if j>1
               if pdir(j,r,2)-pdir(j-1,r,2)>.5*pi
                  pdir(j,r,2)   = pdir(j,r,2)-pi;
               elseif pdir(j,r,2)-pdir(j-1,r,2)<-.5*pi
                  pdir(j,r,2)   = pdir(j,r,2)+pi;
               end
            else
               if pdir(j,r,2)<0
                  pdir(j,r,2) = pdir(j,r,2)+pi;
               end
            end
            %%
         end
      end%%r loop
   end

   %% SAVE RESULTS:
   save([filename,'.mat'],'fraxn_vec','b_vec','theta_vec','phi_vec','m1','m2','pdir');
end

vs    = {theta_vec/180,phi_vec};
xlab  = {'th/pi','phi'};

NN = [length(fraxn_vec) length(b_vec)];
for r=1:2
   jj = 1:NN(r);
   %%
   subplot(2,3,1+3*(r-1));
   GEN_set_plotorder;
   plot(vs{r},m1(:,jj,r),'k');
   GEN_proc_fig(xlab{r},'m1');
   %%
   subplot(2,3,2+3*(r-1));
   GEN_set_plotorder;
   plot(vs{r},m2(:,jj,r),'k');
   GEN_proc_fig(xlab{r},'m2');

   %%
   subplot(2,3,3+3*(r-1));
   GEN_set_plotorder;
   pdir(:,jj,r)  = pdir(:,jj,r)-pi/2;
   if r==2
      jp = [1 2];
      pdir(:,jp,r)  = pdir(:,jp,r)+pi;
      %plot(vs{r},pdir(:,5:6,r)/pi,'r');
      %pdir(1,5:6,r)/pi
      %b_vec(5:6)
      hold on;
   end
   plot(vs{r},pdir(:,jj,r)/pi,'k');
   GEN_proc_fig(xlab{r},'chi/pi');
end
%%
if 1
%  %% POLISH FIGURES AND SAVE THEM;
   lets  = {'(a)','(b)','(c)','(d)','(e)','(f)'};
   Y0 = [.4 .6 .4 .6 .7 .4];
   Y1 = [.55 .85 1 .8 .9 1];
   if 1
      Y1([3 6])   = Y1([3 6]) -.5;
      Y0([3 6])   = Y0([3 6]) -.5;
   end
   X1 = [.5 .5 .5 1 1 1];
   ax = .1;
   ay = [.8 .7 .8 .15 .8 .8];

   for j=1:6
      subplot(2,3,j), hold on;
      box on;
      x  = ax*X1(j);
      y  = Y0(j)+ay(j)*(Y1(j)-Y0(j));
      txt   = text(x,y,lets{j});
      ylim([Y0(j),Y1(j)]);
      set(txt,'FontName','Times','FontSize',16);
      hold off;
   end
%  subplot(1,2,1);
%  xlim([0 .5]);
%  GEN_proc_fig('\ita/\itA','\itm^*_{\ity}');
%  %%
%  subplot(1,2,2);
%  xlim([0 .5]);
%  GEN_proc_fig('\ita/\itA','\itm^*_{\itx}');
   %%
   %saveas(gcf,[filename,'.fig'])
   saveas(gcf,[filename,'.eps'])
   !gv out/crack_2cracks.eps &
end

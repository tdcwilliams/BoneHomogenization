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
filename = 'out/crack_4cracks';
DO2      = 1;

if exist([filename,'.mat'])
   load(filename)
else
   if DO2==0
      load(filename);
   end

   A        = 1;
   B        = 1;
   AB       = [A,B];
   a        = .25;
   a4_vec   = [.4 .3 .2 .1];
   %%
   np       = 50;
   Npolys   = 10;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%prams for rotating crack:
   theta_vec   = (0:np-1)'/(np-1)*90;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   %%prams for sine wave
   phi_vec  = (0:np-1)'/(np-1);
      %% phase of sine wave;
   b        = .25;
   %%
   %  m_y=0*Ach_Li;
   for j=1:np
      disp([num2str(j),' runs done, out of ',num2str(np)]);
      th  = theta_vec(j);
      phi = phi_vec(j);
      %%
      for r=1:length(a4_vec)
         a4 = a4_vec(r);%%length of 4th crack;
         %%
         crk_fxn     = @CURVEprof_strtline;
         crk_prams   = {[]};%% fraction of circle
         radius      = a;%% outer radius;
         srt         = {[a  1],0  ,[0.5 0.5];...
                        [a  1],0  ,[-0.5 0.5];...
                        [a  1],0  ,[-0.5 -0.5];...
                        [a4 1],th,[0.5 -0.5]...
                        };
         %%
         for s = 1:4
            Irr_vars(s,:)  = {crk_fxn,crk_prams,srt(s,:)};
            %Irr_vars{s,1},Irr_vars{s,2},Irr_vars{s,3},pause
         end
         %Irr_vars,disp('paused'),pause
         %%
         C_eff = BONE_GF_cracks_rect_cell(AB,Irr_vars,Npolys);
         Mhat  = C_eff{1};
         %%
         [U,Lam]     = eig(Mhat);
         m1(j,r,1)   = min(diag(Lam));
         m2(j,r,1)   = max(diag(Lam));
         %%
         jj          = find(diag(Lam)==m1(j,r,1));
         jj          = jj(1);
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

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if DO2
            %%
            crk_fxn4 = @CURVEprof_cos;
            crk_prams4  = {1,phi};%% fraction of cos wave & phase;
            a4          = a4_vec(r);
            %%
            srt         = {[a  1],0,[0.5 0.5];...
                           [a  1],0,[-0.5 0.5];...
                           [a  1],0,[-0.5 -0.5];...
                           [a4 b],0,[0.5 -0.5]...
                           };
            %%
            for s = 1:4
               Irr_vars(s,:)  = {crk_fxn,crk_prams,srt(s,:)};
               %if s<4,Irr_vars{s,1},Irr_vars{s,2},Irr_vars{s,3},pause,end;
            end
            Irr_vars{4,1}  = crk_fxn4;
            Irr_vars{4,2}  = crk_prams4;
            %Irr_vars{4,1},Irr_vars{4,2},Irr_vars{4,3},pause
            %Irr_vars,disp('paused'),pause
            %%
            C_eff = BONE_GF_cracks_rect_cell(AB,Irr_vars,Npolys);
            Mhat  = C_eff{1};
            %%
            [U,Lam]     = eig(Mhat);
            m1(j,r,2)   = min(diag(Lam));
            m2(j,r,2)   = max(diag(Lam));
            %%
            jj          = find(diag(Lam)==m1(j,r,2));
            jj          = jj(1);
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
      end
   end

   %% SAVE RESULTS:
   save([filename,'.mat'],'theta_vec','phi_vec','m1','m2','pdir');
end

vs    = {theta_vec/180,phi_vec};
xlab  = {'th/pi','phi'};

for r=1:2
   subplot(2,3,1+3*(r-1));
   GEN_set_plotorder;
   plot(vs{r},m1(:,:,r),'k');
   GEN_proc_fig(xlab{r},'m1');
   %%
   subplot(2,3,2+3*(r-1));
   GEN_set_plotorder;
   plot(vs{r},m2(:,:,r),'k');
   GEN_proc_fig(xlab{r},'m2');

   %%
   subplot(2,3,3+3*(r-1));
   GEN_set_plotorder;
%  if r==2
%     pdir(:,1:2,r)  = pdir(:,1:2,r)+pi;
%  end
   plot(vs{r},pdir(:,:,r)/pi-0.5,'k');
   GEN_proc_fig(xlab{r},'chi/pi');
end
%%
if 1
%  %% POLISH FIGURES AND SAVE THEM;
   lets  = {'(a)','(b)','(c)','(d)','(e)','(f)'};
   Y0 = [.75 .85 .5 .7 .9 .5];
   Y1 = [.9 1.01 .7 .9 .96 .56];
   if 1
      Y1([3 6])   = Y1([3 6]) -.5;
      Y0([3 6])   = Y0([3 6]) -.5;
   end
   X1 = [.5 .5 .5 1 1 1];
   ax = .1;
   ay = [.85 .75 .85 .85 .15 .85];

   for j=1:6
      subplot(2,3,j), hold on;
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
end

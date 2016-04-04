%% fig_crack_circle_lim.m
%% plots (a) m_x and (b) m_y for elliptical arc cracks
%%  as the fraction of the ellipse covered tends to 1;
%% we check that m_x & m_y tend to those of an elliptical cavity
%%  (m_xy=0 here as axes are aligned with sides of unit cell);
%% eccentricity of ellipse ('ecc') is 0.25;

%% RECT STUFF FROM FIG_CRACK_CIRCLE_LIM.M
filename = 'out/crack_circle_lim.mat'
load(filename)
%%
subplot(2,2,1);
GEN_set_plotorder;
%plot(aoA_vec,m_y,'k'), hold on;
plot(aoA_vec,m1,'k'), hold on;
box on;
%GEN_proc_fig('\ita/\itA','\itm^*_{\ity}');
GEN_proc_fig('a/A','my');
%GEN_proc_fig('a/A','m1');
X  = [0 .5];
xlim(X);
Y  = [.75 1.01];
ylim(Y);
ax = .15;
ay = .2;
x  = X(1)+ax*(X(2)-X(1));
%%
hold on;
y     = Y(1)+ay*(Y(2)-Y(1));
txt   = text(x,y,'(a)');
set(txt,'FontName','times','FontSize',16);
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,2);
GEN_set_plotorder;
%plot(aoA_vec,m_x,'k')
plot(aoA_vec,m2,'k')
box on;
%GEN_proc_fig('\ita/\itA','\itm^*_{\itx}');
GEN_proc_fig('a/A','mx');
xlim(X);
Y  = [.92 1.005];
ylim(Y);
%%
hold on;
y     = Y(1)+ay*(Y(2)-Y(1));
txt   = text(x,y,'(b)');
set(txt,'FontName','times','FontSize',16);
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HEX STUFF;
filename2   = 'out/circle_lim_v2';

%% to use a different solver and/or change the unit cell shape;
cav_solver  = @BONE_MPG_cavs_hexcell;
crk_solver  = @BONE_MPG_cracks_hexcell;
%%
B     = .5;
D     = 2*B/sqrt(3);%%equilateral hexagon;
rho   = .5;%%E = rho*D

unitcell_prams = [B D rho];

%% outer 'bubble';
crk_fxn_out    = @CURVEprof_circarc;
crk_prams_out  = {1};%% fraction of circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

np       = 50;
aoB_vec  = (1:np)'/np*.5;
   %% a = major axis of ellipse (parallel to x axis);

if ~exist([filename2,'.mat'])
   ecc  = 0.25;
      %% b/a = eccentricity of ellipse;

   %%
   Npolys      = 10;
   fraxn_vec   = [.9 .7 .5];
   NMP         = 15;
   Nterms      = 15;
   %  m_y=0*Ach_Li;
   for j=1:np
      a              = aoB_vec(j)*B;
      rout           = .5*(a+B);
      srt_out        = {rout*[1 1],0,[0 0]};
      Irr_vars_outer = {crk_fxn_out,crk_prams_out,srt_out};
      disp([num2str(j),' configurations done (out of ',num2str(np),')']);
      %%
      crk_fxn            = @CURVEprof_circarc;
      crk_prams          = {1};%% fraction of circle
      radius             = a;
      srt                = {radius*[1 ecc],0,[0 0]};
      Irr_vars_inner{1}  = {crk_fxn,crk_prams,srt};
      %%
      C_eff     = BONE_MPG_cavs_hexcell(unitcell_prams,Irr_vars_outer,...
                     Irr_vars_inner,NMP,Nterms);
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
         crk_fxn           = @CURVEprof_circarc;
         crk_prams         = {fraxn_vec(r)};%% fraction of circle
         radius            = a;
         srt               = {radius*[1 ecc],0,[0 0]};
         Irr_vars_inner{1} = {crk_fxn,crk_prams,srt};
         %%
         C_eff     = BONE_MPG_cracks_hexcell(unitcell_prams,Irr_vars_outer,...
                        Irr_vars_inner,NMP,Nterms);
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
   save(filename2,'aoB_vec','m1','m2','pdir');
end

subplot(2,2,3);
GEN_set_plotorder;
%plot(aoA_vec,m_y,'k'), hold on;
plot(aoB_vec,m1,'k'), hold on;
box on;
%GEN_proc_fig('\ita/\itA','\itm^*_{\ity}');
GEN_proc_fig('a/B','my');
%GEN_proc_fig('a/A','m1');
xlim(X);
Y  = [.75 1.01];
ylim(Y);
%%
hold on;
y     = Y(1)+ay*(Y(2)-Y(1));
txt   = text(x,y,'(c)');
set(txt,'FontName','times','FontSize',16);
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,4);
GEN_set_plotorder;
%plot(aoA_vec,m_x,'k')
plot(aoB_vec,m2,'k')
box on;
%GEN_proc_fig('\ita/\itA','\itm^*_{\itx}');
GEN_proc_fig('a/B','mx');
xlim(X);
Y  = [.92 1.005];
ylim(Y);
%%
hold on;
y     = Y(1)+ay*(Y(2)-Y(1));
txt   = text(x,y,'(d)');
set(txt,'FontName','times','FontSize',16);
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveas(gcf,'out/crack_circle_lim_v2.eps')

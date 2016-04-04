%% fig_crack_lim_v2.m
%% elliptical cavities (major axis=a)
%% -> straight crack (length 2a) as eccentricity is increased;
%%
%% also can compare Achenbach & Li (1986)
%% & and Nemat-Nasser (**) results for straight cracks;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% rectangular stuff;
filename    = 'out/crack_lim';
filename3   = 'out/crack_lim_v3';
load([filename,'.mat']);

subplot(1,2,1);
GEN_set_plotorder;
plot(aoA_vec,m_y,'k'), hold on;
%GEN_proc_fig('\ita/\itA','\itm^*_{\ity}');
GEN_proc_fig('a/A','my');
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


subplot(1,2,2);
GEN_set_plotorder;
plot(aoA_vec,m_x,'k')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% hex stuff;
%filename2   = 'out/crack_lim_v2';
%
%%% to use a different solver and/or change the unit cell shape;
%cav_solver  = @BONE_MPG_cavs_hexcell;
%crk_solver  = @BONE_MPG_cracks_hexcell;
%%%
%B     = .5;
%D     = 2*B/sqrt(3);%%equilateral hexagon; (sides=.5)
%rho   = .5;%%E = rho*D
%
%unitcell_prams = [B D rho];
%
%%% outer 'bubble';
%crk_fxn        = @CURVEprof_circarc;
%crk_prams      = {1};%% fraction of circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%np       = 50;
%aoB_vec  = (1:np)'/np*.5;
%
%
%if ~exist([filename2,'.mat'])
%   Npolys   = 10;
%   NMP      = 15;
%   ecc_vec  = [.05,.1,.25];
%   Nvec     = [20 25 35];
%   %  m_y=0*Ach_Li;
%   for j=1:np
%      a              = aoB_vec(j)*B;
%      rout           = .5*(a+B);
%      srt_out        = {rout*[1 1],0,[0 0]};
%      Irr_vars_outer = {crk_fxn,crk_prams,srt_out};
%      disp([num2str(j),' configurations done (out of ',num2str(np),')']);
%      %%
%      Irr_vars_inner{1}  = CURVEget_strtline(-a,a);
%      C_eff              = feval(crk_solver,unitcell_prams,Irr_vars_outer,...
%                            Irr_vars_inner,NMP,Npolys);
%      m_y(j,1)           = C_eff{1}(4);%[m_y(j),Ach_Li(j)]
%      m_x(j,1)           = 1;
%      %%
%      for r=1:length(ecc_vec)
%         crk_fxn           = @CURVEprof_circarc;
%         crk_prams         = {1};%% fraction of circle
%         radius            = a;
%         ecc               = ecc_vec(r);
%         srt               = {radius*[1 ecc],0,[0 0]};
%         Irr_vars_inner{1} = {crk_fxn,crk_prams,srt};
%         %%
%         Nterms     = Nvec(r);
%         C_eff      = feval(cav_solver,unitcell_prams,Irr_vars_outer,...
%                            Irr_vars_inner,NMP,Npolys);
%         m_y(j,r+1) = C_eff{1}(4);
%         m_x(j,r+1) = C_eff{1}(1);
%      end
%   end
%   save([filename2,'.mat'],'aoB_vec','unitcell_prams',...
%         'm_x','m_y');
%else
%   load([filename2,'.mat']);
%   %m_y,m_x
%end
%
%%%
%subplot(2,2,3);
%GEN_set_plotorder;
%plot(aoB_vec,m_y,'k'), hold on;
%%GEN_proc_fig('\ita/\itA','\itm^*_{\ity}');
%GEN_proc_fig('a/B','my');
%xlim([0 .5]);
%Y  = [.75 1.01];
%ylim(Y);
%%%
%hold on;
%y     = Y(1)+ay*(Y(2)-Y(1));
%txt   = text(x,y,'(c)');
%set(txt,'FontName','times','FontSize',16);
%hold off;
%%%
%subplot(2,2,4);
%GEN_set_plotorder;
%plot(aoB_vec,m_x,'k')
%%GEN_proc_fig('\ita/\itA','\itm^*_{\itx}');
%GEN_proc_fig('a/B','mx');
%xlim(X);
%Y  = [.92 1.005];
%ylim(Y);
%%%
%hold on;
%y     = Y(1)+ay*(Y(2)-Y(1));
%txt   = text(x,y,'(d)');
%set(txt,'FontName','times','FontSize',16);
%hold off;

%GEN_setsize_eps([],16,7);
saveas(gcf,[filename3,'.eps']);
%!gv out/crack_lim_v2.eps

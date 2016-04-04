%% fig_crack_lim.m
%% elliptical cavities (major axis=a)
%% -> straight crack (length 2a) as eccentricity is increased;
%%
%% also can compare Achenbach & Li (1986)
%% & and Nemat-Nasser (**) results for straight cracks;


filename = 'out/crack_lim';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% to use a different solver and/or change the unit cell shape;
cav_solver  = @BONE_GF_cavs_rect_cell;
crk_solver  = @BONE_GF_cracks_rect_cell;
%%
A              = .5;
B              = .5;
AB             =[A,B];
unitcell_prams = AB;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

np       = 50;
aoA_vec  = (1:np)'/np*.5;

%% Achenbach & Li (1986);
Ach_Li   = 1./(1-2*A/pi/B*log(cos(pi/2*aoA_vec)));

%% Nemat-Nasser et al (1993);
s10   = 0;
for n=1:5000
   v1    = n*pi*aoA_vec;
   v2    = n*pi*B/A;
   t1    = 2*besselj(1,v1)./v1;
   t2    = v2/tanh(v2);

   %% sum defined as over \Zfield-{0}, but
   %% I sum over \Nfield (since it is even in n),
   %% so need to multiply by 2;
   s10   = s10+2*t2*t1.^2;
end
Nem_Nasser  = 1./(1+1./s10);

if 0
   plot(aoA_vec,[Ach_Li,Nem_Nasser]);
   return;
elseif 0
   subplot(1,2,1)
   n0 = 1;
   JJ = n0:2*n0:np;
   plot(aoA_vec(JJ),Ach_Li(JJ),'.'), hold on;
   plot(aoA_vec(JJ),Nem_Nasser(JJ),'o');
end


if ~exist([filename,'.mat'])
   Npolys   = 10;
   ecc_vec  = [.05,.1,.25];
   Nvec     = [20 25 35];
   %  m_y=0*Ach_Li;
   for j=1:np
     a         = aoA_vec(j)*A;
     disp([num2str(j),' configurations done (out of ',num2str(np),')']);
     %%
     Irr_vars  = CURVEget_strtline(-a,a);
     C_eff     = feval(crk_solver,unitcell_prams,Irr_vars,Npolys);
     m_y(j,1)  = C_eff{1}(4);%[m_y(j),Ach_Li(j)]
     m_x(j,1)  = 1;
     %%
     for r=1:length(ecc_vec)
       crk_fxn    = @CURVEprof_circarc;
       crk_prams  = {1};%% fraction of circle
       radius     = a;
       ecc        = ecc_vec(r);
       srt        = {radius*[1 ecc],0,[0 0]};
       Irr_vars   = {crk_fxn,crk_prams,srt};
       %%
       Nterms     = Nvec(r);
       C_eff      = feval(cav_solver,unitcell_prams,Irr_vars,Nterms);
       m_y(j,r+1) = C_eff{1}(4);
       m_x(j,r+1) = C_eff{1}(1);
     end
   end
   save([filename,'.mat'],'aoA_vec','unitcell_prams',...
         'm_x','m_y','Ach_Li','Nem_Nasser');
else
   load([filename,'.mat']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,1);
GEN_set_plotorder;
plot(aoA_vec,m_y,'k'), hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2);
GEN_set_plotorder;
plot(aoA_vec,m_x,'k')
%%
if 1
  subplot(1,2,1);
  xlim([0 .5]);
  %GEN_proc_fig('\ita/\itA','\itm^*_{\ity}');
  GEN_proc_fig('a/A','my');
  %%
  subplot(1,2,2);
  xlim([0 .5]);
  %GEN_proc_fig('\ita/\itA','\itm^*_{\itx}');
  GEN_proc_fig('a/A','mx');
  %%
  %saveas(gcf,[filename,'.fig']);
  GEN_setsize_eps([],16,7);
  saveas(gcf,[filename,'.eps']);
  %%
  %!gv out/crack_lim.eps &
end

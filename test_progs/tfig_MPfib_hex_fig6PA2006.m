%% tfig_MPfib_hex_fig6pA2006.m

rvec  = (.05:.05:.45)';
nr    = length(rvec);
%%
m1vec = [10 5 2 1 .5 .1 .01];%% m1=10 is value in fig 6 of PA06 paper;
nm    = length(m1vec);
mx    = zeros(nr,nm);
%%
srt   = {[],0,[0 0]};
%%

%BDrho  = [sqrt(3)/4, 1/2,.5];
%BDrho  = [sqrt(3)/2,1/2,.5];
BDrho  = BONE_hexconvert(1/sqrt(3));
Nterms = 35;
for j=1:nr
   for s = 1:nm
      radius      = rvec(j);
      m1          = m1vec(s);
      midi        = [m1 3];
      srt{1}      = radius*[1 1];
      Irr_vars    = {@CURVEprof_circarc,{1},srt,midi};
      %%
      Ceff        = BONE_MP_fibres_hexcell(BDrho,Irr_vars,Nterms);
      mx(j,s)     = Ceff{1}(1);[m1,radius],Ceff{1}
   end
   Irr_vars    = Irr_vars(1:3);
   Ceff        = BONE_MP_cavs_hexcell(BDrho,Irr_vars,Nterms);
   mxcav(j,1)  = Ceff{1}(1);
   %%
   disp([j nr]);
end

plot(rvec,mx);
hold on;
plot(rvec,mxcav,'xk');
%%
%%eg points from wp06

%%m_1=10
plot([.2 .3 .4 .45],[1.4 1.8 2.6 4],'x');
RR    = .1*(1:4)';
m_eff = [1.05
         1.3
         1.8
         3];
plot(RR,m_eff,'or');

%%m_1=0.1
RR    = .1*(1:4)';
m_eff = [0.94235
         0.787755
         0.578358
         0.355193];
plot(RR,m_eff,'vr');

%% m=10^(-8) (WP's equivalent "cavity")
RR    = .1*(1:4)';
m_eff = [0.929988
         0.746566
         0.507712
         0.264145];
plot(RR,m_eff,'^r');


GEN_proc_fig('r_1','m^*');
saveas(gcf,'out/tfig_MPfib_hex_fig6PA2006.eps')
%plot(rvec,mx,'k');

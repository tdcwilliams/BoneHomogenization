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
filename = 'out/crack_2ringres';


if exist([filename,'.mat'])
   load(filename)
else

   A  = .5;
   B  = .5;
   AB = [A,B];
   np = 50;
   %%
   aoA_vec  = (1:np)'/np*.5;
      %% a = major axis of ellipse (parallel to x axis);
   ecc      = 1;
      %% b/a = eccentricity of ellipse;
   radratio = .5;
      %% (radius of inner ring)/(radius of outer ring)

   %%
   Npolys      = 10;
   fraxn_vec   = [.9 .8 .7 .6];
   Nterms      = 15;
   %  m_y=0*Ach_Li;
   for j=1:np
     a   = aoA_vec(j)*A;[j np]
     %%
     for r=1:length(fraxn_vec)
        crk_fxn   = @CURVEprof_circarc;
        crk_prams = {fraxn_vec(r)};%% fraction of circle
        radius    = a;%% outer radius;
        rot       = 90;%% rotation of outer circ (opening points up);
        srt       = {radius*[1 ecc],rot,[0 0]};
        srt2      = {radratio*radius*[1 ecc],rot+180,[0 0]};
           %% inner circle has radius radratio*radius & opening points down;

        Irr_vars  = {crk_fxn,crk_prams,srt;
                        crk_fxn,crk_prams,srt2};
        %%
        C_eff     = BONE_GF_cracks_rect_cell(AB,Irr_vars,Npolys);
        m_y(j,r)  = C_eff{1}(4);
        m_x(j,r)  = C_eff{1}(1);
        m_xy(j,r) = C_eff{1}(2);
     end
   end

   %% SAVE RESULTS:
   save([filename,'.mat'],'aoA_vec','m_x','m_y','m_xy');
end

subplot(1,2,1);
GEN_set_plotorder;
plot(aoA_vec,m_y,'k'), hold on;
%%
subplot(1,2,2);
GEN_set_plotorder;
plot(aoA_vec,m_x,'k')
%%
if 1
   %% POLISH FIGURES AND SAVE THEM;
   subplot(1,2,1);
   xlim([0 .5]);
   GEN_proc_fig('\ita/\itA','\itm^*_{\ity}');
   %%
   subplot(1,2,2);
   xlim([0 .5]);
   GEN_proc_fig('\ita/\itA','\itm^*_{\itx}');
   %%
   saveas(gcf,[filename,'.fig'])
   saveas(gcf,[filename,'.eps'])
end

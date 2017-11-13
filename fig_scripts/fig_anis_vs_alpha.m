%% fig_comp_achli_nemnas.m
%% - change 'height' of unit (rectangular) cell;
%% - compare to Delameter et al (1975);
%% - check effect of orthotropic host medium;

clear;
outdir   = 'out'
if ~exist(outdir,'dir')
   mkdir(outdir);
end
filename = 'out/anis_vs_alpha';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% to use a different solver and/or change the unit cell shape;
crk_solver  = @BONE_GF_cracks_rect_cell;
%%
A     = .5;
Bhat  = .5;
ahat  = A/2;
%%
na       = 100;
alp_vec  = linspace(0.1,2,na)';
%%
nt          = 4;
theta_vec   = linspace(0,pi/2,nt)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1%~exist([filename,'.mat'])
   Npolys   = 50;
   m1       = zeros(na,nt);
   m2       = m1;
   pdir     = m1;
   %%
   for j=1:na
      for r=1:nt
         alp      = alp_vec(j);
         th_hat   = theta_vec(r);
         Dalp     = diag([1 1/alp]);
         %%
         B  = Bhat*alp;
         th = atan(alp*tan(th_hat));
         a  = ahat*sqrt(cos(th_hat)^2+alp^2*sin(th_hat)^2);
         z  = a*exp(1i*th);
         AB = [A,B];
         %%
         Irr_vars    = CURVEget_strtline(-z,z);
         C_eff       = feval(crk_solver,AB,Irr_vars,Npolys);
         Mhat        = Dalp*C_eff{1}*Dalp;
         %%
         [U,Lam]  = eig(Mhat);
         m1(j,r)  = min(diag(Lam));
         m2(j,r)  = max(diag(Lam));
         %%
         jj          = find(diag(Lam)==m1(j,r));
         pdir(j,r)   = angle(U(1,jj)+1i*U(2,jj));
         %%
         if j>1
            if pdir(j,r)-pdir(j-1,r)>.5*pi
               pdir(j,r)   = pdir(j,r)-pi;
            elseif pdir(j,r)-pdir(j-1,r)<-.5*pi
               pdir(j,r)   = pdir(j,r)+pi;
            end
         end
      end
      disp([num2str(j),' configurations done (out of ',num2str(na),')']);
   end
   save([filename,'.mat'],'alp_vec','theta_vec','A','Bhat','ahat','m1','m2','pdir');
else
   load([filename,'.mat']);
end

ylab  = {'m1','m2','chi/pi'};
yy    = {m1,m2,pdir/pi};
%%
for j=1:3
   subplot(1,3,j);
   GEN_set_plotorder;
   plot(alp_vec,yy{j},'k');
   GEN_proc_fig('alp',ylab{j});
end
saveas(gcf,[filename,'.eps'],'epsc');

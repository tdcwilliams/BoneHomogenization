%% fig_comp_achli_nemnas.m
%% - change 'height' of unit (rectangular) cell;
%% - compare to Delameter et al (1975);
%% - check effect of orthotropic host medium;

clear;
filename = 'out/anis_vs_rotang';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% to use a different solver and/or change the unit cell shape;
crk_solver  = @BONE_GF_cracks_rect_cell;
%%
A     = .5;
Bhat  = .5;
ahat  = A/2;
%%
nt          = 100;
theta_vec   = linspace(0,pi/2,nt)';
%%
alp_vec  = [.5 .75 1 1.25 1.5 1.75 2];
na       = length(alp_vec);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist([filename,'.mat'])
   Npolys   = 50;
   m1       = zeros(nt,na);
   m2       = m1;
   pdir     = m1;
   %%
   for j=1:nt
   %for j=nt
      for r=1:na
         th_hat   = theta_vec(j);
         alp      = alp_vec(r);
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
         if 0%alp==1
            [th_hat,th]
            [A,Bhat,ahat]
            [A,B,a]
            C_eff{1}
            Mhat
            [m1(j,r),m2(j,r),pdir(j,r)]
            pause
         end
         if j>1
            if pdir(j,r)-pdir(j-1,r)>.5*pi
               pdir(j,r)   = pdir(j,r)-pi;
            elseif pdir(j,r)-pdir(j-1,r)<-.5*pi
               pdir(j,r)   = pdir(j,r)+pi;
            end
         end
      end
      disp([num2str(j),' configurations done (out of ',num2str(nt),')']);
   end
   save([filename,'.mat'],'alp_vec','theta_vec','A','Bhat','ahat','m1','m2','pdir');
else
   load([filename,'.mat']);
end

ylab  = {'m1','m2','chi/pi'};
yy    = {fliplr(m1),fliplr(m2),fliplr(pdir)/pi};
%%
for j=1:3
   subplot(1,3,j);
   GEN_set_plotorder;
   plot(theta_vec/pi,yy{j}(:,1:4),'k');
   hold on;
   plot(theta_vec/pi,yy{j}(:,5:end),'r');
   GEN_proc_fig('th/pi',ylab{j});
end
saveas(gcf,[filename,'.eps'],'epsc');

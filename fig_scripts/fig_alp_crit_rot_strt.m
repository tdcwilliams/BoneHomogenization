%%fig_alp_crit_rot_strt.m

outdir   = 'out'
if ~exist(outdir,'dir')
   mkdir(outdir);
end
fname = [outdir,'/alp_crit_rot_strt.mat'];
if exist(fname)
   load(fname);
   %rot_vec
else
   aoA      = .5;
   rot_vec  = 90:-2:0;
   Nr       = length(rot_vec);
   alpc_vec = 0*rot_vec;
   fmin_vec = 0*rot_vec;

   for j=1
      rot                  = rot_vec(j);
      [alpc,fvalue,Mhat]   = BONE_alpcrit_strt_rot(aoA,rot);
      alpc_vec(j)          = alpc;
      fmin_vec(j)          = fvalue;
   end

   for j=2:Nr
      rot                  = rot_vec(j);
      [alpc,fvalue,Mhat]   = BONE_alpcrit_strt_rot(aoA,rot,alpc);
      alpc_vec(j)          = alpc;
      fmin_vec(j)          = fvalue;
      GEN_progrep([j Nr]);
   end
   save(fname,'alpc_vec','fmin_vec','rot_vec');
end


%xfac     = 1;
xfac     = 1/180;
rot_vec  = xfac*rot_vec;
%%
[ax,hh{1},hh{2}]  = plotyy(rot_vec,alpc_vec,rot_vec,fmin_vec);

if BW==0
   lc = {[0 0 0],[1 0 0]};
else
   lc = {[0 0 0],.35*[1 1 1]};
end
ls = {'-','--'};
lw = {1,1.5};
if 0
   yl = {'ac','Del'};
else
   yl = {'ac',''};
end

if xfac==1
   XX = [0 90];
else
   XX = [0 0.5];
end

for j=1:2
   yyl   = get(ax(j),'Ylabel');
   set(yyl,'String',yl{j},'fontname','times','fontsize',16);
   if j==1
      xxl1  = get(ax(j),'Xlabel');
      set(xxl1,'String','th/pi','fontname','times','fontsize',16);
      xlp1  = get(xxl1,'position');
      set(ax(j),'xtick',[0 .25 .5]);
   else
      set(ax(j),'xtick',[]);
   end
   %set(get(ax(j),'fontname'),'String','times');
   set(ax(j),'fontname','times','fontsize',16,'ycolor',lc{j},...
         'xlim',XX);
   %set(han,'string');
   %set(,'fontname','times','fontsize',16);
   %set(ax(j),'FontName','times','fontsize',16,...
      %'xlabel','theta, deg','ylabel',yl{j});
   set(hh{j},'color',lc{j},'linestyle',ls{j},'linewidth',lw{j});
   %%
   %set(xl,'FontName','times','fontsize',16);
end

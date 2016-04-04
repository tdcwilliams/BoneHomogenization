kh    = 0:.1:2*pi;
aod   = .25;
hod   = 0.5;
kd    = kh/hod;
%%
C        = -pi./kd/log(cos(pi*aod));
cos_qh   = cos(kh)-1./C.*sin(kh);
%%
qh = acos(cos_qh);
slope = (kh(2)-kh(1))/(kh(2)-kh(1))

plot(qh/pi,kh/pi);
saveas(gcf,['out/achli-full-hod',num2str(hod),'.eps'])

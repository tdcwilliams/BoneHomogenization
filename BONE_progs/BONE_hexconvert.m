function BDrho = BONE_hexconvert(ABE)
%% CALL: BDrho = BONE_hexconvert(ABE)
%% periodicity of 2 rectangular lattices is [A,B]
%% side length is 2E

do_test  = 0;
if nargin==0
   %%equilateral of side length sl
   sl = .5;
   %%
   th = (0:5)'/6*2*pi;
   x  = sl*cos(th);
   y  = sl*sin(th);
   plot(x,y,'o');
   %%
   D     = x(1);
   E     = x(2);
   A     = D+E;
   B     = y(2);
   ABE   = [A B E];
   %%
   do_test  = 1;
elseif length(ABE)==1
   %%equilateral of side length sl
   sl = ABE;
   %%
   th = (0:5)'/6*2*pi;
   x  = sl*cos(th);
   y  = sl*sin(th);
   %%
   D     = x(1);
   E     = x(2);
   A     = D+E;
   B     = y(2);
   ABE   = [A B E];
end

B  = ABE(2);
E  = ABE(3);
%%
A     = ABE(1);%[A,sqrt(3)/2],[A,B]
D     = A-E;
rho   = E/D;
%%
BDrho = [B,D,rho];

if do_test
   hold on;
   x1 = [D,E,-E,-D,-E,E,D]';
   y1 = [0,B,B,0,-B,-B,0]';
   plot(x1,y1);
   hold off;
end

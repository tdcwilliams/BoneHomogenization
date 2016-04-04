function F=BONE_multipoles_laplace(...
             X,Y,M,centre)

DO_TEST=0;
if nargin==0%%do test
  DO_TEST=1;
  XTEST=0;
  AB=[.5 .5];
  M=6;
  centre=[0 0];
  %%
  Nint=200;
  [tq,wq]=GEN_numint_exp(Nint);
  if XTEST
    X=tq*AB(1);
    Y=.2+0*X;
  else
    Y=tq*AB(2);
    X=.2+0*Y;
  end
end

Z=X-centre(1)+1i*(Y-centre(2));
jj=1:M;
[JJ,ZZ]=meshgrid(jj,Z);
F=ZZ.^JJ;
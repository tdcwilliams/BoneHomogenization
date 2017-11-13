function addpaths_matlibdir(matlibdir)

if ~exist('matlibdir','var')
   matlibdir   = [getenv('MATLIBDIR')];
   if strcmp('', matlibdir)
      matlibdir   = pwd;
   end
end
matlibdir   = [matlibdir,'/'];

addpath([matlibdir,'BES_progs/'   ])
addpath([matlibdir,'CURVE_progs/' ])
addpath([matlibdir,'GEN_progs/'   ])
addpath([matlibdir,'OP_progs/'    ])
addpath([matlibdir,'SF_progs/'    ])
addpath([matlibdir,'STAT_progs/'  ])

Switch=1;%%use addpath (not rmpath)

%% DEFINE MASTER DIRECTORY:
%% farm:
%  master_directory='/home/matdcw/WORK/programs/matlab';
%  Windows=0;

if 0%% optimus:
  master_directory='/maybehome/twilliams/MATHS2010/programs';
  Windows=0;
  SS='usr';
end

if 1%% Nansen:
  master_directory='/Home/timill/MATHS-NansenWork/programs';
  Windows=0;
  SS='opt';
end

%% SPECIFY WHICH SUBDIRECTORIES ARE REQUIRED:
sub_directories{1}='GEN_progs';
sub_directories{end+1}='ND_progs';
%  sub_directories{end+1}='RTS_progs';
%  sub_directories{end+1}='WT_progs';
%  sub_directories{end+1}='RDG_progs';
%  sub_directories{end+1}='MF_progs';
%  sub_directories{end+1}='BIG_FILES';
%  sub_directories{end+1}='WH_progs';
%  sub_directories{end+1}='GRN_progs';
sub_directories{end+1}='GEN_progs/OP_progs';
sub_directories{end+1}='CURVE_progs';
%  sub_directories{end+1}='EIG_progs';

%% ADD/REMOVE SPECIFIED DIRECTORIES:
if Switch==1
  for j=1:length(sub_directories)
    addpath([master_directory,'/',sub_directories{j}]);
  end
else
  for j=1:length(sub_directories)
    rmpath([master_directory,'/',sub_directories{j}]);
  end
end

if Windows==0;%% print user-defined path directories:
  D=path;
  j=[0,find(D==':')];
  r=1;
%  crit=( D(j(r)+2:j(r)+4)~='usr' );
  crit=( D(j(r)+2:j(r)+4)~=SS );
  if crit==0
    disp('There are no user-defined path directories.')
  else
    disp('User-defined path directories:');
    while crit
      jj=j(r)+1:j(r+1)-1;
      disp(D(jj));
      r=r+1;
      crit=( D(j(r)+2:j(r)+4)~=SS );
    end
  end
end

clear;

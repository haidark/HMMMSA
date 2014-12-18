program hmmstrgca
  !-----------------------------------------------------------------------
  use hmmstr
  use database
  implicit none
  !-----------------------------------------------------------------------
  integer,parameter::atom_max=50000
  real,parameter::prob_cut=1.0
  real,parameter::radius=3.0
  real,parameter::d2_cut=radius*radius
  !-----------------------------------------------------------------------
  character(len=150)::modelfile
  character(len=150)::seqfile
  character(len=150)::profile
  character(len=150)::pdbfile
  character(len=150)::dbfile
  character(len=150)::outfile
  character(len=1)::chain
  character(len=1)::str_flag
  integer::iarg
  integer::natom
  integer::ios
  integer::nres
  logical::use_str
  real,dimension(3,atom_max)::structure
  real,dimension(:,:),allocatable::xgamma
  real,dimension(:,:,:,:,:),allocatable::v_hbonds
  real,dimension(:,:),allocatable  ::calpha
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  integer::ires,t,k,iatom
  integer::xlow,xhigh,ylow,yhigh,zlow,zhigh
  integer,dimension(3)::binmax
  real::x,d2,xscale
  real,dimension(3)::xyzmax
  real,dimension(3)::inframe
  !-----------------------------------------------------------------------
  ! Parse the command line
  !-----------------------------------------------------------------------
  iarg = command_argument_count()
  !-----------------------------------------------------------------------
  if(iarg<3) then
    write(*,*)"-----------------------------------------------------------------------------------------------------------"
    write(*,*) &
    'Usage: xhmmstrgca model_R.hmm seqfile(fasta) gcafile profile(psi-blast) use_structure(HMMstr T/F) [pdbfile chain] '
    write(*,*)"-----------------------------------------------------------------------------------------------------------"
    write(*,*)"==> Output in GCA format. "
    write(*,*)"-----------------------------------------------------------------------------------------------------------"
    stop 'hmmstrgca.f90 v.  Fri Jun 25 12:46:00 EDT 2010'
  endif
  !-----------------------------------------------------------------------
  pdbfile=" "; chain=" ";str_flag="F";profile=" "
  call getarg(1,modelfile)
  call getarg(2,seqfile)
  call getarg(3,outfile)
  if (iarg>=4) &
  call getarg(4,profile)
  if (iarg>=5) &
  call getarg(5,str_flag)
  read(str_flag,'(l1)',iostat=ios) use_str
  if (ios/=0) then
    write(0,*) 'hmmstrgca: str_flag=',trim(str_flag)
    stop 'hmmstrgca: Bad value for arg 5: use_structure'
  endif
  !-----------------------------------------------------------------------
  if (iarg>=6) &
  call getarg(6,pdbfile)
  if (iarg>=7) &
  call getarg(7,chain)
  !-----------------------------------------------------------------------
  call hmmstr_gamma(xgamma,nres,modelfile,seqfile,profile, &
               use_str,pdbfile,chain,gca=outfile)
  !-----------------------------------------------------------------------
  if(allocated(v_hbonds))    deallocate(v_hbonds)
  if(allocated(xgamma))      deallocate(xgamma)
  !-----------------------------------------------------------------------
end program hmmstrgca

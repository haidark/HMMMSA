module hmmstr
  implicit none
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  public::hmmstr_gamma
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  private::idatom
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  private::crossprod
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  private::read_profile
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  private::read_seq
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  private::readmodel
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  private::get_outtrans
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  private::get_intrans
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  private::flagrat_brdc_prod
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  private::bprof_bg
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  private::ramatype
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  private::getalpha
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  private::getbeta
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  private::read_backbone
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  private::release_model
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  private::set_flags
  !-----------------------------------------------------------------------------
  private
   !-------------------------------------------------------------
   integer,parameter::maxseqlen=5000
   integer,parameter::m=20
   integer,parameter::mrama=11
   integer,parameter::mdssp=6
   integer,parameter::mctxt=10
   integer,parameter::smin=3
   integer,parameter::smax=27
   real,parameter::ncount=1.5
   real,parameter::eps=0.01
   real,parameter::alphaeps=0.01
   !-------------------------------------------------------------
   type model
     integer::N
     real::bground(m)
     real::rground(mrama)
     real::dground(mdssp)
     real::cground(mctxt)
     real,dimension(:),allocatable::prior
     real,dimension(:,:),allocatable::dih
     character(len=1),dimension(:),allocatable::ss
     real,dimension(:,:),allocatable::a
     real,dimension(:,:),allocatable::b
     real,dimension(:,:),allocatable::r
     real,dimension(:,:),allocatable::d
     real,dimension(:,:),allocatable::c
     real,dimension(:,:),allocatable::logb
     character(len=1),dimension(6)::ss_map
   end type 
   !-------------------------------------------------------------
   type bb
     real,dimension(3,4)::xyz
   end type
   !-------------------------------------------------------------
   type trans
    integer,dimension(:),allocatable::c
   end type
   !-------------------------------------------------------------
   type(model)::modelR
   type(trans),dimension(:),allocatable::intrans,outtrans
   real::sumd(mdssp),sumr(mrama),sumc(mctxt),gsum
   !-------------------------------------------------------------
   real,dimension(:,:),allocatable  :: calpha
   real,dimension(:,:),allocatable  :: angles  
 contains
!------------------------------------------------------------------------------------------------------------
  subroutine  &
  hmmstr_gamma(xgamma,nres,mfile,seqfile,pfile,use_str,pdbfile,chain,gca)
   !---------------------------------------------
   implicit none
   !---------------------------------------------
   character(len=*),intent(in)::mfile
   character(len=*),intent(in)::seqfile
   character(len=*),intent(in)::pfile
   character(len=*),intent(in)::pdbfile
   character(len=1),intent(inout)::chain
   integer,intent(out)::nres
   logical,intent(in)::use_str
   real,dimension(:,:),allocatable,intent(inout)::xgamma
   character(len=*),optional :: gca
   !---------------------------------------------
   character(len=1),dimension(:),allocatable::ramaseq
   ! character(len=3000)::seq
   integer,dimension(3000) ::seq, aaseq
   integer::t
   integer::i
   integer::drctFlag
   integer::iatom
   integer::ires
   integer::drctlen
   integer,dimension(:),allocatable::iflag
   real,dimension(:),allocatable::ct
   real,dimension(:,:),allocatable::backbone
   real,dimension(:,:),allocatable::profile
   real,dimension(:,:),allocatable,target::angles
   real,dimension(:,:),allocatable::alpha
   real,dimension(:,:),allocatable::beta
   type(bb),dimension(:),allocatable::bbxyz
   !---------------------------------------------
   call readmodel(mfile,modelR,sumd,sumr,sumc)
   call get_intrans(intrans,modelR)
   call get_outtrans(outtrans,modelR)
   !---------------------------------------------
   if(allocated(bbxyz)) deallocate(bbxyz)
   if(allocated(profile)) deallocate(profile)
   if(allocated(backbone)) deallocate(backbone)
   if(allocated(calpha))  deallocate(calpha)
   if(allocated(angles))  deallocate(angles)
   !---------------------------------------------
   aaseq = 0
   !check if we are dealing with a drctfile
   i = index(seqfile, ".")
   if (trim(seqfile(i:))==".drct") then
     !make sure not to get profile again
     drctFlag = 1
     call getdrctlen(seqfile, drctlen)
     allocate(profile(20,drctlen))
     profile = 0.
     call readdrct(seqfile, drctlen, profile, seq, nres)
   else 
     drctFlag = 0
     call read_seq(seqfile,seq,nres)
     allocate(profile(20,nres))
   endif
   !! diagnostic
   !write(*,'(20i4)') seq(1:nres)
   !---------------------------------------------
   allocate(calpha(3,nres)) ; calpha=999.
   allocate(angles(3,nres)) ; angles=0.
   allocate(backbone(3,4*nres))
   allocate(bbxyz(nres))
   
   if(use_str.and.pdbfile/=" ") then
     aaseq = seq
     call read_backbone(pdbfile,chain,bbxyz,seq,nres)
     if (any(seq(1:nres)/=aaseq(1:nres))) then
       write(0,*) "Sequence from fasta file does not match PDB file"
       write(*,*) "============= FASTA sequence ==============="
       write(*,'(20i4)') aaseq(1:nres)
       write(*,*) "============= PDB sequence ==============="
       write(*,'(20i4)') seq(1:nres)
       stop
     endif
     do ires=1,nres
       calpha(1:3,ires) = bbxyz(ires)%xyz(1:3,2)
     enddo
     iatom=1
     do ires=1,nres
       backbone(:,iatom:iatom+3) = bbxyz(ires)%xyz(:,1:4)
       iatom = iatom+4
     enddo
     call get_angles(backbone,nres,angles)
     call read_profile(pfile,profile,seq,nres)
   else
    if (drctFlag/=1) then
     if (pfile/=" ") then
       profile = 0.
       call read_profile(pfile,profile,seq,nres)
     else
       profile = 0.
       do ires=1,nres
         if (seq(ires)<=20) profile(seq(ires),ires) = 1.00
       enddo
     endif
    endif
   endif
   write(*,'(a,11f8.4)')"Background RAMA ",modelR%rground
   write(*,'(a,6f8.4)')"Background DSSP ",modelR%dground
   write(*,'(a,10f8.4)')"Background CTXT ",modelR%cground
   if(allocated(iflag))   deallocate(iflag)
   if(allocated(alpha))   deallocate(alpha)
   if(allocated(ct))      deallocate(ct)
   if(allocated(beta))    deallocate(beta)
   if(allocated(xgamma))  deallocate(xgamma)
   if(allocated(ramaseq)) deallocate(ramaseq)
   allocate(ct(nres))
   allocate(iflag(nres))
   allocate(ramaseq(nres))
   allocate(alpha(modelR%n,nres))
   allocate(beta(modelR%n,nres))
   allocate(xgamma(modelR%n,nres))
   call set_flags(ramaseq,nres,use_str,angles,iflag)
   call getalpha(ct,alpha,nres,modelR,profile,intrans,ramaseq,iflag)
   call getbeta(ct,beta,nres,modelR,profile,outtrans,ramaseq,iflag)
   do t=1,nres
     do i=1,modelR%N
       xgamma(i,t)=alpha(i,t)*beta(i,t)/ct(t)
     enddo
   enddo
   !!------------------------------------------------------------------------------------------
   if (present(gca)) &
     call write_gca(gca,mfile,nres,seq,angles,calpha,modelR%N,xgamma,modelR,profile)
   !!---------------------------------------------------------------------------------------------
   if(allocated(bbxyz))    deallocate(bbxyz)
   if(allocated(iflag))    deallocate(iflag)
   if(allocated(alpha))    deallocate(alpha)
   if(allocated(ct))       deallocate(ct)
   if(allocated(beta))     deallocate(beta)
   if(allocated(ramaseq))  deallocate(ramaseq)
   if(allocated(profile))  deallocate(profile)
   if(allocated(backbone)) deallocate(backbone)
   if(allocated(angles))   deallocate(angles)
   if(allocated(bbxyz)) deallocate(bbxyz)
   !------------------------------------------------------------------------------------------------
   call release_model
   if(allocated(intrans))    deallocate(intrans)
   if(allocated(outtrans))   deallocate(outtrans)
   if(allocated(calpha))  deallocate(calpha)
   if(allocated(angles))  deallocate(angles)
  end subroutine hmmstr_gamma
!-------------------------------------------------------------------------------------------------------------
  subroutine readdrct(filename, length, profile, seq, nres)
     implicit none
     character(*), intent(in)::filename
     integer, intent(in)::length
     real*4, dimension (length, 20) ,intent(out)::profile
     integer, intent(out)::nres
     character(len=3000)::seqch
     integer,dimension(3000),intent(out) ::seq
     character(len=21),parameter :: aa1="ACDEFGHIKLMNPQRSTVWYX"
     integer*4::i, j, status=0
     integer*2::seq_list, resSeq
     integer, parameter::recsize=136
     real*4::phi_list, psi_list, omega_list, csum, gap_freq, ins_freq
     character::structure_name*4, chainID*1, seq_resid*1, ss*1, rotamer*1, domainID*1, contex*1
     integer*2, dimension(8)::clusterID
     real*4, dimension(3)::xyzca
     real*4, dimension(20)::vall_pro
     open(unit=99, file=filename, status="old", form="unformatted", access="direct", recl=recsize)
     i=0
     do
     	i=i+1
     	read(99, rec=i, iostat=status) seq_list, resSeq, structure_name, chainID, seq_resid, ss, rotamer,&
        	(xyzca(j),j=1, 3), phi_list, psi_list, omega_list, csum, (vall_pro(j), j=1, 20), gap_freq, ins_freq,&
        	domainID, contex, (clusterID(j),j=1, 3)
     	if (status/=0) exit
           profile(i, 1:20)=vall_pro
	   seqch(i:i)=seq_resid
	   
     end do
     close(99)
     nres = i-1
     do i=1,nres
        seq(i) = index(aa1,seqch(i:i))
     enddo
     return
     end subroutine readdrct
!-------------------------------------------------------------------------------------------------------------
subroutine getdrctlen(filename, length)
     implicit none
     character(*), intent(in)::filename
     integer, intent(out)::length
     integer*4::i, j, status=0
     integer*2::seq_list, resSeq
     integer, parameter::recsize=136
     real*4::phi_list, psi_list, omega_list, csum, gap_freq, ins_freq
     character::structure_name*4, chainID*1, seq_resid*1, ss*1, rotamer*1, domainID*1, contex*1
     integer*2, dimension(8)::clusterID
     real*4, dimension(3)::xyzca
     real*4, dimension(20)::vall_pro
     open(unit=99, file=filename, status="old", form="unformatted", access="direct", recl=recsize)
     i=0
     do
       i=i+1
       read(99, rec=i, iostat=status) seq_list, resSeq, structure_name, chainID, seq_resid, ss, rotamer,&
           (xyzca(j),j=1, 3), phi_list, psi_list, omega_list, csum, (vall_pro(j), j=1, 20), gap_freq, ins_freq,&
           domainID, contex, (clusterID(j),j=1, 3)
       if (status/=0) exit
     end do
     length=i-1
     close(99)
     return
    end subroutine getdrctlen
 
!-------------------------------------------------------------------------------------------------------------
  subroutine write_gca(gcafile,mfile,nres,seq,angles,calpha,nnode,xgamma,modelR,profile)
     character(len=*),intent(in) :: gcafile, mfile
     integer,intent(in) :: nres, nnode
     integer,dimension(nres),intent(in) :: seq
     real,dimension(3,nres),intent(in) :: angles, calpha
     real,dimension(nnode,nres),intent(in) :: xgamma
     type(model),intent(in)::modelR
     real,dimension(20,nres),intent(inout)::profile
     !---
     character :: aachar
     character(len=20),parameter :: aa="ACDEFGHIKLMNPQRSTVWY"
     integer :: i,ios,gunit=29,ires,iaa,inode
     character(len=6) :: sschar="HEGST_"
     real,dimension(20) :: gprofile
     real,dimension(mdssp) :: ss
     !---
     open(gunit,file=trim(gcafile),status='replace',form='formatted',iostat=ios)
     if (ios/=0) then
       write(0,*) "iostat=",ios
       stop 'hmmstrgca:: error opening output gca file'
     endif
     write(gunit,'(a)') 'HMM: '//trim(mfile)//' prob_flag= 1 num_nodes= 282'
     do ires=1,nres
       !! Write RESIDUE line
       aachar="X"
       if (seq(ires)>0.and.seq(ires)<=20) aachar=aa(seq(ires):seq(ires))
       write(gunit,'(a,2i4,1x,a1,$)') 'RESIDUE: ',ires,ires,aachar
       write(gunit,'(3f9.1,3f9.2,a)') angles(1:3,ires),calpha(1:3,ires),' . . .'
       !! Write SS line
       ss(1:mdssp) = 0.
       do inode=1,nnode
         ss(1:mdssp) = ss(1:mdssp) + xgamma(inode,ires)*modelR%d(inode,1:mdssp)
       enddo
       ss(1:mdssp) = ss(1:mdssp)/sum(ss(1:mdssp))
       iaa = maxloc(ss(1:mdssp),dim=1)
       aachar = sschar(iaa:iaa)
       write(gunit,'("SS: ",a1,6f9.6)') aachar,ss(1:mdssp)
       !! Write PROFILE line
       !profile = 0.
       !profile(seq(ires)) = 1.000
       write(gunit,'("PROFILE:",20f9.6)') profile(1:20,ires)
       !! Write GPROFILE line
       gprofile(1:20) = 0.
       do inode=1,nnode
         gprofile(1:20) = gprofile(1:20) +  &
           xgamma(inode,ires)*modelR%b(inode,1:20)
       enddo
       gprofile(1:20) = gprofile(1:20)/sum(gprofile(1:20))
       write(gunit,'("GPROFILE:",20f9.6)') gprofile
       !! Write GAMMA line
       write(gunit,'("GAMMA:",999f9.6)') xgamma(1:nnode,ires)
     enddo
     close(gunit)
  end subroutine write_gca
!-------------------------------------------------------------------------------------------------------------
  subroutine release_model()
   implicit none
   if(allocated(modelR%prior)) deallocate(modelR%prior)
   if(allocated(modelR%dih))   deallocate(modelR%dih)
   if(allocated(modelR%ss))    deallocate(modelR%ss)
   if(allocated(modelR%a))     deallocate(modelR%a)
   if(allocated(modelR%b))     deallocate(modelR%b)
   if(allocated(modelR%r))     deallocate(modelR%r)
   if(allocated(modelR%d))     deallocate(modelR%d)
   if(allocated(modelR%c))     deallocate(modelR%c)
   if(allocated(modelR%logb))  deallocate(modelR%logb)
  end subroutine release_model
!------------------------------------------------------------------------------------------------------------
  subroutine set_flags(ramaseq,nres,use_str,angles,iflag)
   !----------------------------------------
   implicit none
   !----------------------------------------
   character(len=1),dimension(:),intent(out)::ramaseq
   integer,intent(in)::nres
   integer,dimension(:),intent(out)::iflag
   logical,intent(in)::use_str
   real,dimension(:,:),intent(in)::angles
   !----------------------------------------
   integer::i
   !----------------------------------------
   ramaseq="?"
   iflag=1
   if(use_str) then
     do i=1,nres
        if(all(angles(:,i)/=999.)) then
           iflag(i)=3
        endif
        if(iflag(i)==3) then
           ramaseq(i)=ramatype(angles(1,i),angles(2,i),angles(3,i))
        endif
     enddo
   endif
  end subroutine set_flags
!------------------------------------------------------------------------------------------------------------
  subroutine getbeta(ct,beta,seqres_len,modelR,profile,otrs,ramaseq,iflag)
   !--------------------------------------------
   implicit none
   !--------------------------------------------
   character(len=1),dimension(:),intent(in)::ramaseq
   integer,intent(in)::seqres_len
   integer,dimension(:),intent(in)::iflag
   real,dimension(:),intent(inout)::ct
   real,dimension(:,:),intent(out)::beta
   real,dimension(:,:),intent(in)::profile
   type(model),intent(in)::modelR
   type(trans),dimension(:),intent(in)::otrs
   !--------------------------------------------
   integer::i
   integer::j
   integer::t
   integer::nq
   integer::jnd
   real(kind=8)::drat
   real(kind=8)::rrat
   real(kind=8)::ss
   real(kind=8)::cct
   !--------------------------------------------
   !!beta(nq,nres)  !!outtrans(nr)%c(nonzero)! profile(a(1:20),nres)
   beta=0
   nq=modelR%n
   do i=2,nq
     beta(i,seqres_len)=ct(seqres_len)
   enddo
   ss=0
   do j=1, otrs(1)%c(0)
     jnd = otrs(1)%c(j)
     ss = ss + modelR%a(1,jnd)*beta(jnd,seqres_len)* &
          flagrat_brdc_prod(iflag(seqres_len),jnd,seqres_len,profile,modelR,ramaseq(seqres_len))
   enddo
   beta(1,seqres_len)=ss
   do t=seqres_len-1,1,-1
      cct=0
      do i=2,nq
        ss=0
        do j=1,otrs(i)%c(0)
           jnd = otrs(i)%c(j)
           ss = ss + modelR%a(i,jnd)*beta(jnd,t+1)*flagrat_brdc_prod(iflag(t),jnd,t+1,profile,modelR,ramaseq(t+1))
        enddo
        beta(i,t)=ss*ct(t)
        cct = cct+beta(i,t)
      enddo
      ss=0
      do j=1,otrs(1)%c(0)
        jnd = otrs(1)%c(j)
        ss = ss + modelR%a(1,jnd)*beta(jnd,t)*flagrat_brdc_prod(iflag(t),jnd,t,profile,modelR,ramaseq(t))
      enddo
      beta(1,t)=ss
      if(cct==0)then
         write(*,*)"NO out going transistions!"
         stop
      endif
    enddo
  end subroutine
!-------------------------------------------------------------------------------------------------------------
  subroutine getalpha(ct,alpha,seqres_len,modelR,profile,intrs,ramaseq,iflag)
   !--------------------------------------------------
   implicit none
   !--------------------------------------------------
   character(len=1),dimension(:),intent(in)::ramaseq
   integer,intent(in)::seqres_len
   integer,dimension(:),intent(in)::iflag
   real,dimension(:,:),intent(out)::alpha
   real,dimension(:,:),intent(in)::profile
   real,dimension(:),intent(inout)::ct
   type(model),intent(in)::modelR
   type(trans),dimension(:),intent(in)::intrs
   !--------------------------------------------------
   integer::i
   integer::j
   integer::t
   integer::nq
   integer::ind
   real(kind=8)::drat
   real(kind=8)::rrat
   real(kind=8)::ss
   real(kind=8)::junk
   real(kind=8)::tmp
   !--------------------------------------------------
   !!alpha(nq,nres)  !!intrans(nr)%c(nonzero)! profile(a20,nres)
   alpha=0;ct=0
   drat=1.;rrat=1.0
   ct(1)=0
   nq=modelR%n
   do j=1,nq
     alpha(j,1)=modelR%a(1,j) * flagrat_brdc_prod(iflag(1),j,1,profile,modelR,ramaseq(1))
     ct(1) = ct(1) + alpha(j,1)
   enddo
   ct(1) = 1./ct(1)
   do j=1,nq
     alpha(j,1) = alpha(j,1)*ct(1)
   enddo
   j=1
   ss=0
   do i=1, intrs(j)%c(0)
     ind = intrs(j)%c(i)
     ss = ss + alpha(ind,1)*modelR%a(ind,j)
   enddo
   alpha(j,1)=ss
   do t=2,seqres_len
     ct(t)=0.
     do j=2,nq  !! 1 is not node
       ss=0
       do i=1, intrs(j)%c(0)
         ind = intrs(j)%c(i)
         ss = ss + alpha(ind,t-1)*modelR%a(ind,j)
       enddo
       alpha(j,t) = ss * flagrat_brdc_prod(iflag(t),j,t,profile,modelR,ramaseq(t))
       ct(t) = ct(t) + alpha(j,t)
     enddo
     if(ct(t)==0) then
       write(*,*) "No incoming transistions were found!"
       ct=1.
       return
     else
       ct(t) = 1./ct(t)
       do j=2,nq
         alpha(j,t) = alpha(j,t)*ct(t)
       enddo
       j=1
       ss=0
       do i=1,intrs(j)%c(0)
         ind = intrs(j)%c(i)
         ss = ss + alpha(ind,t) * modelR%a(ind,j)
       enddo
       alpha(j,t)=ss
       !write(*,*)t,sum(alpha(:,t))
     endif   
   enddo
  end subroutine
!-------------------------------------------------------------------------------------------------------------
  character(len=1) function ramatype(ph,ps,om)
   !--------------------------------------------
   implicit none
   !--------------------------------------------
   character(len=mrama)::rama_string="HGBEdbeLlxc"
   real,intent(in)::ph
   real,intent(in)::ps
   real,intent(in)::om
   !--------------------------------------------
   integer::i
   integer::icen
   real::d
   real::dmin
   real::ds
   real::df
   real,dimension(11,2)::ramacen
   !--------------------------------------------
   ramacen(1,1)=-61.91   ;ramacen(1,2)=-45.20
   ramacen(2,1)=-109.78  ;ramacen(2,2)=20.88
   ramacen(3,1)=-70.58   ;ramacen(3,2)=147.22
   ramacen(4,1)=-132.89  ;ramacen(4,2)=142.43
   ramacen(5,1)=-135.03  ;ramacen(5,2)=77.26
   ramacen(6,1)=-85.03   ;ramacen(6,2)=72.26
   ramacen(7,1)=-165.00  ;ramacen(7,2)=175.00
   ramacen(8,1)=55.88    ;ramacen(8,2)=38.62
   ramacen(9,1)=85.82    ;ramacen(9,2)=-.03
   ramacen(10,1)=80.     ;ramacen(10,2)=-170.00
   ramacen(11,1)=-70.0   ;ramacen(11,2)=150.00
   if(ph==999.)then; ramatype="?"; return; endif
   if(ps==999.)then; ramatype="?"; return; endif
   if(om<90.and.om>-90)then; ramatype="c"; return; endif
   dmin=999
   do i=1,mrama-1
     df=abs(ph-ramacen(i,1)); if(df>180.) df= 360.-df  
     ds=abs(ps-ramacen(i,2)); if(ds>180.) ds= 360.-ds  
     d=sqrt(ds*ds+df*df)
     if(d<dmin) then
       icen=i; dmin=d
     endif  
   enddo
   ramatype=rama_string(icen:icen)
  end function ramatype
!-------------------------------------------------------------------------------------------------------------
  real function bprof_bg(iq,t,profile,modelR)
   !---------------------------
   implicit none
   !---------------------------
   integer,intent(in)::t
   integer,intent(in)::iq
   real,dimension(:,:),intent(in)::profile
   type(model),intent(in)::modelR
   !---------------------------
   integer::a
   real(kind=8)::logsum
   real(kind=8)::tmp
   !---------------------------
   logsum=0
   do a=1,20
     if(profile(a,t)/=0)then
       tmp = log((modelR%b(iq,a)+(EPS*modelR%bground(a)))/((1.+EPS)*modelR%bground(a)))    
       logsum = logsum + profile(a,t)*ncount*tmp   
     endif
   enddo
   bprof_bg=exp(logsum)
  end function bprof_bg
!-------------------------------------------------------------------------------------------------------------
  real function flagrat_brdc_prod(jflag,iq,t,profile,modelR,rama_c)
   !-----------------------------------------
   implicit none
   !-----------------------------------------
   character(len=1),intent(in)::rama_c
   integer,intent(in)::jflag
   integer,intent(in)::t,iq
   real,dimension(:,:),intent(in)::profile
   type(model),intent(in)::modelR
   !-----------------------------------------
   integer::k
   real::brat
   real::rrat
   character(len=mrama)::rama_string="HGBEdbeLlxc"
   !-----------------------------------------
   brat=1; rrat=1
   if(iq==1) then
     flagrat_brdc_prod=brat
     return
   endif
   if(jflag==3)then
     do k=1,mrama
       if(rama_string(k:k)==rama_c) exit
     enddo
     if( k < mrama+1) then
       rrat = (modelR%r(iq,k) + EPS*modelR%rground(k))/((1.+EPS)*modelR%rground(k))
     !else
     !  write(*,*)"RAMA character not fount! ",rama_c
     !  stop
     endif
   endif  
   if(jflag<=3)then
     brat = bprof_bg(iq,t,profile,modelR)
   endif
   flagrat_brdc_prod = brat*rrat
  end function flagrat_brdc_prod
!-------------------------------------------------------------------------------------------------------------
  subroutine get_intrans(intrs,mdr)
   !---------------------------------------
   implicit none
   !---------------------------------------
   type(model),intent(in)::mdr
   type(trans),dimension(:),allocatable,intent(out)::intrs
   !---------------------------------------
   integer::nq
   integer::ncount
   integer::i
   integer::j
   integer::ios
   !---------------------------------------
   nq=mdr%N
   allocate(intrs(nq))
   do j=1,nq
     ncount=0 
     do i=1,nq
       if(mdr%a(i,j)/=0) ncount = ncount + 1
     enddo
     allocate(intrs(j)%c(0:ncount))
     ncount=0 
     do i=1,nq
       if(mdr%a(i,j)/=0) then
         ncount = ncount + 1
         intrs(j)%c(ncount) = i
       endif
       intrs(j)%c(0) = ncount
     enddo  
   enddo
  end subroutine get_intrans
!-------------------------------------------------------------------------------------------------------------
  subroutine get_outtrans(outtrs,mdr)
   !--------------------------------------------
   implicit none
   !--------------------------------------------
   type(model),intent(in)::mdr
   type(trans),dimension(:),allocatable,intent(out)::outtrs
   !--------------------------------------------
   integer::nq
   integer::ncount
   integer::ios
   integer::i
   integer::j
   !--------------------------------------------
   nq=mdr%N
   allocate(outtrs(nq))
   do i=1,nq
     ncount=0 
     do j=1,nq
       if(mdr%a(i,j)/=0) ncount = ncount + 1
     enddo
     allocate(outtrs(i)%c(0:ncount))
     ncount=0 
     do j=1,nq
       if(mdr%a(i,j)/=0) then
         ncount = ncount + 1
         outtrs(i)%c(ncount) = j
       endif
       outtrs(i)%c(0) = ncount
     enddo  
   enddo
  end subroutine get_outtrans
!-------------------------------------------------------------------------------------------------------------
  subroutine readmodel(mfile,modelR,sumd,sumr,sumc)
   !---------------------------------------
   implicit none
   !---------------------------------------
   character(len=*),intent(in)::mfile
   real,dimension(:),intent(out)::sumd
   real,dimension(:),intent(out)::sumr
   real,dimension(:),intent(out)::sumc
   type(model),intent(out)::modelR
   !---------------------------------------
   character(len=300)::aline
   character(len=200)::junk
   integer::ios
   integer::i
   integer::j
   integer::nnodes
   integer::n
   real::x
   real::z
   !---------------------------------------
 
   open(1,file=mfile,iostat=ios)
   if(ios/=0) then 
      write(*,*) "Cannot open model file!"
      stop
   endif
   sumr=0; sumd=0; sumc=0
   do 
     read(1,'(a)',iostat=ios) aline
     if(ios/=0) then
       write(*,*) "Improperly formatted modelR.hmm"
       stop
     endif 
     if(aline(1:16) == "default_bg_freqs")then
        read(aline(17:197),*)modelR%bground(1:20)
        write(*,'(a,20f9.5)')"BACKGROUND FREQS ",modelR%bground(1:20)
     endif
     if(aline(1:9) == "num_nodes")then
        read(aline(10:15),*)modelR%N
        write(*,'(a,i5)')"NUM NODES ",modelR%N
        exit
     endif
   enddo
   modelR%N = modelR%N + 1
   n=modelR%N
   allocate(modelR%prior(n))
   allocate(modelR%dih(n,3))
   allocate(modelR%ss(n))
   allocate(modelR%a(n,n))
   allocate(modelR%b(n,m))
   allocate(modelR%r(n,mrama))
   allocate(modelR%d(n,mdssp))
   allocate(modelR%c(n,mctxt))
   allocate(modelR%logb(n,m))
   do 
     read(1,'(a)',iostat=ios) aline
     if(ios/=0) then
       write(*,*) "Improperly formatted modelR.hmm"
       stop
     endif 
     if(aline(1:9) == "unk_node ")then
        read(aline(14:207),*)modelR%prior(1)
        read(1,*)modelR%b(1,1:m)
        modelR%dih(1,1)= -75
        modelR%dih(1,2)= -15
        modelR%dih(1,3)= 180
        modelR%ss(1) = "_"
     endif
     if(aline(1:11) == "unk_node_ss") then
        read(aline(19:63),*)modelR%d(1,1:mdssp)
        sumd = sumd + modelR%prior(1)*modelR%d(1,:)
     endif
     if(aline(1:13) == "unk_node_rama") then
        read(aline(21:105),*)modelR%r(1,1:mrama)
        sumr = sumr + modelR%prior(1)*modelR%r(1,:)
        exit
     endif            
   enddo
   !!! For all known nodes
   do i=2,n
     read(1,'(a)',iostat=ios)aline
     if(ios/=0) then
       write(*,*) "Improperly formatted modelR.hmm"
       stop
     endif 
     read(aline(6:8),*)z
     if(i-1/=z .or. aline(1:4) /= "node") then
       write(*,*) "Improperly formatted modelR.hmm"
       stop
     endif 
     read(aline,*)junk,junk,junk,junk,modelR%prior(i),modelR%dih(i,1:3)
     read(1,*)modelR%b(i,1:10)
     read(1,*)modelR%b(i,11:20)
     read(1,'(a)',iostat=ios)aline
     read(aline(21:66),*)modelR%d(i,1:mdssp)
     sumd = sumd + modelR%prior(i)*modelR%d(i,:)
     read(1,'(a)',iostat=ios)aline
     read(aline(24:108),*)modelR%r(i,1:mrama)
     sumr = sumr + modelR%prior(i)*modelR%r(i,:)
     read(1,'(a)',iostat=ios)aline
     read(aline(27:103),*)modelR%c(i,1:mctxt)
     sumc = sumc + modelR%prior(i)*modelR%c(i,:)
   enddo
   read(1,*)aline
   if(aline(1:13)/="transit_freqs") then
     write(*,*) "Improperly Formatted modelR.hmm"
   endif
   x=0; x = sum(sumr)
   modelR%rground = sumr/x
   x=0; x = sum(sumd)
   modelR%dground = sumd/x
   x=0; x = sum(sumc)
   modelR%cground = sumc/x
   do i=1,modelR%N
     read(1,*)modelR%a(i,1:modelR%n)
     !write(*,*)modelR%a(i,1:n)
   enddo
   close(1)
  end subroutine readmodel
!-------------------------------------------------------------------------------------------------------------
  subroutine read_seq(sfile,seq,seqnres)
   !-------------------------------------
   implicit none
   !-------------------------------------
   integer,intent(out)::seqnres
   character(len=*),intent(in)::sfile
   character(len=3000)::seqch
   integer,dimension(3000),intent(out) ::seq
   !-------------------------------------
   character(len=200)::aline
   character(len=21),parameter :: aa1="ACDEFGHIKLMNPQRSTVWYX"
   integer::ios,ires
   !-------------------------------------
   open(33,file=sfile,iostat=ios)
   if(ios/=0)then
     write(*,*)"Cannot open seq file!"
     stop
   endif
   read(33,*)aline
   if(aline(1:1)/=">")then
      write(*,*)"Sequence file should be fasta format"
      stop
   endif
   read(33,*)aline
   seqch = trim(aline)
   do 
     read(33,*,iostat=ios)aline
     if(ios/=0)exit
     seqch = trim(seqch)//trim(aline)
   enddo
   seqnres=len_trim(seqch)
   do ires=1,seqnres
     seq(ires) = index(aa1,seqch(ires:ires))
   enddo
   close(33)
  end subroutine read_seq
!-------------------------------------------------------------------------------------------------------------
  subroutine read_profile(pfile,profile,seq,nres)
   !--------------------------------------
   implicit none
   character(len=*),intent(in)::pfile
   character(len=3000)::seqch
   integer,dimension(3000),intent(in)::seq
   integer,intent(in)::nres
   real,dimension(20,nres),intent(inout)::profile
   !--------------------------------------
   character(len=300)::aline
   integer::i
   integer::j
   integer::ios, ires
   real,dimension(nres)::tmpvec
   character(len=1),dimension(21)::res1=(/&
     'A','C','D','E','F','G','H','I','K','L','M','N',&
     'P','Q','R','S','T','V','W','Y','X'/)
   character(len=21),parameter :: aa1="ACDEFGHIKLMNPQRSTVWYX"
   !--------------------------------------
   do ires=1,nres
     seqch(ires:ires) = res1(seq(ires))
   enddo
   write(*,*) "Sequence=",seqch(1:nres)
   open(1,file=pfile,form='formatted',status='old',iostat=ios)
   if(ios/=0)then
     write(*,*) "subroutine read_profile (hmmstr.f95) :: cannot open profile!"
     stop
   endif
   read(1,'(a)') aline
   do while(aline(12:21)/="A  R  N  D")
     read(1,'(a)') aline
   enddo
   profile=0.
   do i=1,nres
     read(1,'(a)',iostat=ios)aline
     if(ios/=0) exit
     if (aline(1:5)=="     ")exit
     if(aline(7:7)/=seqch(i:i)) then
        write(*,*)i,aline(7:7),seqch(i:i)
        write(*,*)"subroutine read_profile (hmmstr.f95) :: profile and seqfile nres do not match! (1)"
        write(*,*)"nres i = ",nres, i
        stop
     endif
     read(aline(72:150),*)profile(1:20,i)
   enddo
   if(i>nres+1) then
     write(*,*)"subroutine read_profile (hmmstr.f95) :: profile and seqfile nres do not match! (2)"
     write(*,*)"nres i = ",nres, i
     stop
   endif
   close(1)
   profile=profile/sum(profile)
   tmpvec = profile(2,:); profile(2,:)=profile(5,:); profile(5,:)=tmpvec
   tmpvec = profile(3,:); profile(3,:)=profile(4,:); profile(4,:)=tmpvec
   tmpvec = profile(4,:); profile(4,:)=profile(7,:); profile(7,:)=tmpvec
   tmpvec = profile(5,:); profile(5,:)=profile(14,:); profile(14,:)=tmpvec
   tmpvec = profile(6,:); profile(6,:)=profile(8,:); profile(8,:)=tmpvec
   tmpvec = profile(7,:); profile(7,:)=profile(9,:); profile(9,:)=tmpvec
   tmpvec = profile(8,:); profile(8,:)=profile(10,:); profile(10,:)=tmpvec
   tmpvec = profile(9,:); profile(9,:)=profile(12,:); profile(12,:)=tmpvec
   tmpvec = profile(10,:); profile(10,:)=profile(11,:); profile(11,:)=tmpvec
   tmpvec = profile(11,:); profile(11,:)=profile(13,:); profile(13,:)=tmpvec
   tmpvec = profile(13,:); profile(13,:)=profile(15,:); profile(15,:)=tmpvec
   tmpvec = profile(14,:); profile(14,:)=profile(15,:); profile(15,:)=tmpvec
   tmpvec = profile(18,:); profile(18,:)=profile(19,:); profile(19,:)=tmpvec
   tmpvec = profile(20,:); profile(20,:)=profile(18,:); profile(18,:)=tmpvec
   do i=1,nres
     if(all(profile(1:20,i)==0))then
        do j=1,21
           if(res1(j)==seqch(i:i)) exit
        enddo
        if(j/=21) then
           profile(j,i)=1.0
        else
           profile(1:20,i)= modelR%bground(1:20)
        endif
     endif
     ! write(*,*)profile(1:20,i)
   enddo
  end subroutine read_profile
!-------------------------------------------------------------------------------------------------------------
  subroutine get_angles(backbone,nres,angles)
   !---------------------------------------------------
   implicit none
   integer,intent(in)::nres
   real,dimension(:,:),intent(in)::backbone
   real,dimension(:,:),intent(out)::angles
   !---------------------------------------------------
   integer::i
   integer::j
   integer::iatom
   real,dimension(3,4)::phixyz
   real,dimension(3,4)::psixyz
   real,dimension(3,4)::omgxyz
   !---------------------------------------------------
   angles=999.
   i=1
   psixyz(:,1:3) = backbone(1:3,i:i+2) 
   psixyz(:,4) = backbone(1:3,i+4) 
   omgxyz(:,1:2) = backbone(1:3,i+1:i+2) 
   omgxyz(:,3:4) = backbone(1:3,i+4:i+5) 
   angles(2,1) = cal_tors(psixyz)
   angles(3,1) = cal_tors(omgxyz)
   j=2
   do i=5,(nres-1)*4,4
     phixyz(:,1) = backbone(:,i-2) 
     phixyz(:,2:4) = backbone(:,i:i+2) 
     psixyz(:,1:3) = backbone(1:3,i:i+2) 
     psixyz(:,4) = backbone(1:3,i+4) 
     omgxyz(:,1:2) = backbone(1:3,i+1:i+2) 
     omgxyz(:,3:4) = backbone(1:3,i+4:i+5) 
     angles(1,j) = cal_tors(phixyz)
     angles(2,j) = cal_tors(psixyz)
     angles(3,j) = cal_tors(omgxyz)
     j = j + 1
   enddo
   phixyz(:,1) = backbone(:,i-2) 
   phixyz(:,2:4) = backbone(:,i:i+2)
   angles(1,j) = cal_tors(phixyz)
  end subroutine get_angles
!-------------------------------------------------------------------------------------------------------------
  real function cal_tors(xyz) 
   !----------------------------------
   implicit none
   !----------------------------------
   real,dimension(3,4),intent(in)::xyz
   !----------------------------------
   integer::i
   real::angle
   real::handcheck
   real::len_c13
   real::len_c24
   real,dimension(3)::vec
   real,dimension(3)::c13
   real,dimension(3)::c24
   real,dimension(3)::c13xc24
   real,dimension(3)::rotaxis
   real,dimension(3,4)::tmpxyz
   real,parameter::pi=3.14159265
   real,parameter::rad2deg=180/pi
   !----------------------------------
   vec=xyz(:,2)
   do i=1,4
     tmpxyz(:,i) = xyz(:,i) - vec
   enddo
   call crossprod(tmpxyz(:,1),tmpxyz(:,3),c13)
   rotaxis=tmpxyz(:,3)
   vec=xyz(:,3)
   do i=1,4
     tmpxyz(:,i) = xyz(:,i) - vec
   enddo
   call crossprod(tmpxyz(:,2),tmpxyz(:,4),c24)
   call crossprod(c24,c13,c13xc24)
   handcheck = sum(c13xc24*rotaxis)
   len_c13 = sqrt(sum(c13*c13));  len_c24 = sqrt(sum(c24*c24))
   angle = sum(c13*c24)/(len_c13*len_c24)
   if(angle>1) angle=1;   if(angle < -1) angle = -1 !! Safety first.
   angle = acos(angle)
   if(handcheck>0) then
     angle = -angle 
   endif
   cal_tors = angle*rad2deg
  end function cal_tors
!-------------------------------------------------------------------------------------------------------------
  subroutine read_backbone(pdbfile,chain,bbxyz,seq,nres)
   !----------------------------------------------------
   implicit none
   !----------------------------------------------------
   character(len=*),intent(in)::pdbfile
   integer,dimension(3000),intent(out)::seq
   character(len=3000)::seqch
   character(len=1),intent(inout)::chain
   integer,intent(inout)::nres
   type(bb),dimension(nres),intent(inout)::bbxyz
   !----------------------------------------------------
   character(len=1)::altloc
   character(len=100)::aline
   integer::i
   integer::j
   integer::k
   integer::iatom
   integer::ios
   integer::ires,mres
   integer::atype
   character(len=6)::last
   character(len=9)::curres
   character(len=9)::qres
   integer,parameter :: NAA=21
   character(len=3),dimension(NAA)::three=(/'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE',&
     'LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL',&
     'TRP','TYR','CRO'/)
   character(len=1),dimension(NAA)::res1=(/'A','C','D','E','F','G','H','I','K','L','M','N',&
     'P','Q','R','S','T','V','W','Y','X'/)
   !----------------------------------------------------
   open(2,file=pdbfile,status='old',iostat=ios)
   if(ios/=0)then
     write(*,*) "Cannot open pdbfile ",trim(pdbfile),"!"
     stop
   endif
   ires=0;curres=" ";seqch=" ";seq=0
   do 
     read(2,'(a)',iostat=ios) aline
     if(ios/=0) exit
     if(aline(1:6)=='ENDMDL')exit              !! use first model of NMR structures!
     if(aline(1:3)=="TER".and.ires/=0)exit     !! Atoms with the correct chain id but not part of the chain! 
     if(aline(1:6)/='ATOM  '.and.aline(1:6)/='HETATM') cycle
     if(aline(22:22)/=chain.and.aline(22:22)/=" ")cycle
     if(aline(22:22)==" ")chain="_"
     ! if(aline(13:16)/=' CA ')cycle
     if(aline(18:26)==curres)cycle             !! Only one CA per residue allowed!
     !! diagnostic
     ! write(*,'(a)') trim(aline)
     curres=aline(18:26)
     ires = ires + 1
     k=NAA+1
     do j=1,NAA
       if(aline(18:20)==three(j)) then
          k=j; exit
       endif
     enddo
     seqch(ires:ires) = res1(k)
     seq(ires) = k
   enddo
   mres=ires; iatom=0; ires=0;
   do ires=1,mres
     do iatom=1,4
       bbxyz(ires)%xyz(1:3,iatom)=999.
     enddo
   enddo
   rewind(2)
   curres=" ";altloc=" "
   ires=0; 
   do 
     iatom=0
     read(2,'(a)',iostat=ios) aline
     if(ios/=0) exit
     if(aline(1:6)=='ENDMDL')exit  
     if(aline(1:3)=="TER".and.ires/=0)exit                  !! Atoms with the correct chain id but not part of the chain! 
     if(aline(1:6)/='ATOM  '.and.aline(1:6)/='HETATM') cycle
     if(aline(22:22)/=chain.and.aline(22:22)/=" ")cycle
     if(aline(18:26)==curres)cycle
     atype = idatom(aline(13:16),aline(20:22))
     if(atype==-1)cycle                                      !! This is not a backbone atom. Reject.
     ires = ires + 1; iatom = iatom + 1                      !! Attempting to read this residue.
     if(ires>nres)then
        write(*,*)"MISSING BACKBONE ATOMS?(1) ",curres," ",aline(18:26)," ",trim(pdbfile)," ",chain
        stop
     endif
     read(aline(31:54),'(3f8.3))')bbxyz(ires)%xyz(1:3,atype) !! Atype==0 for N.
     curres=aline(18:26)                                     !! Curently reading...
     do   
       altloc=" "
       if(iatom==4)exit                                      !! All four backbone atoms have been read!
       read(2,'(a)',iostat=ios)aline
       if(aline(1:6)/='ATOM  '.and.aline(1:6)/='HETATM') cycle
       qres=aline(18:26)                                     !! Couln't find all four atoms.
       if(curres/=qres)then
         write(*,*)"MISSING BACKBONE ATOMS?(2) ",curres," ",aline(18:26)," ",trim(pdbfile)," ",chain
         stop
       endif
       atype = idatom(aline(13:16),aline(20:22))
       if(atype==-1)cycle
       altloc=aline(17:17)                              
       if(all(bbxyz(ires)%xyz(1:3,atype)==999.))then
         read(aline(31:54),'(3f8.3))')bbxyz(ires)%xyz(1:3,atype) !! Atype==0 for N.
         iatom = iatom + 1
       elseif(altloc/=" ")then
         cycle                                                 !! You already have a copy of this atom.
       else
         write(*,*)"ATTEMPTING TO OVERWRITE ATOM COORDINATES!",ires,atype,trim(pdbfile)," ,",chain
         write(*,*)aline(18:26)
         stop
       endif
     enddo
   enddo
   close(2)
  end subroutine read_backbone
!-------------------------------------------------------------------------------------------------------------
  subroutine crossprod(v1,v2,v3)
   !------------------------------
   implicit none
   !------------------------------
   real,intent(in)::v1(3)
   real,intent(in)::v2(3)
   real,intent(out)::v3(3)
   !------------------------------
   v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
   v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
   v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
  end subroutine crossprod
!-------------------------------------------------------------------------------------------------------------
  integer function idatom(atom,res)
   !----------------------------------
   implicit none
   !----------------------------------
   character(len=4),intent(in)::atom
   character(len=3),intent(in)::res
   !----------------------------------
   integer::i
   character(len=4),dimension(1:4)::atype=(/" N  "," CA "," C  "," O  "/)
   character(len=3),dimension(20)::three=(/'ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE',&
     'LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL',&
     'TRP','TYR'/)
   character(len=1),dimension(21)::res1=(/'A','C','D','E','F','G','H','I','K','L','M','N',&
     'P','Q','R','S','T','V','W','Y','X'/)
   !----------------------------------
   do i=1,4
     if(atom==atype(i))exit
   enddo
   if(i==5)then
     idatom = -1
   else
     idatom=i
   endif
  end function idatom
!-------------------------------------------------------------------------------------------------------------
end module hmmstr


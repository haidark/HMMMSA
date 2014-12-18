module database
  !-------------------------------------
  implicit none
  !-------------------------------------
  public:: database_read_bsf
  !-------------------------------------
  private:: chain_break
  !-------------------------------------
  type frag
    type(struct_node),pointer::istart
    real::dme,score
    type(frag),pointer::next
  end type
  !-------------------------------------
  type gamma_node
    integer::i ! hmmstate
    real::g ! non-zero gamma value
    type(gamma_node),pointer::next
  end type
  !-------------------------------------
  type struct_node
    type(gamma_node),pointer:: gam ! hmm_gamma profile
    integer::n_gamma  !! number of non-zero hmm_gamma values
    character(len=1)::restype,dssp
    real::delta,theta
    real::phi,psi,omg
    real,dimension(3)::calpha
    character(len=5)::code
    type(struct_node),pointer::next
  end type
  !-------------------------------------
 contains
!!!---------------------------------------------------------------------------------------------------------
  subroutine database_read_bsf(bsf,bsf_root,nnodes)
   implicit none
   type(struct_node),pointer,intent(inout)::bsf_root
   type(struct_node),pointer::bsf_ptr
   type(gamma_node),pointer::gam_ptr
   character(len=*),intent(in)::bsf
   integer::ios,ii,i
   integer,optional,intent(inout)::nnodes
   integer::n_gamma
   real::delta,theta,phi,psi,omg
   real,dimension(3)::calpha
   character(len=1)::restype,dssp
   character(len=5)::code

   open(33,file=bsf,form='unformatted',status='old',action='read',iostat=ios)
   if(ios/=0) stop 'Build_DBase :: Error opening database binary structure file!'
   allocate(bsf_root,stat=ios); if (ios/=0) stop 'Build_Dbase :: Error allocating Database root data!'
   bsf_ptr => bsf_root
   ii = 0
   do
     read(33,iostat=ios)n_gamma,code,delta,theta,restype,dssp,phi,psi,omg,calpha
     if (ios/=0) exit
     allocate(bsf_ptr%next,stat=ios) ; if (ios/=0) stop 'Build_DBase :: Error allocating a database node!'
     bsf_ptr => bsf_ptr%next; nullify(bsf_ptr%next)
     bsf_ptr%n_gamma=n_gamma; bsf_ptr%code=code; bsf_ptr%delta=delta; bsf_ptr%theta=theta; bsf_ptr%restype=restype
     bsf_ptr%dssp=dssp; bsf_ptr%phi=phi; bsf_ptr%psi=psi; bsf_ptr%omg=omg; bsf_ptr%calpha=calpha
     allocate(bsf_ptr%gam,stat=ios) ; if (ios/=0) stop 'Build_DBase :: Error allocating gamma node head!'
     gam_ptr => bsf_ptr%gam
     do i=1,bsf_ptr%n_gamma
       read(33) gam_ptr%i, gam_ptr%g
       if (ios/=0) stop 'Build_trgt :: Error missing gamma data!'
       if (i<bsf_ptr%n_gamma) then
         allocate(gam_ptr%next,stat=ios) ; if (ios/=0) stop 'Build_DBase :: Error allocating gamma node!'
         gam_ptr => gam_ptr%next
       endif
     enddo
     ii = ii + 1
   enddo
   bsf_ptr=>bsf_root%next
   deallocate(bsf_root)
   bsf_root=>bsf_ptr
   if(present(nnodes))nnodes=ii
   write(0,*)"Database Nodes ",ii
   close(33)
 end subroutine
!!!---------------------------------------------------------------------------------------------------------
 subroutine database_build_fraglib(bsf_root,frags,nfrag,frag_len)
  implicit none
  type(struct_node),pointer,intent(in)::bsf_root
  integer,intent(in)::frag_len
  integer,intent(out)::nfrag
  type(struct_node),pointer::bsf_ptr
  type(frag),dimension(:),allocatable,intent(out)::frags
  type(frag),pointer::tmplib_root,frag_ptr,frag_ptr2
  logical::broke
  integer::ifrag

   bsf_ptr=>bsf_root
   allocate(tmplib_root)
   frag_ptr=>tmplib_root
   nfrag=0
   !! Build a linked list of pointers to fragment beginnings first
   do
     broke = chain_break(bsf_ptr,frag_len)
     if(.not.broke)then
       allocate(frag_ptr%next)
       frag_ptr=>frag_ptr%next
       frag_ptr%istart=>bsf_ptr
       nfrag=nfrag+1
     endif
     if(associated(bsf_ptr%next))then; bsf_ptr=>bsf_ptr%next; else; exit; endif
   enddo
   !! Remove list head
   frag_ptr=>tmplib_root%next
   deallocate(tmplib_root)
   tmplib_root=>frag_ptr
   write(*,*)"Database Frags ",nfrag
   allocate(frags(nfrag))
   frag_ptr=>tmplib_root
   frag_ptr2=>tmplib_root
   !! Transfer list to an array of fragment pointers
   do ifrag=1,nfrag
      frags(ifrag)=frag_ptr
      frag_ptr=>frag_ptr%next
      deallocate(frag_ptr2)
      frag_ptr2=>frag_ptr
   enddo
  end subroutine
!!!---------------------------------------------------------------------------------------------------------
  function chain_break(bsf_ptr,frag_len)
   implicit none
   logical::chain_break
   type(struct_node),intent(in),pointer::bsf_ptr
   type(struct_node),pointer::bsf_ptr2
   real,dimension(3)::v,c1,c2
   integer,intent(in)::frag_len
   character(len=5)::last_code
   real::dist
   integer::i

   bsf_ptr2=>bsf_ptr
   chain_break=.false.
   c1=bsf_ptr2%calpha
   last_code=bsf_ptr2%code
   do i=1,frag_len-1
     if(associated(bsf_ptr2%next)) then; bsf_ptr2=>bsf_ptr2%next; else; chain_break=.true.; return; endif
     if(last_code/=bsf_ptr2%code) then; chain_break=.true.; return; endif
     c2=bsf_ptr2%calpha
     v = c1(:)-c2(:)
     dist = sqrt(sum(v*v))
     if(dist > 4.0 .or. dist < 2.4 )then
       chain_break = .true.
     endif
     c1=c2
     last_code=bsf_ptr2%code
   enddo
  end function
!!!---------------------------------------------------------------------------------------------------------
  subroutine database_free_bsf(root)
   implicit none
   type(struct_node),intent(inout),pointer::root
   type(struct_node),pointer::db_ptr,db_ptr2
   type(gamma_node),pointer::gam_ptr,gam_ptr2
   integer::i

   db_ptr => root
   do
     db_ptr2 => db_ptr
     gam_ptr => db_ptr%gam
     do i=1,db_ptr%n_gamma
        gam_ptr2 => gam_ptr
        if (i<db_ptr%n_gamma) then
          gam_ptr => gam_ptr%next
        endif
        deallocate(gam_ptr2)
     enddo
     if(associated(db_ptr%next))then; db_ptr => db_ptr%next; else; exit; endif
     deallocate(db_ptr2)
   enddo
   deallocate(db_ptr)
  end subroutine
!!---------------------------------------------------------------------
  subroutine database_dump_xyz(ifrag,xyz,frag_len)
   implicit none
   type(struct_node),intent(in),pointer::ifrag
   type(struct_node),pointer::frag_ptr
   real,dimension(3,frag_len),intent(inout)::xyz
   integer,intent(in)::frag_len
   integer::i

   frag_ptr=>ifrag
   do i=1,frag_len
     xyz(:,i)=frag_ptr%calpha
     if(associated(frag_ptr%next))frag_ptr=>frag_ptr%next;
   enddo
  end subroutine
!!---------------------------------------------------------------------------------------------------------
end module


        program SEQ2DRCT
!!* v.24-sep-99 C.Bystroff byte-size recs
!!* v.17-JUL-98 C.Bystroff
!cb Temporary, I hope, link converting a single sequence
!cb to a drct format file for input to xPSSMpred
!cb Eventually, this should be done within one program...
!cb 30-OCT-96
!!* Convert a readable set of profiles from xmakeVALL
!!* into a single binary direct access file.
!!* Leave space at the end for 8 clusterIDs
!!*
        implicit none
        integer,parameter ::  MAXLINE=1000
        character    :: code*4,ocode*4,chainID*1,ochain*1,sschar*1,aa*1
        character    :: rotamer*1 ! future use
        character    :: domain*1 ! future use
        integer(kind=2)    :: I,iseq
        integer(kind=4)    :: nrec,prec,vrecsize,J,K
        integer(kind=2)    :: resseq, clusterID(3)
        real(kind=4)       :: calpha(3),dih(3),profile(20),csum,Fi(20)
        character(len=80) :: filename
        character(len=MAXLINE)    :: aline
        character(len=23) :: res1
        character(len=1)  :: caps
        logical      :: gapchar,putchar
        integer      :: iargc, jarg, ios=0
!!
        res1='ACDEFGHIKLMNPQRSTVWYBXZ'
!!*
        Fi=(/ 0.08279,0.01937,0.05855,0.05992,0.04014,0.08089, &
        0.02275,0.05552,0.05959,0.08020,0.02070,0.04729,0.04599, &
        0.03728,0.04640,0.06246,0.05888,0.06866,0.01507,0.03756/)  ! overall AA frequencies.
!!* functions
!!*
        vrecsize = 136 ! 136 bytes  machine dependent!!
        ! vrecsize = 34 ! 136 bytes
        jarg = iargc()
        rotamer = ' '
        do I=1,3
          clusterID(I) = 0
        enddo
!!*
        if (jarg.ge.1) then
          call getarg(1,filename)  ! output file
        else
          stop 'Usage: xseq2drct outputfile [label] < sequence'
        endif
        if (jarg.ge.2) then
          call getarg(2,aline)  ! input file
          ocode = aline(1:4)
          ochain = '_'
        else
          ocode = ' '
          ochain = '_'
        endif
!!*
        open(2,file=filename,status='unknown',form='unformatted', &
        access='direct',recl=vrecsize)
        code = filename(1:4)
        if (ocode.ne.' ') code = ocode
        chainID = ochain
        nrec = 0
        prec = 0
        csum = 1.0
        sschar = 'U'
        rotamer = ' '
        do I=1,3
          calpha(I) = 0.0
          dih(I) = 0.0
        enddo
        iseq = 0
        nrec = 0
        write(*,'(''Type or paste sequence (^D when finished): '',$)')
      do
        read(*,'(a)',iostat=ios) aline
        if (ios/=0) exit
        if (aline(1:1).eq.'#') cycle     ! comment
        if (aline(1:1).eq.'>') cycle     ! comment
        I = 1
        do while (I.le.MAXLINE)
          do while (aline(I:I).eq.' '.or.aline(I:I).eq.'.')
            I = I + 1
            if (I.ge.MAXLINE) exit
          enddo
          if (aline(I:I).eq.'#') exit     ! appended comment
          if (aline(I:I).eq.' ') exit     ! appended comment
          iseq = iseq + 1
          aa = caps(aline(I:I))
          !! diagnostic
          write(*,'(a1,1x,a1)') aline(I:I), aa
          J = 1
          do while (aa.ne.res1(J:J))
            J = J + 1
            if (J.gt.23) then
              write(*,'('' Error reading letter! '',a1)') aa
              stop 'seq2drct.f90:: Bad end.2000'
            endif
          enddo
          do K=1,20
            profile(K) = 0.0
          enddo
          if (J.gt.20) then
            if (J.eq.21) then ! B = N or D
              profile(3) = 0.50*Fi(3)
              profile(12) = 0.50*Fi(12)
            elseif (J.eq.22) then ! X = anything
              do K=1,20
                profile(K) = Fi(K)
              enddo
            elseif (J.eq.23) then ! Z = Q or E
              profile(4) = 0.50*Fi(4)
              profile(14) = 0.50*Fi(14)
            else ! error
              write(*,'('' Error reading letter! '',a1)') aa
              stop 'seq2drct.f90:: Bad end.2001'
            endif
          else
            profile(J) = 1.0
          endif
          nrec = nrec + 1
!cb Updated to new vall format 17-OCT-00 
          write(2,rec=nrec) iseq,iseq,code,chainID,aa,sschar,rotamer, &
          calpha,dih,csum,profile,0.,0.,'_','?',clusterID
          I = I + 1
!cb        putchar = (gapchar.and..not.putchar)  
!cb        if (putchar) then
!cb          nrec = nrec + 1
!cb          write(2,rec=nrec) iseq,iseq,code,chainID,'-',sschar,rotamer,
!cb   &      calpha,dih,csum,profile,clusterID
!cb        endif
        enddo
      enddo
        write(*,'('' Sequence length: '',i10)') iseq
        write(*,'('' Last record: '',i10)') nrec
        prec = 4*nrec
        write(*,'('' Total size: '',i10,'' bytes.'')') prec
        close(2)
        stop 'SEQ2DRCT: Normal end'
        end program seq2drct
!!**********************************************************************
        character(len=1) function caps(card)
        character(len=1) :: card
        logical :: gapchar
        integer :: I,icap
        icap = ichar('A') - ichar('a')
        ! do I=1,len(card)
        I = 1
          gapchar = ('a'.le.card.and.card.le.'z')
          if (gapchar) card = char(ichar(card)+icap)
        ! enddo
        caps = card
        return
        end function caps
!!**********************************************************************


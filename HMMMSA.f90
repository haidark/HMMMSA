!-------------------------------------------------------------------------------
!Heading:           MSA based on the HMMSUM-D model
!Date:              2014/12/8
!Version:           0.0.2
!Authors:           Yao-ming Huang & Haidar Khan
!Purpose:           A multiple sequence alignment program based on HMMSUM-D model
!Usage:
!File Format:       Input: a FASTA file with multi-sequences(12), observed counts, expected counts(13)
!                   Temp: tempseqfile(14)
!                   Output: MSA file(17)
!Note:              # Smith-Waterman algorithm, Local alignment, Affine penalty, No end gap penalty
!                   Use the LLR of weighted probabilities obtained from frequencies directly
!                   Use pseudocount=0.000000001
!                   Use half-bits score
!                   Use 281 states in gamma values
!                   # Bayesian adaptive sequence alignment algorithm (matrix bias is NOT a posterior)
!                   Use HMMSUM-D model
!                   Set bias to 0
!Reference:
!Restriction:
!Revision History:
!
!Copyright (C) 2007 Yao-ming Huang.
!All Rights Reserved.
!Permission to use, copy, modify, and distribute this software and its documentation for non-profit purposes
!and without fee is hereby granted provided that this copyright notice appears in all copies.
!-------------------------------------------------------------------------------

program HMMMSA
  
  implicit none
  integer::memstatus=0, openstatus=0, readstatus=0, sysstatus=0, system, argnum, totalseqnum, totalpair, seqnum, i, j, &
           x, y, w, seqI, seqJ, a, b, p1, p2, c, ss1, ss2, k, m, n, block, L, seqL1, seqL2, len1, len2, child, &
           DPscore, Mnode, aa, counts
  integer(KIND=1)::iargc
  real::e, etime, time(2), bias, pseudo, shift, gapO, gapE, iden, sumM, sumD, sumL, sumE
  character(1)::signal
  character(27)::alphabet="ACDEFGHIKLMNPQRSTVWYBZX.-~*"
  character(1000)::temp1, temp2, aline, infile, outfile, alifile, obsfile, expfile, tempseqfile, tempdrctfile, &
                   tempgcafile, num, string, seqname1, seqname2, parentdrct, parentgca, childdrct, childgca
  character(9999)::seq, seq1, seq2
  logical::alive
  integer, dimension(2)::pos
  integer, dimension(13)::sarray
  integer, allocatable::tag(:)
  real, allocatable::matrix1(:,:,:), matrix2(:,:), gammaM1(:,:), gammaM2(:,:), alignM(:,:), dist(:,:), aliscore(:,:), &
                     pro1(:,:), pro2(:,:), HMM(:,:,:), emit(:,:)
  character(1), allocatable::tback(:), tbback(:)
  character(13), allocatable::seqID(:)
  character(1000), allocatable::nodetable(:)
  character(9999), allocatable::seqchar(:), seqpair(:)

  !-------------------------------------------------------------------------------
  !open and generate drct/gca files from FASTA file
  !-------------------------------------------------------------------------------
  e=etime(time)
  argnum=iargc()
  call unlink("LOG")
  if (argnum==0 .or. mod(argnum, 2)/=0) then
    write (*, *) "xHMMMSA2   V. 0.0.2 Dec 2014"
    write (*, *) "A MSA program based on the HMMSUM-D model"
    write (*, *) "Usage:"
    write (*, *) "  xHMMMSA INPUT OUTPUT -s {Shift} -b {Bias}"
    write (*, *) "  -> INPUT: a FASTA file with multi-sequences"
    write (*, *) "  -> OUTPUT: a file name for the output file"
    write (*, *) "  -> Shift: Matrix bias in DP, default -> 0"
    write (*, *) "  -> Bias: Matrix bias in BAA, default -> 0"
    stop
  end if
  bias=0
  shift=0
  gapO=21
  gapE=0.4
  call getarg(1, temp1)
  call getarg(2, temp2)
  infile=trim(temp1)
  outfile=trim(temp2)
  do i=3, argnum, 2
    call getarg(i, temp1)
    call getarg(i+1, temp2)
    if (temp1=="-s") then
      read (temp2, *, iostat=readstatus) shift
      if (readstatus/=0) stop "Wrong argument in shift!!!"
    elseif (temp1=="-b") then
      read (temp2, *, iostat=readstatus) bias
      if (readstatus/=0) stop "Wrong argument in bias!!!"
    else
      stop "Wrong argument(s)!!!"
    end if
  end do
  open (unit=17, file=outfile, status="replace", iostat=openstatus)

  !-----open FASTA file and count sequences-----!
  open (unit=12, file=infile, status="old", iostat=openstatus)
  if (openstatus/=0) then
    write (*, "(A, $)") "ERROR: File not found -"
    write (*, *) trim(infile)
    stop
  end if
  totalseqnum=0
  do
    read (12, "(A)", iostat=readstatus) temp1
    if (readstatus/=0) exit
    if (temp1(1:1)==">") totalseqnum=totalseqnum+1
  end do
  write (*, *) 
  write (*, *) "Sequence number in the FASTA file:", totalseqnum
  write (*, *)
  allocate (seqID(totalseqnum), stat=memstatus)
  if (memstatus/=0) stop 'ERROR: Allocating seqID!!!'
  seqID=" "
  rewind(12)
  i=1
  do
    read (12, "(A)", iostat=readstatus) temp1
    if (readstatus/=0) exit
    if (temp1(1:1)==">") then
      seqID(i)=temp1(2:11)
      i=i+1
    end if
  end do
  rewind(12)
  read (12, "(A)", iostat=readstatus) temp1
  allocate (seqchar(totalseqnum), stat=memstatus)
  if (memstatus/=0) stop 'ERROR: Allocating seqchar!!!'
  seqchar=" "
  seqnum=0
  seq=" "
  do while (readstatus==0)
    if (temp1(1:1)==">") seqnum=seqnum+1
    write (num, *) seqnum
    tempseqfile=trim(adjustL(num)) // '.seq'
    tempdrctfile=trim(adjustL(num)) // '.drct'
    tempgcafile=trim(adjustL(num)) // '.gca'
    if (stat(tempseqfile, sarray)==0) call unlink(tempseqfile)
    if (stat(tempdrctfile, sarray)==0) call unlink(tempdrctfile)
    if (stat(tempgcafile, sarray)==0) call unlink(tempgcafile)
    write (*, *) "...parsing sequence ", trim(adjustL(num))

    !-----generate sequence file-----!
    open (unit=14, file=tempseqfile, status="replace", iostat=openstatus)
    !write (14, "(A1, $)") ">"
    close(14)
    string=' cat ' // trim(infile) // ' | grep ">" | head -' // trim(adjustL(num)) // &
           ' | tail -1 | awk ''{print $2}'' | sed -e "s/(//" -e "s/)//" >> ' // trim(tempseqfile)
    !write (*, *) "Running:", trim(string)
    sysstatus=system(string)
    if (sysstatus /= 0) stop "ERROR: System command (sequence file)!!!"
    read (12, "(A)", iostat=readstatus) temp1
    do while (temp1(1:1)/=">")
      seq=trim(seq) // trim(temp1)
      read (12, "(A)", iostat=readstatus) temp1
      if (readstatus/=0) exit
    end do
    seqchar(seqnum)=trim(seq)
    open (unit=14, file=tempseqfile, status="old", iostat=openstatus)
    read (14, *)
    write (14, "(A)") ">fillertext"
    write (14, "(A)") trim(seq)
    seq=" "
    close(14)
    !-----generate drct file-----!
    call seq2drct(tempseqfile)
    
    !-----generate gca file-----!
    call drct2gca(tempdrctfile)
    
    call unlink(tempseqfile)
  end do
  close(12)
  write (*, *)
  
  !-------------------------------------------------------------------------------
  !perform pairwise alignments based on HMMSUM-D
  !-------------------------------------------------------------------------------
  !Load observed and expected frequencies from HMMSTR into matrix1 & matrix2
  obsfile="observed"
  expfile="expected"
  allocate (matrix1(20, 20, 281), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating matrix1 for observed frequencies!!!"
  allocate (matrix2(20, 281), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating matrix2 for expected frequencies!!!"
  matrix1=0
  matrix2=0
  call readmatrix(obsfile, matrix1, 20, 20, 281)
  open (unit=13, file=expfile, status="old", iostat=openstatus)
  if (openstatus/=0) then
    write (*, "(A, $)") "ERROR: File not found -"
    write (*, *) trim(expfile)
    stop
  end if
  w=1
  do
    read (13, "(A)", iostat=readstatus) aline
    if (readstatus/=0) exit
    if (aline==" ") cycle
    read (aline, "(20F15.7)", iostat=readstatus) (matrix2(x, w), x=1, 20)
    if (readstatus/=0) stop "ERROR: Bad matrix2!!!"
    w=w+1
    if (w>281) exit
  end do
  close(13)

  !-----add pseudocounts-----!
  pseudo=0.000000001
  matrix1=matrix1+(pseudo**2) !observed counts + pseudocounts**2
  matrix2=matrix2+pseudo      !background counts + pseudocounts

  !-----allocate distance matrix-----!
  allocate (dist(totalseqnum, totalseqnum), stat=memstatus)
  if (memstatus/=0) stop 'ERROR: Allocating dist!!!'      
  dist=-1

  !-----allocate alignment score matrix-----!
  allocate (aliscore(totalseqnum, totalseqnum), stat=memstatus)
  if (memstatus/=0) stop 'ERROR: Allocating aliscore!!!'
  aliscore=0
  DPscore=0

  !-----read gamma values into memory-----!
  totalpair=(totalseqnum*(totalseqnum-1))/2
  counts=0
  write (*, *) "...performing all possible pairwise alignments based on HMMSUM-D(DP)"
  do seqI=1, totalseqnum
    do seqJ=seqI+1, totalseqnum
      counts=counts+1
      !progress bar code
      write (*, "(A, $)") "\r"
      write (*, "(I3, A2, $)") floor((real(counts)/real(totalpair))*100.0), "% "
      write (*, "(A1, $)") "|"
      do i=1, floor((real(counts)/real(totalpair))*100.0)/2
        write (*, "(A1, $)") "="
      end do
      if (floor((real(counts)/real(totalpair))*100.0)==100) then
        write (*, "(A1, $)") "|"
        write (*, *)
      else
        write (*, "(A1, $)") ">"
      end if
      !end progress bar code
      
      !read gamma values from .gca files
      write (num, *) seqI
      seqname1=trim(adjustL(num)) // '.gca'
      write (num, *) seqJ
      seqname2=trim(adjustL(num)) // '.gca'
      seq1=seqchar(seqI)
      seq2=seqchar(seqJ)
      a=len_trim(seq1)
      b=len_trim(seq2)
      allocate (gammaM1(a, 0:281), stat=memstatus)
      if (memstatus/=0) stop "ERROR: Allocating gammaM1!!!"
      allocate (gammaM2(b, 0:281), stat=memstatus)
      if (memstatus/=0) stop "ERROR: Allocating gammaM2!!!"
      gammaM1=0
      gammaM2=0
      call readgamma(seqname1, a, gammaM1)
      call readgamma(seqname2, b, gammaM2)

      !-----generate Mij-----!
      allocate (alignM(a, b), stat=memstatus)
      if (memstatus/=0) stop "ERROR: Allocating align matrix!!!"
      alignM=0
      !returns alignment Matrix between the two sequences
      call makeMij(a, b, seq1, seq2, gammaM1, gammaM2, matrix1, matrix2, alignM)
      alignM=2*alignM+shift
      !maxS=Maxval(alignM)
      !minS=Minval(alignM)
      !averS=sum(alignM)/(a*b)

      !-----run Dynamic programming-----!
      c=a+b
      allocate (tback(c), stat=memstatus)
      if (memstatus/=0) stop "ERROR: Allocating tback!!!"
      tback=" "
      allocate (tbback(c), stat=memstatus)
      if (memstatus/=0) stop "ERROR: Allocating tbback!!!"
      tbback=" "
      call DP(a, b, gapO, gapE, alignM, tback, DPscore)
      aliscore(seqI, seqJ)=DPscore

      !-----write pairwise alignments to outfile(17)-----!
      ss1=a-1   !In order to match the position of sequences
      ss2=b-1
      do i=1, c
        tbback(i)=tback(c+1-i)
      end do
      tback=" "
      k=1
      do i=1, c
        if (tbback(i)==" ") cycle
        if (tbback(i)=="E") cycle
        tback(k)=tbback(i)
        k=k+1
      end do
      a=len_trim(seq1)
      b=len_trim(seq2)
      block=ceiling(real(k-2)/50)
      m=2
      n=2
      allocate (seqpair(2), stat=memstatus)
      if (memstatus/=0) stop 'ERROR: Allocating seqpair!!!'
      seqpair=" "
      k=1
      L=1
      write (17, "(A, $)") "Initial pairwise alignment: "
      write (17, "(A, $)") trim(seqname1)
      write (17, "(A, $)") ", "
      write (17, "(A, $)") trim(seqname2)
      write (17, *)
      do i=1, block
        write (17, "(A, $)") "QUERY        "
        write (17, "(I4, $)") ss1
        write (17, "(A, $)") "   "
        do j=1, 50
          if (m>c) cycle
          if (tback(m)=="M" .OR. tback(m)=="I") then
            write (17, "(A1, $)") seq1(ss1:ss1)
            seqpair(1)(k:k)=seq1(ss1:ss1)
            ss1=ss1+1
            k=k+1
          elseif (tback(m)=="D") then
            write (17, "(A1, $)") "."
            seqpair(1)(k:k)="."
            k=k+1
          end if
          m=m+1
        end do
        write (17, "(A, $)") "   "
        write (17, "(I4, $)") ss1-1
        write (17, *)
        write (17, "(A, $)") "TARGET       "
        write (17, "(I4, $)") ss2
        write (17, "(A, $)") "   "
        do j=1, 50
          if (m>c) cycle
          if (tback(n)=="M" .OR. tback(n)=="D") then
            write (17, "(A1, $)") seq2(ss2:ss2)
            seqpair(2)(L:L)=seq2(ss2:ss2)
            ss2=ss2+1
            L=L+1
          elseif (tback(n)=="I") then
            write (17, "(A1, $)") "."
            seqpair(2)(L:L)="."
            L=L+1
          end if
          n=n+1
        end do
        write (17, "(A, $)") "   "
        write (17, "(I4, $)") ss2-1
        write (17, *)
        write (17, *)
      end do
      write (17, *)
      call flush(17)
   
      !-----percent identity-----!
      ! iden = # of matches between the two sequences (exclude gaps)
      ! seqL* = total # of residues in each sequence
      ! iden = iden/ max(seqL*)
      iden=0
      seqL1=0
      seqL2=0
      do i=1, len_trim(seqpair(1))
        if (index(alphabet, seqpair(1)(i:i))<=20) seqL1=seqL1+1
        if (index(alphabet, seqpair(2)(i:i))<=20) seqL2=seqL2+1
        if (index(alphabet, seqpair(1)(i:i))<=20 .AND. index(alphabet, seqpair(2)(i:i))<=20) then
          if (seqpair(1)(i:i)==seqpair(2)(i:i)) then
            iden=iden+1
          end if
        end if
      end do
      if (seqL1>=seqL2) then
        iden=iden/seqL2
      else
        iden=iden/seqL1
      end if
      !percent identity becomes distance between the two sequences
      dist(seqI, seqJ)=iden
      ! deallocate resources used in this loop
      deallocate(tback)
      deallocate(tbback)
      deallocate(alignM)
      deallocate(seqpair)
      deallocate(gammaM1)
      deallocate(gammaM2)
    end do
  end do
  deallocate(seqchar)
  write (*, *)

  !-------------------------------------------------------------------------------
  !generate the parent profile using BAA
  !-------------------------------------------------------------------------------

  !-----define order of pairs----!
  write (*, *) "...generating the parent profile using HMMSUM-D(BAA)"
  allocate (nodetable((totalseqnum-1)*totalseqnum/2), stat=memstatus)
  if (memstatus/=0) stop 'ERROR: Allocating nodetable!!!'
  ! order the pairs of sequences by similarity
  nodetable=" "
  j=1
  do
    temp1=" "
    temp2=" "
    if (maxval(dist)==-1) exit
    pos=maxloc(dist)
    p1=pos(1)
    p2=pos(2)
    dist(p1, p2)=-1
    do i=1, totalseqnum
      if (i/=p1 .and. i/=p2) then
        if (i>p1 .and. i>p2) then
          if (dist(p1, i)>=dist(p2, i)) then
            dist(p1, i)=(dist(p1, i)+dist(p2, i))/2
            dist(p2, i)=-1
          else
            dist(p2, i)=(dist(p1, i)+dist(p2, i))/2
            dist(p1, i)=-1
          end if
        elseif (i>p1 .and. i<p2) then
          if (dist(p1, i)>=dist(i, p2)) then
            dist(p1, i)=(dist(p1, i)+dist(i, p2))/2
            dist(i, p2)=-1
          else
            dist(i, p2)=(dist(p1, i)+dist(i, p2))/2
            dist(p1, i)=-1
          end if
        elseif (i<p1 .and. i>p2) then
          if (dist(i, p1)>=dist(p2, i)) then
            dist(i, p1)=(dist(i, p1)+dist(p2, i))/2
            dist(p2, i)=-1
          else
            dist(p2, i)=(dist(i, p1)+dist(p2, i))/2
            dist(i, p1)=-1
          end if
        else
          if (dist(i, p1)>=dist(i, p2)) then
            dist(i, p1)=(dist(i, p1)+dist(i, p2))/2
            dist(i, p2)=-1
          else
            dist(i, p2)=(dist(i, p1)+dist(i, p2))/2
            dist(i, p1)=-1
          end if
        end if
      end if
    end do

    do i=1, j
      write (num, *) p1
      if (index(nodetable(i), trim(adjustL(num)))<=0) then
        aline=trim(adjustL(num))
        if (len_trim(aline)>=len_trim(temp1)) temp1=trim(aline)
      elseif (index(nodetable(i), trim(adjustL(num)))>0) then
        aline=trim(nodetable(i))
        if (len_trim(aline)>=len_trim(temp1)) temp1=aline
      end if
      write (num, *) p2
      if (index(nodetable(i), trim(adjustL(num)))<=0) then
        aline=trim(adjustL(num))
        if (len_trim(aline)>=len_trim(temp2)) temp2=trim(aline)
      elseif (index(nodetable(i), trim(adjustL(num)))>0) then
        aline=trim(nodetable(i))
        if (len_trim(aline)>=len_trim(temp2)) temp2=aline
      end if
    end do
    if (trim(temp1)==trim(temp2)) exit
    nodetable(j)=trim(temp1) // "-" // trim(temp2)
    j=j+1
    tempgcafile=trim(temp1) // '.gca'
    seqname1=trim(temp1) // '.drct'
    inquire(file=tempgcafile, exist=alive)
    if (alive) then
      write (*, *) trim(tempgcafile), " existed"
    else
      write (*, *) "...making ", trim(tempgcafile)
      call drct2gca(seqname1)
    end if
    tempgcafile=trim(temp2) // '.gca'
    seqname2=trim(temp2) // '.drct'
    inquire(file=tempgcafile, exist=alive)
    if (alive) then
      write (*, *) trim(tempgcafile), " existed"
    else
      write (*, *) "...making ", trim(tempgcafile)
      call drct2gca(seqname2)
    end if
    parentdrct=trim(temp1) // "-" // trim(temp2) // '.drct'
    parentgca=trim(temp1) // "-" // trim(temp2) // '.gca'
    call getdrctlen(seqname1, len1)
    call getdrctlen(seqname2, len2)
    
    call BAA(seqname1, seqname2, len1, len2, matrix1, matrix2, bias)
  end do
  deallocate(dist)
  deallocate(nodetable)
  
  !-------------------------------------------------------------------------------
  !generate MSA using DP (parent profile V.S. each sequence)
  !-------------------------------------------------------------------------------

  !-----generate gca file from parent drct file-----!
  write (*, *) "...making initial multiple sequence alignment"
  write (*, *)
  !allocate SEQPAIR here
  allocate (seqpair(2*totalseqnum), stat=memstatus)
  if (memstatus/=0) stop 'ERROR: Allocating seqpair!!!'
  seqpair=" "
  allocate (tag(totalseqnum), stat=memstatus)
  if (memstatus/=0) stop 'ERROR: Allocating tag!!!'
  tag=0
  inquire(file=parentgca, exist=alive)
  if (alive) then
    write (*, *) trim(parentgca), " existed"
    write (*, *)
  else
    write (*, *) "...making ", trim(parentgca)
    call drct2gca(parentdrct)
  end if
  call getdrctlen(parentdrct, len1)
  allocate (pro1(len1, 20), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating pro1!!!"
  pro1=0
  call readdrct(parentdrct, len1, pro1)
  allocate (gammaM1(len1, 0:281), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating gammaM1!!!"
  gammaM1=0
  call readgamma(parentgca, len1, gammaM1)
  call readgcaseq(parentgca, seq1)
  !DEBUG
  !write(*,*) trim(parentgca) // ": " // trim(seq1)
  counts=0
  ! compare each child to parent profile generated from BAA
  do child=1, totalseqnum     
    !progress bar code
    counts=counts+1
    write (*, "(A, $)") "\r"
    write (*, "(I3, A2, $)") floor((real(counts)/real(totalseqnum))*100.0), "% "
    write (*, "(A1, $)") "|"
    do i=1, floor((real(counts)/real(totalseqnum))*100.0)/2
      write (*, "(A1, $)") "="
    end do
    if (floor((real(counts)/real(totalseqnum))*100.0)==100) then
      write (*, "(A1, $)") "|"
      write (*, *)
    else
      write (*, "(A1, $)") ">"
    end if
    !end progress bar code
    
    ! get child drctfile and gca file and extract data from it
    write (num, *) child
    childdrct=trim(adjustL(num)) // '.drct'
    childgca=trim(adjustL(num)) // '.gca'
    call getdrctlen(childdrct, len2)
    
    allocate (pro2(len2, 20), stat=memstatus)
    if (memstatus/=0) stop "ERROR: Allocating pro2!!!"
    pro2=0
    call readdrct(childdrct, len2, pro2)
    allocate (gammaM2(len2, 0:281), stat=memstatus)
    if (memstatus/=0) stop "ERROR: Allocating gammaM2!!!"
    gammaM2=0
    call readgamma(childgca, len2, gammaM2)
    call readgcaseq(childgca, seq2)
    !DEBUG
    !write(*,*) trim(childgca) // ": " // trim(seq2)
    !-----generate Mij-----!
    allocate (alignM(len1, len2), stat=memstatus)
    if (memstatus/=0) stop "ERROR: Allocating align matrix!!!"
    alignM=0
    call makedrctMij(len1, len2, pro1, pro2, gammaM1, gammaM2, matrix1, matrix2, alignM)
    alignM=2*alignM+shift

    !-----run Dynamic programming-----!
    c=len1+len2

    allocate (tback(c), stat=memstatus)
    if (memstatus/=0) stop "ERROR: Allocating tback!!!"
    tback=" "
    allocate (tbback(c), stat=memstatus)
    if (memstatus/=0) stop "ERROR: Allocating tbback!!!"
    tbback=" "
    call DP(len1, len2, gapO, gapE, alignM, tback, DPscore)

   
    !-----output pairwise alignments-----!
    ss1=len1-1   !In order to match the position of sequences
    ss2=len2-1
    do i=1, c
      tbback(i)=tback(c+1-i)
    end do
    tback=" "
    k=1
    do i=1, c
      if (tbback(i)==" ") cycle
      if (tbback(i)=="E") cycle
      tback(k)=tbback(i)
      k=k+1
    end do
    !DEBUG
    !write(*,*) tback
    call getdrctlen(parentdrct, len1)
    call getdrctlen(childdrct, len2)
    block=ceiling(real(k-2)/50)
    m=2
    n=2
    k=1
    L=1
    write (17, "(A, $)") "Multiple sequence alignment (Iteration 1): "
    write (17, "(A, $)") trim(parentdrct)
    write (17, "(A, $)") ", "
    write (17, "(A, $)") trim(childdrct)
    write (17, *)
    tag(child)=ss2
    do i=1, block
      write (17, "(A, $)") "PARENT       "
      write (17, "(I4, $)") ss1
      write (17, "(A, $)") "   "
      do j=1, 50
        if (tback(m)=="M" .OR. tback(m)=="I") then
          write (17, "(A1, $)") seq1(ss1:ss1)
          seqpair(2*child-1)(k:k)=seq1(ss1:ss1)
          ss1=ss1+1
          k=k+1
        elseif (tback(m)=="D") then
          write (17, "(A1, $)") "."
          seqpair(2*child-1)(k:k)="."
          k=k+1
        end if
        m=m+1
      end do
      write (17, "(A, $)") "   "
      write (17, "(I4, $)") ss1-1
      write (17, *)
      write (17, "(A, $)") "CHILD        "
      write (17, "(I4, $)") ss2
      write (17, "(A, $)") "   "
      do j=1, 50
        if (tback(n)=="M" .OR. tback(n)=="D") then
          write (17, "(A1, $)") seq2(ss2:ss2)
          !DEBUG
          !write (*,"(A, $)") seq2(ss2:ss2)
          seqpair(2*child)(L:L)=seq2(ss2:ss2)
          ss2=ss2+1
          L=L+1
        elseif (tback(n)=="I") then
          write (17, "(A1, $)") "."
          !DEBUG
          !write (*, "(A, $)") "."
          seqpair(2*child)(L:L)="."
          L=L+1
        end if
        n=n+1
      end do
      write (17, "(A, $)") "   "
      write (17, "(I4, $)") ss2-1
      write (17, *)
      write (17, *)
    end do
    write (17, *)
    call flush(17)
    deallocate(pro2)
    deallocate(gammaM2)
    deallocate(alignM)
    deallocate(tback)
    deallocate(tbback)
  end do
  deallocate(matrix1)
  deallocate(matrix2)
  deallocate(pro1)
  deallocate(gammaM1)
  write (*, *)

  !-----output multiple sequence alignments to outfile (17)-----!
  call makeMSA(totalseqnum, seqpair)
  !do i=1, totalseqnum
  !  write (*, *) trim(seqpair(2*i))
  !end do
  c=len_trim(seqpair(1))
  write (*, *) c
  block=ceiling(real(c)/50)
  !write (*, *) block
  write (17, "(A, $)") "//"
  write (17, *)
  write (17, *)
  write (17, *)
  !write (17, "(A, $)") "Multiple sequence alignment"
  write (17, *)
  do i=1, block
    do j=1, totalseqnum
      !Write the sequence ID
      write (17, "(A15, $)") adjustL(seqID(j))
      !write the start number (remove?)
      !write (17, "(I4, $)") tag(j)
      !write the sequence (need blocks of 10)
      write (17, "(A, $)") " "
      do k = 0, 4
         write (17, "(A, $)") seqpair(2*j)((50*(i-1)+1)+10*k:(50*(i-1)+1)+10*(k+1)-1)
         write (17, "(A, $)") " "
      end do
      !write (17, "(A, $)") seqpair(2*j)((50*(i-1)+1):50*i)
      !write (17, "(A, $)") " "
      !count the number of residues
      !do k=1, 50
      !  if (index(alphabet, seqpair(2*j)((50*(i-1)+k):(50*(i-1)+k)))>=1 .and. &
      !      index(alphabet, seqpair(2*j)((50*(i-1)+k):(50*(i-1)+k)))<=20) tag(j)=tag(j)+1
      !end do
      !write the end number
      !write (17, "(I4, $)") tag(j)-1
      !advance to next line
      write (17, *)
    end do
    write (17, *)
    write (17, *)
  end do
  deallocate(tag)
  call flush(17)
  
  !-------------------------------------------------------------------------------
  !generate the profile HMM from initial MSA
  !-------------------------------------------------------------------------------

  !-----pick parent sequence from MSA-----!
  x=0
  y=0
  do i=1, totalseqnum
    if (sum(aliscore(i, :))+sum(aliscore(:, i))>=x) then
      x=sum(aliscore(i, :))+sum(aliscore(:, i))
      y=i
    end if
  end do
  write (*, *) trim(seqID(y)), " has been selected as the parent sequence"
  write (*, *)
  deallocate(aliscore)

  !-----define the number of match states-----!
  i=1
  Mnode=0
  do while (index(alphabet, seqpair(2*y)(i:i))>=1)
    if (index(alphabet, seqpair(2*y)(i:i))<=20) Mnode=Mnode+1
    i=i+1
  end do
  write (*, *) "...define ", Mnode, " states"
  write (*, *)

  !-----build a HMM profile-----!
   allocate (HMM(3, 3, Mnode+1), stat=memstatus)
  if (memstatus/=0) stop 'ERROR: Allocating HMM!!!'
  HMM=0
  signal=" "
  do i=2, 2*totalseqnum, 2
    if (i==2*y) cycle     ! don't count parent sequence

    !-----beginning run-----!
    j=1
    do
      if (index(alphabet, seqpair(i)(j:j))<=20 .and. index(alphabet, seqpair(2*y)(j:j))<=20) then        ! M state
        HMM(1, 1, 1)=HMM(1, 1, 1)+1
        j=j+1
        signal="M"
        exit
      elseif (index(alphabet, seqpair(i)(j:j))>20 .and. index(alphabet, seqpair(2*y)(j:j))<=20) then     ! D state
        HMM(2, 1, 1)=HMM(2, 1, 1)+1
        j=j+1
        signal="D"
        exit
      else
        j=j+1
      end if
    end do
    k=1

    !-----continuous run-----!
    do while (index(alphabet, seqpair(i)(j:j))>=1)
      !write (17, *) i, j, k, seqpair(i)(j:j), " ",  seqpair(2*y)(j:j), " ",  signal
      if (index(alphabet, seqpair(i)(j:j))<=20 .and. index(alphabet, seqpair(2*y)(j:j))<=20) then      ! M state
        if (signal=="M") then         ! M -> M
          k=k+1
          HMM(1, 1, k)=HMM(1, 1, k)+1
          j=j+1
          signal="M"
        elseif (signal=="D") then     ! D -> M
          k=k+1
          HMM(1, 2, k)=HMM(1, 2, k)+1
          j=j+1
          signal="M"
        elseif (signal=="L") then     ! L -> M
          k=k+1
          HMM(1, 3, k)=HMM(1, 3, k)+1
          j=j+1
          signal="M"
        end if
      elseif (index(alphabet, seqpair(i)(j:j))>20 .and. index(alphabet, seqpair(2*y)(j:j))<=20) then   ! D state
        if (signal=="M") then         ! M -> D
          k=k+1
          HMM(2, 1, k)=HMM(2, 1, k)+1
          j=j+1
          signal="D"
        elseif (signal=="D") then     ! D -> D
          k=k+1
          HMM(2, 2, k)=HMM(2, 2, k)+1
          j=j+1
          signal="D"
        elseif (signal=="L") then     ! L -> D
          k=k+1
          HMM(2, 3, k)=HMM(2, 3, k)+1
          j=j+1
          signal="D"
        end if
      elseif (index(alphabet, seqpair(i)(j:j))<=20 .and. index(alphabet, seqpair(2*y)(j:j))>20) then   ! L state
        if (signal=="M") then         ! M -> L
          k=k
          HMM(3, 1, k)=HMM(3, 1, k)+1
          j=j+1
          signal="L"
        elseif (signal=="D") then     ! D -> L
          k=k
          HMM(3, 2, k)=HMM(3, 2, k)+1
          j=j+1                       
        elseif (signal=="L") then     ! L -> L
          k=k
          HMM(3, 3, k)=HMM(3, 3, k)+1
          j=j+1                       
        end if 
      elseif (index(alphabet, seqpair(i)(j:j))>20 .and. index(alphabet, seqpair(2*y)(j:j))>20) then
        j=j+1
      end if
      if(k == Mnode) exit
    end do
    
    !-----ending run-----!
    if (signal=="M") then
      k=k+1
      HMM(1, 1, k)=HMM(1, 1, k)+1
      signal="E"
    elseif (signal=="D") then
      k=k+1
      HMM(2, 1, k)=HMM(2, 1, k)+1
      signal="E"
    elseif (signal=="L") then
      k=k+1
      HMM(3, 1, k)=HMM(3, 1, k)+1
      signal="E"
    end if
  end do

  !-----convert to transition probability-----!
  k=1
  sumM=HMM(1, 1, k)+HMM(2, 1, k)
  if (sumM>0) then
    HMM(1, 1, k)=HMM(1, 1, k)/sumM
    HMM(2, 1, k)=HMM(2, 1, k)/sumM
  end if
  sumD=HMM(1, 2, k)+HMM(2, 2, k)
  if (sumD>0) then
    HMM(1, 2, k)=HMM(1, 2, k)/sumD
    HMM(2, 2, k)=HMM(2, 2, k)/sumD
  end if
  sumL=HMM(1, 3, k)+HMM(2, 3, k)
  if (sumL>0) then
    HMM(1, 3, k)=HMM(1, 3, k)/sumL
    HMM(2, 3, k)=HMM(2, 3, k)/sumL
  end if
  do k=2, Mnode+1
    sumM=0
    sumD=0
    sumL=0
    sumM=HMM(1, 1, k)+HMM(2, 1, k)+HMM(3, 1, k-1)
    if (sumM>0) then
      HMM(1, 1, k)=HMM(1, 1, k)/sumM
      HMM(2, 1, k)=HMM(2, 1, k)/sumM
      HMM(3, 1, k-1)=HMM(3, 1, k-1)/sumM
    end if
    sumD=HMM(1, 2, k)+HMM(2, 2, k)+HMM(3, 2, k-1)
    if (sumD>0) then
      HMM(1, 2, k)=HMM(1, 2, k)/sumD
      HMM(2, 2, k)=HMM(2, 2, k)/sumD
      HMM(3, 2, k-1)=HMM(3, 2, k-1)/sumD
    end if
    sumL=HMM(1, 3, k)+HMM(2, 3, k)+HMM(3, 3, k-1)
    if (sumL>0) then
      HMM(1, 3, k)=HMM(1, 3, k)/sumL
      HMM(2, 3, k)=HMM(2, 3, k)/sumL
      HMM(3, 3, k-1)=HMM(3, 3, k-1)/sumL
    end if
  end do
 
  !-----output transition probabilities-----!
  !write (17, *) "M state:"
  !do i=1, 3
  !  do j=1, Mnode+1
  !    write (17, "(F5.2, $)") HMM(1, i, j)
  !  end do
  !  write (17, *)
  !end do
  !write (17, *)
  !write (17, *) "D state:"
  !do i=1, 3
  !  do j=1, Mnode+1
  !    write (17, "(F5.2, $)") HMM(2, i, j)
  !  end do
  !  write (17, *)
  !end do
  !write (17, *)
  !write (17, *) "L state:"
  !do i=1, 3
  !  do j=1, Mnode+1
  !    write (17, "(F5.2, $)") HMM(3, i, j)
  !  end do
  !  write (17, *)
  !end do

  !-----calculate emission probabilities-----!
  allocate (emit(Mnode, 20), stat=memstatus)
  if (memstatus/=0) stop 'ERROR: Allocating emit!!!'
  emit=0
  i=1
  k=1
  do while (index(alphabet, seqpair(2*y)(i:i))>=1)
    if (index(alphabet, seqpair(2*y)(i:i))<=20) then
      sumE=0
      do j=2, 2*totalseqnum, 2
        !if (j==2*y) cycle     ! don't count parent sequence
        aa=index(alphabet, seqpair(j)(i:i)) 
        if (aa<=20 .AND. aa > 0) then
          sumE=sumE+1
          emit(k, aa)=emit(k, aa)+1
        end if
      end do
      if (sumE>0) emit(k, :)=emit(k, :)/sumE
      i=i+1
      k=k+1
    else
      i=i+1
    end if
  end do
  !dont need emission probs right now
  !-----output emission probabilities-----!
  !write (17, *)
  !write (17, *) "Emission probabilities:"
  !do i=1, Mnode
  !  write (17, "(20F5.2)") emit(i, :)
  !end do

  close(17)
  string=' rm *.drct *.gca ' 
  sysstatus=system(string)
 
  if(allocated(seqID)) deallocate(seqID)
  !write(*,*) "After deallocating seqID"
  if(allocated(seqpair)) deallocate(seqpair)
  !write(*,*) "After deallocating seqpair"
  if(allocated(HMM)) deallocate(HMM)
  !write(*,*) "After deallocating HMM"
  if(allocated(emit)) deallocate(emit)
  write(*,*) "done deallocating!"
  e=etime(time)
  write (*, *) "Elapsed:", e, "User:", time(1), "System:", time(2)
end program HMMMSA

!-----subroutine seq2drct-----!
subroutine seq2drct(seqfile)
  implicit none
  character(*), intent(in)::seqfile
  integer::i, sysstatus=0, system
  character(1000)::string, drctfile
  i=index(seqfile, ".")
  drctfile=seqfile(1:i-1) // ".drct"
  !write (*, *) "...converting sequence ", trim(seqfile(1:i-1)), " to drct file"
  string=' (./xseq2drct ' // trim(drctfile) // ' < ' // trim(seqfile) // " >> " // ' LOG) '! // ' >& /dev/null '
  sysstatus=system(string)
  if (sysstatus /= 0) stop "ERROR: System command for making drct files!!!"
  return
end subroutine seq2drct
!-----subroutine seq2gca------!
subroutine seq2gca(seqfile)
  use hmmstr
  implicit none
  character(*), intent(in)::seqfile
  integer::i, nres
  character(1000)::gcafile
  character(1)::chain

  real,dimension(:,:),allocatable::xgamma

  i = index(seqfile, ".")
  gcafile = seqfile(1:i-1) // ".gca"
      !hmmstr_gamma(xgamma,nres,mfile,seqfile,pfile,use_str,pdbfile,chain,gca)
  call hmmstr_gamma(xgamma, nres, "model_R.hmm", seqfile, " ", .FALSE., " ", chain, gcafile)
  return
end subroutine seq2gca
!-----subroutine drct2gca-----!
subroutine drct2gca(drctfile)
  use hmmstr
  implicit none
  character(*), intent(in)::drctfile
  integer::i, sysstatus=0, system, nres
  character(1000)::string, gcafile
  character(1)::chain

  real,dimension(:,:), allocatable::xgamma
  
  i=index(drctfile, ".")
  gcafile=drctfile(1:i-1) // ".gca"
      !hmmstr_gamma(xgamma,nres,mfile,seqfile,pfile,use_str,pdbfile,chain,gca)
  call hmmstr_gamma(xgamma, nres, "model_R.hmm", drctfile, " ", .FALSE., " ", chain, gcafile)
  !write (*, *) "...converting sequence ", trim(drctfile(1:i-1)), " to gca file"
  !string=' ./xhmmstr_gamma model_R.hmm ' // trim(drctfile) // " " // trim(gcafile) //  ' >> ' // ' LOG '
  !sysstatus=system(string)
  !if (sysstatus /= 0) stop "ERROR: System command for making gca files!!!"
  return
end subroutine drct2gca

!-----subroutine readgamma-----!
subroutine readgamma(filename, length, gammamatrix)
  implicit none
  character(*), intent(in)::filename
  integer, intent(in)::length
  real, dimension(length, 0:281), intent(out)::gammamatrix
  integer::x, y, status=0
  character(3000)::gammaline
  gammamatrix=0
  gammaline=" "
  x=0
  open (unit=99, file=filename, status="old")
  do
    read (99, "(A)", iostat=status) gammaline
    if (status/=0) exit
    if (gammaline(1:5)/="GAMMA") cycle
    x=x+1
    read (gammaline, "(7x, 282F9.6)", iostat=status) (gammamatrix(x, y), y=0, 281)
    if (status/=0) stop "ERROR: Bad gammaline!!!"
  end do
  close (99)
  return
end subroutine readgamma

!-----subroutine readdrtc-----!
subroutine readdrct(filename, length, profile)
  implicit none
  character(*), intent(in)::filename
  integer, intent(in)::length
  real*4, dimension (length, 20) ,intent(out)::profile
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
  end do
  close(99)
  return
end subroutine readdrct

!-----subroutine readmatrix-----!
subroutine readmatrix (filename, submatrix, a, b ,c)
  implicit none
  character(*), intent(in)::filename
  integer, intent(in)::a, b, c
  real, dimension(a, b, c), intent(out)::submatrix
  character(1000)::line
  integer::x, y, w, status=0
  open (unit=99, file=filename, status="old", iostat=status)
  if (status/=0) then
    write (*, "(A, $)") "ERROR: File not found -"
    write (*, *) trim(filename)
    stop
  end if
  submatrix=0
  w=1
  line=" "
  do
    do x=1, 20
      read (99, "(A)", iostat=status) line
      if (status/=0) exit
      if (line==" ") exit
      read (line, "(20F15.7)", iostat=status) (submatrix(x, y, w), y=1, 20)
      if (status/=0) stop "ERROR: Bad substitution matrix!!!"
    end do
    read (99, "(A)", iostat=status) line
    read (99, "(A)", iostat=status) line
    if (line==" ") w=w+1
    if (w>c) exit
  end do
  close (99, status="keep")
  return
end subroutine readmatrix

!-----subroutine makeMij-----!
subroutine makeMij(seqlen1, seqlen2, sequence1, sequence2, gmatrix1, gmatrix2, m1, m2, Mij)
  implicit none
  integer, intent(in)::seqlen1, seqlen2
  character(*), intent(in)::sequence1, sequence2
  real, dimension(seqlen1, 0:281), intent(inout)::gmatrix1
  real, dimension(seqlen2, 0:281), intent(inout)::gmatrix2
  real, dimension(20, 20, 281), intent(inout)::m1
  real, dimension(20, 281), intent(inout)::m2
  real, dimension(seqlen1, seqlen2), intent(out)::Mij
  real, dimension(20)::Pi
  integer::i, j, w, p1, p2
  real::sumM, sum1, sum2
  character(27)::alphabet="ACDEFGHIKLMNPQRSTVWYBZX.-~*"
  gmatrix1=gmatrix1+0.000000001
  gmatrix2=gmatrix2+0.000000001
  !normalize m1 and m2?
  do w=1, 281
    sum1=sum(m1(:, :, w))
    if (sum1>0) then
      m1(:, :, w)=m1(:, :, w)/sum1
    end if
  end do
  sum2=sum(m2(:, :))
  m2=m2/sum2
  !prob(residue) = sum over expected freq?
  Pi=0
  do i=1, 20
    do w=1, 281
      Pi(i)=Pi(i)+m2(i, w)
    end do
  end do
  !make substitution matrix for these sequences?
  Mij=0
  do i=1, seqlen1
    do j=1, seqlen2
      p1=index(alphabet, sequence1(i:i))
      p2=index(alphabet, sequence2(j:j))
      if (p1>20 .OR. p2>20 .OR. p1<1 .OR. P2<1) cycle
      sumM=0
      do w=1 ,281
        if (p1>=p2) then
          Mij(i, j)=Mij(i, j)+gmatrix1(i, w)*m1(p2, p1, w)*gmatrix2(j, w)
        else
          Mij(i, j)=Mij(i, j)+gmatrix1(i, w)*m1(p1, p2, w)*gmatrix2(j, w)
        end if
        sumM=sumM+gmatrix1(i, w)*gmatrix2(j, w)
      end do
      Mij(i, j)=Mij(i, j)/sumM
      if (p1==p2) then
        Mij(i, j)=log10(Mij(i, j)/(Pi(p1)*Pi(p2)))/log10(2.0)
      else
        Mij(i, j)=log10(Mij(i, j)/(2*Pi(p1)*Pi(p2)))/log10(2.0)
      end if
    end do
  end do
  return
end subroutine makeMij

!-----subroutine makedrctMij-----!
subroutine makedrctMij(seqlen1, seqlen2, profile1, profile2, gmatrix1, gmatrix2, m1, m2, Mij)
  implicit none
  integer, intent(in)::seqlen1, seqlen2
  real, dimension(seqlen1, 20), intent(inout)::profile1
  real, dimension(seqlen2, 20), intent(inout)::profile2
  real, dimension(seqlen1, 0:281), intent(inout)::gmatrix1
  real, dimension(seqlen2, 0:281), intent(inout)::gmatrix2
  real, dimension(20, 20, 281), intent(inout)::m1
  real, dimension(20, 281), intent(inout)::m2
  real, dimension(seqlen1, seqlen2), intent(out)::Mij
  real, dimension(20)::Pi
  integer::i, j, aa1, aa2, w, p1, p2
  real::sumM, sum1, sum2, score
  gmatrix1=gmatrix1+0.000000001
  gmatrix2=gmatrix2+0.000000001
  do w=1, 281
    sum1=sum(m1(:, :, w))
    if (sum1>0) then
      m1(:, :, w)=m1(:, :, w)/sum1
    end if
  end do
  sum2=sum(m2(:, :))
  m2=m2/sum2
  Pi=0
  do i=1, 20
    do w=1, 281
      Pi(i)=Pi(i)+m2(i, w)
    end do
  end do
  Mij=0
  do i=1, seqlen1
    do j=1, seqlen2
      do aa1=1, 20
        do aa2=1, 20
          if (profile1(i, aa1)>0 .and. profile2(j, aa2)>0) then
            sumM=0
            score=0
            do w=1 ,281
              if (aa1>=aa2) then
                score=score+gmatrix1(i, w)*m1(aa2, aa1, w)*gmatrix2(j, w)
              else
                score=score+gmatrix1(i, w)*m1(aa1, aa2, w)*gmatrix2(j, w)
              end if
              sumM=sumM+gmatrix1(i, w)*gmatrix2(j, w)
            end do
            score=score/sumM
            if (aa1==aa2) then
              score=log10(score/(Pi(aa1)*Pi(aa2)))/log10(2.0)
            else
              score=log10(score/(2*Pi(aa1)*Pi(aa2)))/log10(2.0)
            end if
            Mij(i, j)=Mij(i, j)+profile1(i, aa1)* profile2(j, aa2)*score
          end if
        end do
      end do
    end do
  end do
  return
end subroutine makedrctMij

!-----subroutine DP-----!
subroutine DP(seqlen1, seqlen2, opening, extension, Mij, traceback, score)
  implicit none
  integer, intent(inout)::seqlen1, seqlen2
  integer, intent(out)::score
  real, intent(in)::opening, extension
  real, dimension(1:seqlen1, 1:seqlen2), intent(inout)::Mij
  character(1), dimension(seqlen1+seqlen2), intent(out)::traceback
  integer::memstatus=0, i, j, len
  real::vm, vi, vd, maxM, maxI, maxD, maxV 
  integer, dimension(2)::pos
  real, allocatable::MS(:,:), IS(:,:), DS(:,:)
  character(1), allocatable::ML(:,:), IL(:,:), DL(:,:), trace(:,:)

  !-----make match, insertion, deletion, and trace matrix-----!
  seqlen1=seqlen1+1
  seqlen2=seqlen2+1
  allocate (MS(seqlen1, seqlen2), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating match score matrix!!!"
  allocate (ML(seqlen1, seqlen2), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating match letter matrix!!!"
  allocate (IS(seqlen1, seqlen2), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating insertion score matrix!!!"
  allocate (IL(seqlen1, seqlen2), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating insertion letter matrix!!!"
  allocate (DS(seqlen1, seqlen2), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating deletion score matrix!!!"
  allocate (DL(seqlen1, seqlen2), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating deletion letter matrix!!!"
  allocate (trace(seqlen1, seqlen2), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating traceback matrix!!!"
  MS=0
  IS=0
  DS=0
  ML=" "
  IL=" "
  DL=" "
  trace=" "
  MS(1, 1)=0
  MS(2:, 1)=0
  MS(1, 2:)=0
  IS(1, 1)=0
  IS(1, 2:)=-huge(1)
  IS(2:, 1)=-huge(1)
  DS(1, 1)=0
  DS(2:, 1)=-huge(1)
  DS(1, 2:)=-huge(1)
  vm=0
  vi=0
  vd=0
  maxM=0
  maxI=0
  maxD=0
  do j=2, seqlen2
    do i=2, seqlen1
      vm=MS(i-1, j-1)+Mij(i-1, j-1)
      vi=IS(i-1, j-1)+Mij(i-1, j-1)
      vd=DS(i-1, j-1)+Mij(i-1, j-1)
      maxM=Max(vm, vi, vd, 0.0)
      MS(i, j)=maxM
      if (maxM==vm) then
        ML(i, j)="M"
      elseif (maxM==vi) then
        ML(i, j)="I"
      elseif (maxM==vd) then
        ML(i, j)="D"
      else
        ML(i, j)="E"
      end if
      vm=MS(i-1, j)-opening
      vi=IS(i-1, j)-extension
      maxI=Max(vm, vi)
      IS(i, j)=maxI
      if (maxI==vm) then
        IL(i, j)="M"
      else
        IL(i, j)="I"
      end if
      vm=MS(i, j-1)-opening
      vi=IS(i, j-1)-opening
      vd=DS(i, j-1)-extension
      maxD=Max(vm, vi, vd)
      DS(i, j)=maxD
      if (maxD==vm) then
        DL(i, j)="M"
      elseif (maxD==vi) then
        DL(i, j)="I"
      else
        DL(i, j)="D"
      end if
    end do
  end do
  !-----Alignment, using loca matrix to record traceback (local alignment)-----!
  len=seqlen1+seqlen2
  traceback=" "
  pos=0
  maxV=0
  maxM=maxval(MS)
  maxI=maxval(IS)
  maxD=maxval(DS)
  maxV=Max(maxM, maxI, maxD)
  if (maxV==maxM) then
    traceback(1)="M"
    pos=maxloc(MS)
    seqlen1=pos(1)     !seqlen1 and seqlen2 become end positions
    seqlen2=pos(2)
  end if
  if (maxV==maxI) then
    traceback(1)="I"
    pos=maxloc(IS)
    seqlen1=pos(1)     !seqlen1 and seqlen2 become end positions
    seqlen2=pos(2)
  end if
  if (maxV==maxD) then
    traceback(1)="D"
    pos=maxloc(DS)
    seqlen1=pos(1)     !seqlen1 and seqlen2 become end positions
    seqlen2=pos(2)
  end if
  if (traceback(1)=="M") traceback(2)=ML(seqlen1, seqlen2)
  if (traceback(1)=="I") traceback(2)=IL(seqlen1, seqlen2)
  if (traceback(1)=="D") traceback(2)=DL(seqlen1, seqlen2)
  trace(seqlen1, seqlen2)=traceback(2)
  do i=2, len
    if (traceback(i)=="E") then
      if (traceback(i-2)=="M") then
        seqlen1=seqlen1+1
        seqlen2=seqlen2+1
        exit
      elseif (traceback(i-2)=="I") then
        seqlen1=seqlen1+1
        seqlen2=seqlen2
        exit
      elseif (traceback(i-2)=="D") then
        seqlen1=seqlen1
        seqlen2=seqlen2+1
        exit
      else
        exit
      end if
    end if
    if (traceback(i-1)=="M") then
      if (seqlen1<=2 .OR. seqlen2<=2) exit
      seqlen1=seqlen1-1
      seqlen2=seqlen2-1
      if (traceback(i)=="M") traceback(i+1)=ML(seqlen1, seqlen2)
      if (traceback(i)=="I") traceback(i+1)=IL(seqlen1, seqlen2)
      if (traceback(i)=="D") traceback(i+1)=DL(seqlen1, seqlen2)
      trace(seqlen1, seqlen2)=traceback(i+1)
    elseif (traceback(i-1)=="I") then
      if (seqlen1<=2) exit
      seqlen1=seqlen1-1
      seqlen2=seqlen2
      if (traceback(i)=="M") traceback(i+1)=ML(seqlen1, seqlen2)
      if (traceback(i)=="I") traceback(i+1)=IL(seqlen1, seqlen2)
      if (traceback(i)=="D") traceback(i+1)=DL(seqlen1, seqlen2)
      trace(seqlen1, seqlen2)=traceback(i+1)
    elseif (traceback(i-1)=="D") then
      if (seqlen2<=2) exit
      seqlen1=seqlen1
      seqlen2=seqlen2-1
      if (traceback(i)=="M") traceback(i+1)=ML(seqlen1, seqlen2)
      if (traceback(i)=="I") traceback(i+1)=IL(seqlen1, seqlen2)
      if (traceback(i)=="D") traceback(i+1)=DL(seqlen1, seqlen2)
      trace(seqlen1, seqlen2)=traceback(i+1)
    else
      cycle
    end if
  end do
  score=maxV
  deallocate(MS)
  deallocate(IS)
  deallocate(DS)
  deallocate(ML)
  deallocate(IL)
  deallocate(DL)
  deallocate(trace)
  return
end subroutine DP

!-----subroutine getdrctlen-----!
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

!-----subroutine BAA-----!
subroutine BAA(drct1, drct2, seqlen1, seqlen2, m1, m2, bias)
  implicit none
  character(*), intent(in)::drct1, drct2
  integer, intent(in)::seqlen1, seqlen2
  real, intent(in)::bias
  real, dimension(20, 20, 281), intent(inout)::m1
  real, dimension(20, 281), intent(inout)::m2
  integer::memstatus=0, aa, bb, block, k, i, j, sample, samplesize=1000, a, b, currentB, signal, l
  real(KIND=8)::logPmax, sumofscale, shiftscale, pS1S2, currentpK, MP, IP, DP, sumP, sumofhis, sum
  real::r
  character(27)::alphabet="ACDEFGHIKLMNPQRSTVWYBZX.-~*"
  character(1000)::gca1, gca2, outfile
  integer, dimension(1)::pos
  real, dimension(seqlen1, seqlen2)::alignM
  real, allocatable::pro1(:,:), pro2(:,:), gmatrix1(:,:), gmatrix2(:,:), Match(:,:)
  real(KIND=8), allocatable::matrixM(:,:,:), matrixD(:,:,:), matrixI(:,:,:), totalcount(:), count(:), scale(:), &
                             pK(:), logP(:), outprofile(:,:)
  integer*4::x, y
  integer*2::seq_list, resSeq
  integer, parameter::vrecsize=136
  real*4::phi_list, psi_list, omega_list, csum, gap_freq, ins_freq
  character::infile*1000, structure_name*4 , chainID*1, seq_resid*1, ss*1, rotamer*1, domainID*1, contex*1
  real*4, dimension(3)::xyzca
  real*4, dimension(20)::vall_pro
  integer*2, dimension(8)::clusterID

  write (*, *) "...aligning (BAA): ", trim(drct1), " ", trim(drct2)

  !-----make Mij-----!
  allocate (pro1(seqlen1, 20), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating pro1!!!"
  allocate (pro2(seqlen2, 20), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating pro2!!!"
  pro1=0
  pro2=0
  call readdrct(drct1, seqlen1, pro1)
  call readdrct(drct2, seqlen2, pro2)
    
  aa=index(drct1, ".")
  bb=index(drct2, ".")
  gca1=drct1(1:aa-1) // ".gca"
  gca2=drct2(1:bb-1) // ".gca"
  allocate (gmatrix1(seqlen1, 0:281), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating gmatrix1!!!"
  allocate (gmatrix2(seqlen2, 0:281), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating gmatrix2!!!"
  gmatrix1=0
  gmatrix2=0
  call readgamma(gca1, seqlen1, gmatrix1)
  call readgamma(gca2, seqlen2, gmatrix2)
  call makedrctMij(seqlen1, seqlen2, pro1, pro2, gmatrix1, gmatrix2, m1, m2, alignM)
  alignM=2*alignM+bias

  !-----define maximum number of block-----!
  block=min(seqlen1/10, seqlen2/10, 20)
  write (*, "(A, $)") " Maximum number of block:"
  write (*, *) block
  write (*, *)

  !-----algorithm: posterior probabilities of blocks and score matrices-----!
  !-----algorithm(1): allocate M, I, D matrices (3-dimension) and Match matrix)-----!
  allocate (matrixM(0:seqlen1, 0:seqlen2, 0:block), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating matrixM!!!"
  allocate (matrixD(0:seqlen1, 0:seqlen2, 0:block), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating matrixD!!!"
  allocate (matrixI(0:seqlen1, 0:seqlen2, 0:block), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating matrixI!!!"
  matrixM=0
  matrixD=0
  matrixI=0

  !-----algorithm(2): allocate totalcount, count, scale matrices-----!
  allocate (totalcount(0:block), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating totalcount!!!"
  totalcount=0
  totalcount(0)=1
  allocate (count(0:block), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating count!!!"
  count=0
  count(0)=1
  allocate (scale(0:block), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating scale!!!"
  scale=0
  scale(0)=1

  !-----algorithm(3): initialize  M, I and D matrices-----!
  matrixM(:, 0, 0)=0
  matrixI(:, 0, 0)=1
  matrixD(:, 0, 0)=0
  matrixM(:, 1:seqlen2, 0)=0
  matrixI(:, 1:seqlen2, 0)=0
  matrixD(:, 1:seqlen2, 0)=1
  do k=1, block
    matrixM(:, 0, k)=0
    matrixI(:, 0, k)=0
    matrixD(:, 0, k)=0
    matrixM(0, 1:seqlen2, k)=0
    matrixI(0, 1:seqlen2, k)=0
    matrixD(0, 1:seqlen2, k)=0
  end do

  !-----algorithm(4): count the number of alignments with k blocks -- totalcount-----!
  do k=1, block
    do i=1, seqlen1
      do j=1, seqlen2
        matrixM(i, j, k)=matrixM(i-1, j-1, k)+matrixI(i-1, j-1, k-1)+matrixD(i-1, j-1, k-1)
        matrixI(i, j, k)=matrixM(i-1, j, k)+matrixI(i-1, j, k)
        matrixD(i, j, k)=matrixM(i, j-1, k)+matrixI(i, j-1, k)+matrixD(i, j-1, k)
      end do
    end do
    totalcount(k)=matrixM(seqlen1, seqlen2, k)+matrixI(seqlen1, seqlen2, k)+matrixD(seqlen1, seqlen2, k)
  end do
  !write (*, *) "totalcount:", totalcount

  !-----algorithm(5): sum of the matrices with substitution scores -- count-----!
  !-----(consider matrix bias as a parameter)-----!
  matrixM=0
  matrixD=0
  matrixI=0
  matrixM(:, 0, 0)=0
  matrixI(:, 0, 0)=1
  matrixD(:, 0, 0)=0
  matrixM(:, 1:seqlen2, 0)=0
  matrixI(:, 1:seqlen2, 0)=0
  matrixD(:, 1:seqlen2, 0)=1
  do k=1, block
    matrixM(:, 0, k)=0
    matrixI(:, 0, k)=0
    matrixD(:, 0, k)=0
    matrixM(0, 1:seqlen2, k)=0
    matrixI(0, 1:seqlen2, k)=0
    matrixD(0, 1:seqlen2, k)=0
  end do
  do k=1, block
    do i=1, seqlen1
      do j=1, seqlen2
        matrixM(i, j, k)=(matrixM(i-1, j-1, k)+matrixI(i-1, j-1, k-1)+matrixD(i-1, j-1, k-1))*(2**alignM(i, j))
        matrixI(i, j, k)=matrixM(i-1, j, k)+matrixI(i-1, j, k)
        matrixD(i, j, k)=matrixM(i, j-1, k)+matrixI(i, j-1, k)+matrixD(i, j-1, k)
      end do
    end do
    count(k)=matrixM(seqlen1, seqlen2, k)+matrixI(seqlen1, seqlen2, k)+matrixD(seqlen1, seqlen2, k)
    scale(k)=1/count(k)
    !write (*, *) count(k), scale(k)
    matrixM(:, :, k)=matrixM(:, :, k)*scale(k)
    matrixI(:, :, k)=matrixI(:, :, k)*scale(k)
    matrixD(:, :, k)=matrixD(:, :, k)*scale(k)
    count(k)=matrixM(seqlen1, seqlen2, k)+matrixI(seqlen1, seqlen2, k)+matrixD(seqlen1, seqlen2, k)
  end do
  !write (*, *) "scale:", scale
  !write (*, *) "count:", count

  !-----posterior probability: probability of choosing the number of blocks - P(k|s1,s2)-----!
  !-----P(k|s1,s2)=(count(k)/totalcount(k))/(count(all)/totalcount(all))-----!
  !-----(scaling)-----!
  allocate (pK(0:block), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating pK!!!"
  pK=0
  allocate (logP(0:block), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating logP!!!"
  logP=0
  logPmax=0
  sumofscale=0
  do k=0, block
    sumofscale=sumofscale+log(scale(k))
    logP(k)=log(count(k)/totalcount(k))-sumofscale
    if (logP(k)>logPmax) logPmax=logP(k)
  end do
  shiftscale=logPmax-4
  pS1S2=0
  do k=0, block
    if (logP(k)<shiftscale) then
      pK(k)=0
    else
      pK(k)=exp(logP(k)-shiftscale)
    end if
    pS1S2=pS1S2+pK(k)
  end do
  do k=0, block
    pK(k)=pK(k)/pS1S2
  end do
  !write (*, *) "pS1S2:", pS1S2
  !write (*, *) "pK:", pK
  !write (*, *) "sum of pK:", sum(pK)

  !-----sampleback: loop SAMPLESIZE times, choose matrix biases => choose number of blocks-----!
  allocate (Match(0:seqlen1, 0:seqlen2), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating Match!!!"
  Match=0
  call INITRAND()
  do sample=1, samplesize
    a=seqlen1
    b=seqlen2
    call random_number(r)  !choose number of blocks
    currentpK=0
    currentB=0
    do k=0, block
      currentpK=currentpK+pK(k)
      if (r<currentpK) then
        currentB=k
        exit
      end if
    end do
    sumP=matrixM(a, b, currentB)+matrixI(a, b, currentB)+matrixD(a, b, currentB)  !sampleback
    MP=matrixM(a, b, currentB)/sumP
    IP=matrixI(a, b, currentB)/sumP
    DP=matrixD(a, b, currentB)/sumP
    signal=0
    do while (a>0 .AND. b>0 .AND. currentB>0)
      call random_number(r)
      if (r<MP) then
        Match(a, b)=Match(a, b)+1
        a=a-1
        b=b-1
        signal=1
        sumP=matrixM(a, b, currentB)+(matrixI(a, b, currentB-1)*scale(currentB))+(matrixD(a, b, currentB-1)*scale(currentB))
        MP=matrixM(a, b, currentB)/sumP
        IP=(matrixI(a, b, currentB-1)*scale(currentB))/sumP
        DP=(matrixD(a, b, currentB-1)*scale(currentB))/sumP
        !write (*, *) MP, IP, DP
      else if (r<MP+IP) then
        if (signal==1) currentB=currentB-1
        a=a-1
        signal=0
        sumP=matrixM(a, b, currentB)+matrixI(a, b, currentB)
        MP=matrixM(a, b, currentB)/sumP
        IP=matrixI(a, b, currentB)/sumP
        !write (*, *) MP, IP
      else
        if (signal==1) currentB=currentB-1
        b=b-1
        signal=0
        sumP=matrixM(a, b, currentB)+matrixI(a, b, currentB)+matrixD(a, b, currentB)
        MP=matrixM(a, b, currentB)/sumP
        IP=matrixI(a, b, currentB)/sumP
        DP=matrixD(a, b, currentB)/sumP
        !write (*, *) MP, IP, DP
      end if
    end do
  end do 
  Match=Match/samplesize
  allocate (outprofile(max(seqlen1, seqlen2), 20), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating outprofile!!!"
  sumofhis=0
  sum=0
  if (seqlen1>=seqlen2) then
    do i=1, seqlen1
      do j=1, 20
        do k=1, seqlen2
          !do l=1, 20
            !sumofhis=sumofhis+Match(i, k)*pro2(k, l)
            sumofhis=sumofhis+Match(i, k)*pro2(k, j)
          !end do
        end do
        outprofile(i, j)=0.5*pro1(i, j)+0.5*sumofhis
        sum=sum+outprofile(i, j)
        sumofhis=0
      end do
      outprofile(i, :)=outprofile(i, :)/sum
      sumofhis=0
      sum=0
    end do 
  else
    do i=1, seqlen2
      do j=1, 20
        do k=1, seqlen1
          !do l=1, 20
            !sumofhis=sumofhis+Match(k, i)*pro1(k, l)
            sumofhis=sumofhis+Match(k, i)*pro1(k, j)
          !end do
        end do
        outprofile(i, j)=0.5*pro2(i, j)+0.5*sumofhis
        sum=sum+outprofile(i, j)
        sumofhis=0
      end do
      outprofile(i, :)=outprofile(i, :)/sum
      sumofhis=0
      sum=0
    end do
  end if
  outfile=trim(drct1(1:aa-1)) // "-" // trim(drct2(1:bb-1)) // ".drct"
  open(unit=99, file=outfile, status="replace", form="unformatted", access="direct", recl=vrecsize)
  seq_list=0
  resSeq=0
  phi_list=0
  psi_list=0
  omega_list=0
  csum=1
  gap_freq=0 
  ins_freq=0
  structure_name=" "
  chainID=" "
  seq_resid=" "
  ss=" "
  rotamer=" "
  domainID=" "
  contex=" "
  xyzca=0
  vall_pro=0
  clusterID=0
  do x=1, max(seqlen1, seqlen2)
    pos=maxloc(outprofile(x, :))
    seq_resid=alphabet(pos(1):pos(1))
    vall_pro=outprofile(x, :)
    seq_list=x
    resSeq=x
    write (99, rec=x) seq_list, resSeq, structure_name, chainID, seq_resid, ss, rotamer, (xyzca(y), y=1, 3), &
          phi_list, psi_list, omega_list, csum, (vall_pro(y), y=1, 20), gap_freq, ins_freq, domainID, contex, &
          (clusterID(y), y=1, 3)
  end do
  close(99)
  !if (seqlen1>=seqlen2) then 
  !  do i=1, seqlen1
  !    write (*, "(20F13.10)") (outprofile(i, j), j=1, 20)
  !  end do
  !else
  !  do i=1, seqlen2
  !    write (*, "(20F13.10)") (outprofile(i, j), j=1, 20)
  !  end do
  !end if
  deallocate(pro1)
  deallocate(pro2)
  deallocate(gmatrix1)
  deallocate(gmatrix2)
  deallocate(matrixM)
  deallocate(matrixD)
  deallocate(matrixI)
  deallocate(Match)
  deallocate(totalcount)
  deallocate(count)
  deallocate(pK)
  deallocate(scale)
  deallocate(logP)
  deallocate(outprofile)
  return
end subroutine BAA

!-----subroutine INITRAND-----!
subroutine INITRAND()
  implicit none
  integer::i, j, msec
  real::x
  integer, dimension(8)::dt
  call date_and_time(values=dt)
  msec=1000*dt(7)+dt(8)
  call random_seed(size=i)
  call random_seed(put=(/(msec, j=1, i)/))
  call random_number(x)
end subroutine INITRAND

!-----subroutine readgcaseq-----!
subroutine readgcaseq(filename, sequence)
  implicit none
  character(*), intent(in)::filename
  character(*), intent(out)::sequence
  integer::i, p, status=0
  character(27)::alphabet="ACDEFGHIKLMNPQRSTVWYBZX.-~*"
  character(1000)::line
  sequence=" "
  line=" "
  open (unit=16, file=filename, status="old", iostat=status)
  if (status/=0) then
    write (*, "(A, $)") "ERROR: gca file not found -"
    write (*, *) trim(filename)
    stop
  end if
  do
    read (16, "(A)", iostat=status) line
    if (status/=0) exit
    if (line(1:7)/="RESIDUE") cycle
    do i=9, 20
      p=index(alphabet, line(i:i))
      if (p>=1 .AND. p<=20) sequence=trim(sequence) // trim(line(i:i))
    end do
  end do
  close (16)
  return
end subroutine readgcaseq

!-----subroutine makeMSA-----!
subroutine makeMSA(seqnumber, pairs)
  implicit none
  integer, intent(in)::seqnumber
  !tried this: 
  !character(9999), dimension(seqnumber*2),allocatable,  intent(inout)::pairs(:)
  !old code:
  character(9999), dimension(seqnumber*2), intent(inout)::pairs
  integer::i, j, m, n, len1, len2, L1, L2, g, h, x, y, s1, s2, t1, t2, k
  character(27)::alphabet="ACDEFGHIKLMNPQRSTVWYBZX.-~*"
  character(9999)::seq1, seq2, aline
  do i=2, seqnumber
    j=0
    seq1=" "
    j=len_trim(pairs(1))
    seq1=pairs(1)
    m=1
    do n=1, j
      if (seq1(n:n)/=" " .and. seq1(n:n)/=".") then
        seq1(m:m)=seq1(n:n)
        m=m+1
      end if
    end do
    seq1(m:)=" "
    seq2=" "
    j=len_trim(pairs(2*i-1))
    seq2=pairs(2*i-1)
    m=1
    do n=1, j
      if (seq2(n:n)/=" " .and. seq2(n:n)/=".") then
        seq2(m:m)=seq2(n:n)
        m=m+1
      end if
    end do
    seq2(m:)=" "
    len1=len_trim(pairs(1))
    len2=len_trim(pairs(2*i-1))
    L1=len_trim(seq1)
    L2=len_trim(seq2)
    aline=" "
    do g=1, L1-9
      aline=seq1(g:g+9)
      do j=1, Len2-9
        x=0
        y=1
        do while (y<=10)
          h=j+x
          x=x+1
          if (pairs(2*i-1)(h:h)==".") cycle
          if (pairs(2*i-1)(h:h)/=aline(y:y)) exit
          y=y+1
        end do
        if (y==11) exit
      end do
      if (y==11) exit
    end do
    do
      if (pairs(2*i-1)(j:j)==".") then
        j=j+1
      else
        if (g>1) then
          s1=0
          exit
        else
          s1=j
          exit
        end if
      end if
    end do
    aline=" "
    do g=1, L2-9
      aline=seq2(g:g+9)
      do j=1, Len1-9
        x=0
        y=1
        do while (y<=10)
          h=j+x
          x=x+1
          if (pairs(1)(h:h)==".") cycle
          if (pairs(1)(h:h)/=aline(y:y)) exit
          y=y+1
        end do
        if (y==11) exit
      end do
      if (y==11) exit
    end do
    do
      if (pairs(1)(j:j)==".") then
        j=j+1
      else
        if (g>1) then
          s2=0
          exit
        else
          s2=j
          exit
        end if
      end if
    end do
    aline=" "
    do g=L1-9, 1, -1
      aline=seq1(g:g+9)
      do j=Len2, 10, -1
        x=0
        y=10
        do while (y>=1)
          h=j+x
          x=x-1
          if (pairs(2*i-1)(h:h)==".") cycle
          if (pairs(2*i-1)(h:h)/=aline(y:y)) exit
          y=y-1
        end do
        if (y==0) exit
      end do
      if (y==0) exit
    end do
    do
      if (pairs(2*i-1)(j:j)==".") then
        j=j-1
      else
        if (g<L1-9) then
          t1=0
          exit
        else
          t1=j
          exit
        end if
      end if
    end do
    aline=" "
    do g=L2-9, 1, -1
      aline=seq2(g:g+9)
      do j=Len1, 10, -1
        x=0
        y=10
        do while (y>=1)
          h=j+x
          x=x-1
          if (pairs(1)(h:h)==".") cycle
          if (pairs(1)(h:h)/=aline(y:y)) exit
          y=y-1
        end do
        if (y==0) exit
      end do
      if (y==0) exit
    end do
    do
      if (pairs(1)(j:j)==".") then
        j=j-1
      else
        if (g<L2-9) then
          t2=0
          exit
        else
          t2=j
          exit
        end if
      end if
    end do
    !write (*, *) s1, s2
    !write (*, *) t1, t2
    if (s1<=s2) then
      do j=1, s2-1
        pairs(2*i-1)(2:)=trim(pairs(2*i-1)(1:))
        pairs(2*i-1)(1:1)="."
        pairs(2*i)(2:)=trim(pairs(2*i)(1:))
        pairs(2*i)(1:1)="."
      end do
    else
      do j=1, s1-1
        pairs(1)(j+1:)=trim(pairs(1)(j:))
        pairs(1)(j:j)=pairs(2*i-1)(j:j)
        do k=2, 2*i-2
          pairs(k)(2:)=trim(pairs(k)(1:))
          pairs(k)(1:1)="."
        end do
      end do
    end if
    j=max(s1, s2)
    do 
      !write (*, *) j, index(alphabet, pairs(1)(j:j)), index(alphabet, pairs(2*i-1)(j:j))
      if (index(alphabet, pairs(1)(j:j))<1) then
        if (index(alphabet, pairs(2*i-1)(j:j))>=1) then
          pairs(1)(j:j)=pairs(2*i-1)(j:j)  
          do k=2, 2*i-2
            pairs(k)(j:j)="."
          end do
          j=j+1
          cycle
        elseif (index(alphabet, pairs(2*i-1)(j:j))<1) then
          exit
        end if
      elseif (index(alphabet, pairs(2*i-1)(j:j))<1) then
        if (index(alphabet, pairs(1)(j:j))>=1) then
          pairs(2*i-1)(j:j)="."
          pairs(2*i)(j:j)="."
          j=j+1
          cycle
        elseif (index(alphabet, pairs(1)(j:j))<1) then
          exit
        end if
      elseif (index(alphabet, pairs(1)(j:j))==index(alphabet, pairs(2*i-1)(j:j))) then
        j=j+1
        cycle
      elseif (index(alphabet, pairs(1)(j:j))>20 .and. index(alphabet, pairs(2*i-1)(j:j))<=20) then
        pairs(2*i-1)(j+1:)=trim(pairs(2*i-1)(j:))
        pairs(2*i-1)(j:j)="."
        pairs(2*i)(j+1:)=trim(pairs(2*i)(j:))
        pairs(2*i)(j:j)="."
        j=j+1
        cycle
      elseif (index(alphabet, pairs(1)(j:j))<=20 .and. index(alphabet, pairs(2*i-1)(j:j))>20) then
        do k=1, 2*i-2
          pairs(k)(j+1:)=trim(pairs(k)(j:))
          pairs(k)(j:j)="."
        end do
        j=j+1
        cycle
      end if
    end do
  end do
end subroutine makeMSA

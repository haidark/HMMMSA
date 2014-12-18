!-------------------------------------------------------------------------------
!Heading:           Bayesian adaptive sequence alignment algorithm (matrix bias is NOT a posterior)
!Date:              2005/09/19
!Version:           1.0.0
!Author:            Yao-ming huang
!Purpose:           Running HMMSUM-D model on Bayesian adaptive sequence alignment algorithm
!Usage:
!File Format:       Input: gca/drct sequence files for query and target, observed counts, expected counts(18)
!                   Output: ali file(19), histogram file(20), histogram ps file(21)
!Note:              Use HMMSUM-D model
!                   Set bias to 0
!Reference:
!Restriction:
!Revision History:  2005/11/16: V. 2.0.0
!                     This program can take two drct profiles as inputs.
!                       -Add three subroutines: drct2gca, readdrct, makedrctMij
!                   2006/09/18
!                     Bug fix in subroutine "makeMij" and "makedrctMij"
!Copyright (C) 2005 Yao-ming Huang.
!All Rights Reserved.
!Permission to use, copy, modify, and distribute this software and its documentation for non-profit purposes
!and without fee is hereby granted provided that this copyright notice appears in all copies.
!-------------------------------------------------------------------------------

program bayesalign
  implicit none
  integer::openstatus=0, readstatus=0, memstatus=0, argnum, s1, s2, i, j, k, a, b, g, h, x, y, w, &
           seqlen1, seqlen2, block, samplesize=1000, sample, currentB, pos, signal, blo
  integer(KIND=1)::iargc
  real::e, etime, time(2), r, bias, score, tmp, color, maxv, minv, pseudo=0.000000001
  real(KIND=8)::pS1S2, currentpK, sumP, MP, IP, DP, logPmax, sumofscale, shiftscale
  character(27)::aa="ACDEFGHIKLMNPQRSTVWYBZX-~.*"
  character(1000)::temp1, temp2, infile1, infile2, alifile, hisfile, his2psfile, aline, gcafile1, gcafile2
  character(9999)::seq1, seq2
  real, dimension(20,20,281)::matrix1
  real, dimension(20,281)::matrix2
  character(9999), dimension(2)::ali
  real, allocatable::gammaM1(:,:), gammaM2(:,:), alignM(:,:), Match(:,:), profile1(:,:), profile2(:,:)
  real(KIND=8), allocatable::totalcount(:), count(:), scale(:), matrixM(:,:,:), matrixD(:,:,:), &
                             matrixI(:,:,:), pK(:), logP(:)
  !integer, allocatable::alignA(:,:)

  !-------------------------------------------------------------------------------
  !open gca files and read sequence from gca files
  !-------------------------------------------------------------------------------
  e=etime(time)
  argnum=iargc()
  if (argnum==0 .or. mod(argnum, 2)/=0) then
    write (*, *) "xbayesalign_D   V. 2.0.0 Nov 2005"
    write (*, *) "Usage:"
    write (*, *) "  xbayesalign_D Query Target -b {Bias}"
    write (*, *) "  -> Query / Target: gca/drct files for query and target"
    write (*, *) "  -> Bias: Matrix bias, default -> 0"
    stop
  end if
  bias=0
  call getarg(1, temp1)
  call getarg(2, temp2)
  infile1=trim(temp1)
  infile2=trim(temp2)
  do i=3, argnum, 2
    call getarg(i, temp1)
    call getarg(i+1, temp2)
    if (temp1=="-b") then
      read (temp2, *, iostat=readstatus) bias
      if (readstatus/=0) stop "Wrong argument in bias!!!"
    else
      stop "Wrong argument(s)!!!"
    end if
  end do
  s1=index(infile1, ".")
  s2=index(infile2, ".")
  alifile=infile1(1:s1-1) // "-" // infile2(1:s2-1) // ".ali"
  hisfile=infile1(1:s1-1) // "-" // infile2(1:s2-1) // ".his"
  his2psfile=infile1(1:s1-1) // "-" // infile2(1:s2-1) // ".ps"
  seq1=" "
  seq2=" "
  if (index(infile1, ".gca")/=0 .and. index(infile2, ".gca")/=0) then
    call readgcaseq(infile1, seq1)
    call readgcaseq(infile2, seq2)
    seqlen1=len_trim(seq1)
    seqlen2=len_trim(seq2)
    write (*, "(A, $)") "Sequence of"
    write (*, *) trim(infile1)
    write (*, "(A, $)") trim(seq1)
    write (*, *)
    write (*, "(A, $)") "length:" 
    write (*, *) seqlen1
    write (*, *)
    write (*, "(A, $)") "Sequence of"
    write (*, *) trim(infile2)
    write (*, "(A, $)") trim(seq2)
    write (*, *)
    write (*, "(A, $)") "length:"
    write (*, *) seqlen2
    write (*, *)
    write (*, *) "...aligning!!!"
    write (*, *)
  elseif (index(infile1, ".drct")/=0 .and. index(infile2, ".drct")/=0) then
    call drct2gca(infile1)
    call drct2gca(infile2)
    gcafile1=infile1(1:s1-1) // ".gca"
    gcafile2=infile2(1:s2-1) // ".gca"
    call readgcaseq(gcafile1, seq1)
    call readgcaseq(gcafile2, seq2)
    seqlen1=len_trim(seq1)
    seqlen2=len_trim(seq2)
    allocate (profile1(seqlen1, 20), stat=memstatus)
    if (memstatus/=0) stop "ERROR: Allocating profile1!!!"
    allocate (profile2(seqlen2, 20), stat=memstatus)
    if (memstatus/=0) stop "ERROR: Allocating profile2!!!"
    call readdrct(infile1, seqlen1, profile1)
    call readdrct(infile2, seqlen2, profile2)
    write (*, *) "...aligning!!!"
    write (*, *)
  else
    stop "Wrong input file(s)!!!"
  end if

  !-------------------------------------------------------------------------------
  !open and read observed(matrix1) & background(matrix2) counts
  !-------------------------------------------------------------------------------
  matrix1=0
  matrix2=0
  call readmatrix("observed", matrix1, 20, 20, 281)
  open (unit=18, file="expected", status="old", iostat=openstatus)
  if (openstatus/=0) then
    write (*, "(A, $)") "ERROR: File for expected counts not found"
    stop
  end if
  w=1
  do
    read (18, "(A)", iostat=readstatus) aline
    if (readstatus/=0) exit
    if (aline==" ") cycle
    read (aline, "(20F15.7)", iostat=readstatus) (matrix2(x, w), x=1, 20)
    if (readstatus/=0) stop "ERROR: Bad matrix2!!!"
    w=w+1
    if (w>281) exit
  end do
  close(18)

  !-------------------------------------------------------------------------------
  !add pseudocounts
  !-------------------------------------------------------------------------------
  matrix1=matrix1+(pseudo**2) !observed counts + pseudocounts**2
  matrix2=matrix2+pseudo      !background counts + pseudocounts

  !-------------------------------------------------------------------------------
  !open gcb files and reading gamma values into memory
  !-------------------------------------------------------------------------------
  allocate (gammaM1(seqlen1, 0:281), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating gammaM1!!!"
  allocate (gammaM2(seqlen2, 0:281), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating gammaM2!!!"
  gammaM1=0
  gammaM2=0
  if (index(infile1, ".gca")/=0 .and. index(infile2, ".gca")/=0) then
    call readgamma(infile1, seqlen1, gammaM1)
    call readgamma(infile2, seqlen2, gammaM2)
  else
    call readgamma(gcafile1, seqlen1, gammaM1)
    call readgamma(gcafile2, seqlen2, gammaM2)
  end if

  !-------------------------------------------------------------------------------
  !produce Mij (substitution matrix weighted by gamma values)
  !-------------------------------------------------------------------------------
  allocate (alignM(seqlen1, seqlen2), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating Mij matrix!!!"
  alignM=0
  if (index(infile1, ".gca")/=0 .and. index(infile2, ".gca")/=0) then
    call makeMij(seqlen1, seqlen2, seq1, seq2, gammaM1, gammaM2, matrix1, matrix2, alignM)
  else
    call makedrctMij(seqlen1, seqlen2, profile1, profile2, gammaM1, gammaM2, matrix1, matrix2, alignM)
  end if
  alignM=2*alignM+bias
  !write (*, *) "Max score in Mij:", Maxval(alignM)
  !write (*, *) "Mix score in Mij:", Minval(alignM)
  !write (*, *) "Average score in Mij:", Sum(alignM)/(seqlen1*seqlen2)

  !-------------------------------------------------------------------------------
  !produce Aij (match=1, mismatch=0 -- for debuging only)
  !-------------------------------------------------------------------------------
  !allocate (alignA(seqlen1, seqlen2), stat=memstatus)
  !if (memstatus/=0) stop "ERROR: Allocating Aij matrix!!!"
  !alignA=0
  !call makeAij(seqlen1, seqlen2, seq1, seq2, alignA)
  !write (*, "(A, $)") " "
  !do i=1, seqlen2
  !  write (*, "(A1, $)") seq2(i:i)
  !end do
  !write (*, *)
  !k=1
  !do i=1, seqlen1
  !  write (*, "(A1, $)") seq1(k:k)
  !  do j=1, seqlen2
  !    write (*, "(I1, $)") alignA(i, j)
  !  end do
  !  write (*, *)
  !  k=k+1
  !end do

  !-------------------------------------------------------------------------------
  !define maximum number of block
  !-------------------------------------------------------------------------------
  block=min(seqlen1/10, seqlen2/10, 20)
  write (*, "(A, $)") " Maximum number of block:"
  write (*, *) block
  write (*, *)

  !-------------------------------------------------------------------------------
  !algorithm: posterior probabilities of blocks and score matrices
  !algorithm(1): allocate M, I, D matrices (3-dimension) and Match matrix
  !-------------------------------------------------------------------------------
  allocate (matrixM(0:seqlen1, 0:seqlen2, 0:block), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating matrixM!!!"
  allocate (matrixD(0:seqlen1, 0:seqlen2, 0:block), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating matrixD!!!"
  allocate (matrixI(0:seqlen1, 0:seqlen2, 0:block), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating matrixI!!!"
  matrixM=0
  matrixD=0
  matrixI=0

  !-------------------------------------------------------------------------------
  !algorithm(2): allocate totalcount, count, scale matrices
  !totalcount: record the number of alignments with k blocks
  !count: record sum of the matrices (with substitution scores)
  !scale: record scaling factors
  !-------------------------------------------------------------------------------
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

  !-------------------------------------------------------------------------------
  !algorithm(3): initialize  M, I and D matrices
  !-------------------------------------------------------------------------------
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

  !-------------------------------------------------------------------------------
  !algorithm(4): count the number of alignments with k blocks -- totalcount
  !-------------------------------------------------------------------------------
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

  !-------------------------------------------------------------------------------
  !algorithm(5): sum of the matrices with substitution scores -- count
  !(consider matrix bias as a parameter)
  !-------------------------------------------------------------------------------
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
    !write (*, *) "count(k) and scale(k): ", count(k), scale(k)
    matrixM(:, :, k)=matrixM(:, :, k)*scale(k)
    matrixI(:, :, k)=matrixI(:, :, k)*scale(k)
    matrixD(:, :, k)=matrixD(:, :, k)*scale(k)
    count(k)=matrixM(seqlen1, seqlen2, k)+matrixI(seqlen1, seqlen2, k)+matrixD(seqlen1, seqlen2, k)
  end do
  !write (*, *) "scale:", scale
  !write (*, *) "count:", count

  !-------------------------------------------------------------------------------
  !posterior probability: probability of choosing the number of blocks - P(k|s1,s2)
  !P(k|s1,s2)=(count(k)/totalcount(k))/(count(all)/totalcount(all))
  !(scaling)
  !-------------------------------------------------------------------------------
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

  !-------------------------------------------------------------------------------
  !sampleback: loop SAMPLESIZE times, choose matrix biases => choose number of blocks
  !-------------------------------------------------------------------------------
  open (unit=19, file=alifile, status="replace", iostat=openstatus)
  allocate (Match(0:seqlen1, 0:seqlen2), stat=memstatus)
  if (memstatus/=0) stop "ERROR: Allocating Match!!!"
  Match=0
  call INITRAND()
  do sample=1, samplesize
    score=0
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
    ali=" "
    pos=1
    signal=0
    do while (a>0 .AND. b>0 .AND. currentB>0)
      call random_number(r)
      if (r<MP) then
        score=score+alignM(a, b)+bias
        ali(1)(pos:pos)=seq1(a:a)
        ali(2)(pos:pos)=seq2(b:b)
        Match(a, b)=Match(a, b)+1
        a=a-1
        b=b-1
        pos=pos+1
        signal=1
        sumP=matrixM(a, b, currentB)+(matrixI(a, b, currentB-1)*scale(currentB))+(matrixD(a, b, currentB-1)*scale(currentB))
        MP=matrixM(a, b, currentB)/sumP
        IP=(matrixI(a, b, currentB-1)*scale(currentB))/sumP
        DP=(matrixD(a, b, currentB-1)*scale(currentB))/sumP
        !write (*, *) MP, IP, DP
      else if (r<MP+IP) then 
        if (signal==1) currentB=currentB-1
        ali(1)(pos:pos)=seq1(a:a)
        ali(2)(pos:pos)="-"
        a=a-1
        pos=pos+1
        signal=0
        sumP=matrixM(a, b, currentB)+matrixI(a, b, currentB)
        MP=matrixM(a, b, currentB)/sumP
        IP=matrixI(a, b, currentB)/sumP
        !write (*, *) MP, IP
      else
        if (signal==1) currentB=currentB-1
        ali(1)(pos:pos)="-"
        ali(2)(pos:pos)=seq2(b:b)
        b=b-1
        pos=pos+1
        signal=0
        sumP=matrixM(a, b, currentB)+matrixI(a, b, currentB)+matrixD(a, b, currentB)
        MP=matrixM(a, b, currentB)/sumP
        IP=matrixI(a, b, currentB)/sumP
        DP=matrixD(a, b, currentB)/sumP
        !write (*, *) MP, IP, DP
      end if
    end do

  !-------------------------------------------------------------------------------
  !alignments output -- ali file
  !-------------------------------------------------------------------------------
    tmp=pos
    blo=ceiling(tmp/50)
    !write (19, "(A, $)") "SCORE "
    !write (19, "(F, $)") score
    !write (19, "(A, $)") "  "
    write (19, "(A, $)") "Bias "
    write (19, "(F4.1, $)") bias
    write (19, "(A, $)") "  "
    write (19, "(A, $)") "Block# "
    write (19, "(I2, $)") k
    write (19, *)
    do g=1, seqlen1-19, 1
      aline=seq1(g:g+19)
      x=0
      y=1
      do while (y<=20)
        x=x+1
        if (ali(1)(pos-x:pos-x)=="-" .OR. ali(1)(pos-x:pos-x)==" ") cycle
        if (aline(y:y)/=ali(1)(pos-x:pos-x)) exit
        y=y+1
      end do
      if (y==21) exit
    end do
    do h=1, seqlen2-19, 1
      aline=seq2(h:h+19)
      x=0
      y=1
      do while (y<=20)
        x=x+1
        if (ali(2)(pos-x:pos-x)=="-" .OR. ali(2)(pos-x:pos-x)==" ") cycle
        if (aline(y:y)/=ali(2)(pos-x:pos-x)) exit
        y=y+1
      end do
      if (y==21) exit
    end do
    do i=1, blo
      write (19, "(A, $)") "QUERY "
      if (g>seqlen1) then 
        write (19, "(I4, $)") seqlen1
      else
        write (19, "(I4, $)") g
      end if
      write (19, "(A, $)") "   "
      do j=1, 50
        if (pos-j>0) then
          write (19, "(A, $)") ali(1)(pos-j:pos-j)
          if (ali(1)(pos-j:pos-j)/="-") g=g+1
        else
          write (19, "(A, $)") " "
        end if
      end do
      g=g-1
      write (19, "(I4, $)") g
      write (19, *)
      g=g+1
      write (19, "(A, $)") "TEMP  "
      if (h>seqlen2) then
        write (19, "(I4, $)") seqlen2
      else
        write (19, "(I4, $)") h
      end if
      write (19, "(A, $)") "   "
      do j=1, 50
        if (pos-j>0) then
          write (19, "(A, $)") ali(2)(pos-j:pos-j)
          if (ali(2)(pos-j:pos-j)/="-") h=h+1
        else 
          write (19, "(A, $)") " "
        end if
      end do
      h=h-1
      write (19, "(I4, $)") h
      write (19, *)
      h=h+1
      pos=pos-50
    end do
    write (19, *)
    !write (*, *) trim(ali(1)(:))
    !write (*, *) trim(ali(2)(:))
  end do
  close(19)

  !-------------------------------------------------------------------------------
  !histogram output - marginal probability; P(Ai,j=1|R) -- histogram file
  !-------------------------------------------------------------------------------
  open (unit=20, file=hisfile, status="replace", iostat=openstatus)
  Match=Match/samplesize
  do i=1, seqlen1
    do j=1, seqlen2
      write (20, "(I4, $)") i
      write (20, "(A, $)") ", "
      write (20, "(I4, $)") j
      write (20, "(A, $)") ", "
      write (20, "(F9.7, $)") Match(i, j)
      write (20 ,*)
    end do
  end do
  close(20)

  !-------------------------------------------------------------------------------
  !histogram ps output -- histogram ps file
  !-------------------------------------------------------------------------------
  x=0
  y=0
  maxv=maxval(Match)
  minv=minval(Match)
  if ((minv-maxv)==0) stop "No alignment in Bayes part!!!"
  open (unit=21, file=his2psfile, status="replace", iostat=openstatus)
  if (openstatus/=0) stop "ERROR: Opening his2ps file!!!"
  write (21, '("%!PS")')
  write (21, '("% ---- define box procedure ----")')
  write (21, '("/hline")')
  write (21, '("{ 9 0 rlineto")')
  write (21, '("  0 2 rlineto")')
  write (21, '("  stroke} def")')
  write (21, '("/vline")')
  write (21, '("{ 0 9 rlineto")')
  write (21, '("  2 0 rlineto")')
  write (21, '("  stroke} def")')
  write (21, '("/box")')
  write (21, '("{ 0 9 rlineto")')
  write (21, '("  9 0 rlineto")')
  write (21, '("  0 -9 rlineto")')
  write (21, '("  closepath } def")')
  write (21, '("/graybox")')
  write (21, '("{ newpath setgray moveto box fill } def")')
  write (21, '("/hsbbox")')
  write (21, '("{ newpath moveto box sethsbcolor fill } def")')
  write (21, '("% ---- begin program ---- %")')
  write (21, '("/Courier findfont")')
  write (21, '(" 10 scalefont setfont")')
  write (21, '(" 10 10 translate")')
  write (21, '(" 0 0 moveto")')
  do i=1, seqlen1
    write (21, '("hline newpath")')
    x=i*9
    y=0
    write (21, '(2I8, " moveto")') x-7, -7
    write (21, '("(", A1, ")", " show ")') seq1(i:i)
    write (21, '(2I8, " moveto")') x, y
  end do
  x=x+9
  write (21, '(2I8, " moveto")') x, y
  write (21, '("(", A20, ")", " show ")') trim(infile1) 
  write (21, '(2I8, " moveto")') 0, 0
  do j=1, seqlen2
    write (21, '("vline newpath")')
    x=0
    y=j*9
    write (21, '(2I8, " moveto")') -7, y-7
    write (21, '("(", A1, ")", " show ")') seq2(j:j)
    write (21, '(2I8, " moveto")') x, y
  end do
  y=y+9
  write (21, '(2I8, " moveto")') x, y
  write (21, '("(", A20, ")", " show ")') trim(infile2)
  do i=1, seqlen1
    do j=1, seqlen2
      write (21, "(A, 2I4)") "% i, j", i, j
      write (21, "(A, 1F9.7)") "% Probability", Match(i, j)
    end do
  end do
  do i=1, seqlen1
    do j=1, seqlen2
      if (Match(i, j)/=0) then
        x=i*9-9
        y=j*9-9
        !color=0.67-(0.67*Match(i, j)/maxv)
        color=0.67*(Match(i, j)-maxv)/(minv-maxv)
        write (21, '(5F8.2, " hsbbox ")') color, 1.0, 1.0, real(x), real(y)
      end if
    end do
  end do
  write (21, '(F8.2, " setgray ")') 0.0
  do i=1, seqlen1
    do j=1, seqlen2
      if (Match(i, j)/=0) then
        x=i*9-9
        y=j*9-9
        write (21, '(2I8, " moveto ")') x, y
        write (21, '(" box ")')
        write (21, '(" stroke ")')
      end if
    end do
  end do
  write (21, '("%%%")')
  write (21, '("showpage")')
  close (21)
  !deallocate(alignA)
  deallocate(alignM)
  deallocate(gammaM1)
  deallocate(gammaM2)
  deallocate(matrixM)
  deallocate(matrixD)
  deallocate(matrixI)
  deallocate(Match)
  deallocate(totalcount)
  deallocate(count)
  deallocate(pK)
  deallocate(scale)
  deallocate(logP)
  if (allocated(profile1)) deallocate(profile1)
  if (allocated(profile2)) deallocate(profile2)
  e=etime(time)
  write (*, *) "Elapsed:", e, "User:", time(1), "System:", time(2)
end program bayesalign
  
!-----subroutine readgcaseq-----!
subroutine readgcaseq(filename, seq)
  implicit none
  character(*), intent(in)::filename
  character(*), intent(out)::seq
  integer::i, p, status=0
  character(27)::aa="ACDEFGHIKLMNPQRSTVWYBZX-~.*"
  character(1000)::line
  seq=" "
  line=" "
  open (unit=16, file=filename, status="old", iostat=status)
  if (status/=0) then
    write (*, "(A, $)") "ERROR: gcb file not found -"
    write (*, *) trim(filename)
    stop
  end if
  do
    read (16, "(A)", iostat=status) line
    if (status/=0) exit
    if (line(1:7)/="RESIDUE") cycle
    do i=9, 20
      p=index(aa, line(i:i))
      if (p>=1 .AND. p<=20) seq=trim(seq) // trim(line(i:i))
    end do
  end do
  close (16)
  return
end subroutine readgcaseq

!-----subroutine drct2gca-----!
subroutine drct2gca(drctfile)
  implicit none
  character(*), intent(in)::drctfile
  integer::i, sysstatus=0, system
  character(1000)::string, gcafile
  i=index(drctfile, ".")
  gcafile=drctfile(1:i-1) // ".gca"
  write (*, *) "...converting ", trim(drctfile), " to gca file."
  write (*, *)
  string=' xhmmstr_gamma model_R.hmm ' // trim(drctfile) // " " // trim(gcafile) // ' 0 ' // ' >> ' // ' LOG '
  sysstatus=system(string)
  if (sysstatus /= 0) stop "ERROR: System command for making gca files!!!"
  return
end subroutine drct2gca

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
  open(unit=15, file=filename, status="old", form="unformatted", access="direct", recl=recsize)
  i=0
  do
    i=i+1
    read(15, rec=i, iostat=status) seq_list, resSeq, structure_name, chainID, seq_resid, ss, rotamer,&
        (xyzca(j),j=1,3), phi_list, psi_list, omega_list, csum, (vall_pro(j), j=1, 20), gap_freq, ins_freq,&
        domainID, contex, (clusterID(j),j=1,3)
    if (status/=0) exit
    profile(i, 1:20)=vall_pro
  end do
  close(15)
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
  open (unit=14, file=filename, status="old", iostat=status)
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
      read (14, "(A)", iostat=status) line
      if (status/=0) exit
      if (line==" ") exit
      read (line, "(20F15.7)", iostat=status) (submatrix(x, y, w), y=1, 20)
      if (status/=0) stop "ERROR: Bad substitution matrix!!!"
    end do
    read (14, "(A)", iostat=status) line
    read (14, "(A)", iostat=status) line
    if (line==" ") w=w+1
    if (w>c) exit
  end do
  close (14, status="keep")
  return
end subroutine readmatrix
  
!-----subroutine readgamma-----!
subroutine readgamma(filename, length, gammamatrix)
  implicit none
  integer, intent(in)::length
  character(*), intent(in)::filename
  real, dimension(length, 0:281), intent(out)::gammamatrix
  integer::x, y, status=0
  character(3000)::gammaline
  gammamatrix=0
  gammaline=" "
  x=0
  open (unit=15, file=filename, status="old", iostat=status)
  if (status/=0) then
    write (*, "(A, $)") "ERROR: gcb file not found -"
    write (*, *) trim(filename)
    stop
  end if
  do
    read (15, "(A)", iostat=status) gammaline
    if (status/=0) exit
    if (gammaline(1:5)/="GAMMA") cycle
    x=x+1
    read (gammaline, "(7x, 282F9.6)", iostat=status) (gammamatrix(x, y), y=0, 281)
    if (status/=0) stop "ERROR: Bad gammaline in subroutine: readgamma!!!"
  end do
  close (15)
  return
end subroutine readgamma
  
!-----subroutine makeMij-----!
subroutine makeMij(len1, len2, sequence1, sequence2, gmatrix1, gmatrix2, m1, m2, Mij)
  implicit none
  integer, intent(in)::len1, len2
  character(*), intent(in)::sequence1, sequence2
  real, dimension(len1, 0:281), intent(inout)::gmatrix1
  real, dimension(len2, 0:281), intent(inout)::gmatrix2
  real, dimension(20,20,281), intent(inout)::m1
  real, dimension(20,281), intent(inout)::m2
  real, dimension(len1, len2), intent(out)::Mij
  real, dimension(20)::Pi
  integer::i, j, w, p1, p2
  real::sumM, sum1, sum2
  character(27)::aa="ACDEFGHIKLMNPQRSTVWYBZX-~.*"
  gmatrix1=gmatrix1+0.000000001
  gmatrix2=gmatrix2+0.000000001
  do w=1, 281
    sum1=sum(m1(:,:,w))
    if (sum1>0) then
      m1(:,:,w)=m1(:,:,w)/sum1
    end if
  end do
  sum2=sum(m2(:,:))
  m2=m2/sum2
  Pi=0
  do i=1, 20
    do w=1, 281
      Pi(i)=Pi(i)+m2(i,w)
    end do
  end do
  Mij=0
  do i=1, len1
    do j=1, len2
      p1=index(aa, sequence1(i:i))
      p2=index(aa, sequence2(j:j))
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
subroutine makedrctMij(len1, len2, pro1, pro2, gmatrix1, gmatrix2, m1, m2, Mij)
  implicit none
  integer, intent(in)::len1, len2
  real, dimension(len1, 20), intent(inout)::pro1
  real, dimension(len2, 20), intent(inout)::pro2
  real, dimension(len1, 0:281), intent(inout)::gmatrix1
  real, dimension(len2, 0:281), intent(inout)::gmatrix2
  real, dimension(20,20,281), intent(inout)::m1
  real, dimension(20,281), intent(inout)::m2
  real, dimension(len1, len2), intent(out)::Mij
  real, dimension(20)::Pi
  integer::i, j, aa1, aa2, w, p1, p2
  real::sumM, sum1, sum2, score
  character(27)::aa="ACDEFGHIKLMNPQRSTVWYBZX-~.*"
  gmatrix1=gmatrix1+0.000000001
  gmatrix2=gmatrix2+0.000000001
  do w=1, 281
    sum1=sum(m1(:,:,w))
    if (sum1>0) then
      m1(:,:,w)=m1(:,:,w)/sum1
    end if
  end do
  sum2=sum(m2(:,:))
  m2=m2/sum2
  Pi=0
  do i=1, 20
    do w=1, 281
      Pi(i)=Pi(i)+m2(i,w)
    end do
  end do
  Mij=0
  do i=1, len1
    do j=1, len2
      do aa1=1, 20
        do aa2=1, 20
          if (pro1(i, aa1)>0 .and. pro2(j, aa2)>0) then
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
            Mij(i, j)=Mij(i, j)+pro1(i, aa1)* pro2(j, aa2)*score
          end if
        end do
      end do
    end do
  end do
  return
end subroutine makedrctMij

!-----subroutine makeAij-----!
subroutine makeAij(len1, len2, sequence1, sequence2, Aij)
  implicit none
  integer, intent(in)::len1, len2
  character(*), intent(in)::sequence1, sequence2
  integer, dimension(len1, len2), intent(out)::Aij
  integer::i, j, p1, p2
  character(27)::aa="ACDEFGHIKLMNPQRSTVWYBZX-~.*"
  Aij=0
  do i=1, len1
    do j=1, len2
      p1=index(aa, sequence1(i:i))
      p2=index(aa, sequence2(j:j))
      if (p1>20 .OR. p2>20 .OR. p1<1 .OR. P2<1) cycle
      if (p1==p2) then
        Aij(i, j)=1
      else
        Aij(i, j)=0
      end if
    end do
  end do
  return
end subroutine makeAij
  
!-----subroutine random-----!
subroutine INITRAND()
  implicit none
  integer::i, j
  real::x
  integer, dimension(8)::dt
  i=2
  dt=1
  call date_and_time(values=dt)
  call random_seed(size=i,put=dt(5:8))
  call random_number(x)
end subroutine INITRAND

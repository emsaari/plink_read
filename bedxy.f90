program bedtimesy 
implicit none
integer                      :: n,m,n4,i,j,j2,ioerr
integer(kind=1), allocatable,dimension(:,:)    :: bed
integer(kind=1), allocatable,dimension(:,:)    :: Mat
real            , allocatable,dimension(:)     :: y,sol
integer(kind=8), allocatable,dimension(:)      :: ID
character(len=20)            :: bedfile,mapfile,famfile,outfile,infile,vecfile
integer            :: inunit=10, bedunit=11, mapunit=12, famunit=13,outunit=141
integer            :: vecunit=14


! name the files
call getarg(1,infile)
print *,"infile ",infile
call getarg(2,vecfile)
print *,"vecfile ",vecfile

bedfile=  infile
mapfile= trim(bedfile)//".bim"
famfile= trim(bedfile)//".fam"
outfile= trim(bedfile)//".txt"
bedfile= trim(bedfile)//".bed"
! count the rows to get parameters n,m
open(mapunit,file=mapfile)
call row_count(mapunit,m)
!rewind(mapunit)
open(famunit,file=famfile)
call row_count(famunit,n)
rewind(famunit)

open(bedunit,file=bedfile,form="unformatted",access="stream")
open(outunit, file=outfile,form="formatted")
open(vecunit, file=vecfile,form="formatted")

n4=ceiling(n/4.)
print *,mapfile,m
print *,famfile,n,"rows",n4

allocate(ID(n))
allocate(bed(n4,m))
allocate(Mat(1,1),y(m),sol(4*n4))


read(famunit,*)ID

call bimread(bedunit,bed,n4,m,Mat,0)
!
!   
do i=1,m
  read(vecunit,*,iostat=ioerr)y(i)
  if (ioerr/=0) then; print *,"error in read Mat ",ioerr; stop; endif
enddo

!do i=1,m
!    print '(5b8)',(bed(j,i),j=1,n4)
!enddo

call bimxy(bed,n4,m,y,sol,"n")
print *,
print *,"sol"
do i=1,4*n4
  print '(80000f6.2)',sol(i)
enddo
print *,"y"
call bimxy(bed,n4,m,y,sol,"t")
do i=1,m
  print '(80000f7.2)',y(i)
enddo


end

subroutine bimxy(bed,n4,m,snps,inds,info)
implicit none

character 			    :: info  ! n: Zd &  t: Z'*d
integer                             :: i,j,ij(4),m,n4
integer(kind=1)                     :: bed(n4,m),uneven(4),even(4)
real                                :: snps(m),inds(4*n4)
real(kind=4), dimension(4)          :: wrk
real(kind=4)                        :: nol,wsol,wsol4(4)

nol=0.
uneven = (/ 1, 3, 5, 7 /)
even   = (/ 0, 2, 4, 6 /)


if (info.ne."t") then
inds=0.
  do j=1,m
    ij=(/1,2,3,4/);
    wsol=snps(j)
    do i=1,n4
      wrk     =               merge(wsol,nol,btest(bed(i,j),uneven))
      inds(ij)= inds(ij)+ wrk+ merge(wsol,nol,btest(bed(i,j),  even))
     ij=ij+4 
    enddo
  enddo
else    ! transpose
  snps=0.
  do j=1,m
    ij=(/1,2,3,4/);
    do i=1,n4
      wsol4=inds(ij)
      snps(j)  = snps(j)+sum(wsol4,btest(bed(i,j),uneven))  &
                        + sum(wsol4,btest(bed(i,j),  even))
     ij=ij+4 
    enddo
  enddo

endif  ! info
!y=sol(1:4*n4)
end ! subroutine bimxy


subroutine bimread(bedunit,bed,n4,m,Mat,ninfo)
implicit none

integer                             :: i,j,ij(4),bedunit,m,ninfo,n4,ioerr
integer(kind=1), dimension(n4,m)    :: bed
integer(kind=1), dimension(3)       :: magic
integer(kind=1), dimension(4*n4,m)     :: Mat
integer(kind=1), dimension(4)       :: nol,one,two,wrk,uneven,even

nol=0; one=1; two=2
uneven = (/ 1, 3, 5, 7 /)
even   = (/ 0, 2, 4, 6 /)
Mat=0

read(bedunit,iostat=ioerr)magic,bed
if (ioerr/=0) then; print *,"error in read bed ",ioerr; stop; endif
close(bedunit)
!print *,"magics ",magic
if (ninfo.ne.0) then
  do j=1,m
    ij=(/1,2,3,4/);
    do i=1,n4
       wrk       = merge(one,nol,btest(bed(i,j),uneven))
       Mat(ij,j) = merge(two,wrk,btest(bed(i,j),  even))
     ij=ij+4 
    enddo
  enddo
endif  ! ninfo==0

end ! subroutine bimread

subroutine bimwrite(bedunit,bed,n4,m,Mat,n)
implicit none

integer				    :: i,j,ij,bedunit,m,n,n4,ioerr	
integer(kind=1), dimension(n4,m)    ::  bed
integer(kind=1), dimension(3)       :: magic
integer(kind=1), dimension(4*n4,m)     ::  Mat

magic=(/108, 27, 1   /)
bed=0

!CALL MVBITS(FROM, FROMPOS, LEN, TO, TOPOS)
do j=1,m
  where (Mat(:,j).gt.0) Mat(:,j)=Mat(:,j)+1
  ij=1
  do i=1,n4
    call mvbits(Mat(ij,j),  0,2,bed(i,j),0)
    call mvbits(Mat(ij+1,j),0,2,bed(i,j),2)
    call mvbits(Mat(ij+2,j),0,2,bed(i,j),4)
    call mvbits(Mat(ij+3,j),0,2,bed(i,j),6)
   ij=ij+4 
  enddo
enddo

write(bedunit,iostat=ioerr)magic,bed
if (ioerr/=0) then; print *,"error in write bed",ioerr; stop; endif
close(bedunit)
print *,"bed written"
end !subroutine bimwrite

program bed2M 
implicit none
integer                      :: n,m,n4,i,j,j2,ioerr
integer(kind=1), allocatable,dimension(:,:)    :: bed
integer(kind=1), allocatable,dimension(:,:)    :: Mat
integer(kind=8), allocatable,dimension(:)      :: ID
character(len=20)            :: bedfile,mapfile,famfile,outfile,infile
integer            :: inunit=10, bedunit=11, mapunit=12, famunit=13,outunit=141


! name the files
call getarg(1,infile)
print *,"infile ",infile

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

n4=ceiling(n/4.)
print *,mapfile,m
print *,famfile,n,"rows",n4

allocate(ID(n))
allocate(bed(n4,m))
allocate(Mat(4*n4,m))

read(famunit,*)ID

call bimread(bedunit,bed,n4,m,Mat,n)
!
!   To get the ID numbers:
do i=1,n
  write(outunit,'(i5,50000i2)',iostat=ioerr)ID(i),(Mat(i,j),j=1,m)
  if (ioerr/=0) then; print *,"error in read Mat ",ioerr; stop; endif
enddo

!do i=1,m
!    print '(5b8)',(bed(j,i),j=1,n4)
!enddo

print *,
print *,"Mat"
do i=1,n
  print '(i7,80000i2)',ID(i),(Mat(i,j),j=1,m)
enddo
!call bimwrite(outunit,bed,n4,m,Mat,n)

end

subroutine bimread(bedunit,bed,n4,m,Mat,n)
implicit none

integer                             :: i,j,ij(4),bedunit,m,n,n4,ioerr
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
do j=1,m
  ij=(/1,2,3,4/);
  do i=1,n4
   wrk       = merge(one,nol,btest(bed(i,j),uneven))
   Mat(ij,j) = merge(two,wrk,btest(bed(i,j),  even))
   ij=ij+4 
  enddo
enddo

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

!do j=1,m
!    print '(5(1x,b8))',(bed(i,j),i=1,n4)
!enddo
write(bedunit,iostat=ioerr)magic,bed
if (ioerr/=0) then; print *,"error in write bed",ioerr; stop; endif
close(bedunit)
print *,"bed written"
end !subroutine bimwrite

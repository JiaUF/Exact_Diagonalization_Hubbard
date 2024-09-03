module NumOfOrbitalAndElectrons
   integer*8 :: N3d, N, Ndim
   integer*8 :: nup, ndn, ntol
endmodule NumOfOrbitalAndElectrons

module ModelParas
   double precision :: U, t, tt, ttt, miu
   double precision :: phase(30,30)
   double precision :: disp(0:15), lsf(0:15)
endmodule ModelParas

module MPIParas
   integer  ::   comm, myid, nprocs
   integer  ::   source, dest, tag, ierr, rc
   integer*8 :: nloc, localstart, localend
endmodule MPIParas

module ConstantParas
   integer*8 :: HsizeEst=72000, SprSizeEst=72000*200
   integer*8 :: Hsize_3dEst=72000, SprSize_3dEst=72000*200
   integer*8 :: NNN=1000, niter_CL=100, niter_CFE = 100, niter_CG=100
   double precision :: pi=3.1415926535, sqrthalf = 0.707106781
endmodule ConstantParas

module BettsCluster
   character*1 :: Nmore
   integer*8 :: rpoint(0:15,2), runit
   integer*8 :: kpoint(0:15,2), kunit
   integer*8 :: site(0:15,4), gk, gk_opt
   integer*8 :: nQpoints, Kmap(1:16)
   double complex :: QPhase(0:15,0:8)
   double precision :: Qpoint(0:15,2)
!   integer*8 :: N
endmodule BettsCluster


Program Main
use NumOfOrbitalAndElectrons; use ModelParas; use BettsCluster
use ConstantParas; 
implicit none
integer*8 :: ksize_0, SprSize_0
integer*8 :: ii, jj, kk, mm, nn
integer*8 :: qq, qqk, mqq
integer*8, allocatable :: Hsp_0(:)
integer*8, external :: factorial, ksubtract, ksum
integer*8, allocatable :: IndexI_0(:), IndexJ_0(:)
double precision :: E_0, E_0_opt
double precision, allocatable :: sparseH_0(:)
double precision, allocatable :: H_0(:), E_all(:)
double complex :: z
character*1 :: onechar


open(unit=10101,file='input_SpinCrlt',Status='old');
read(10101,'(2I8)')   N3d, Ndim
read(10101,'(2I8)')   nup, ndn
close(10101)


write(*,*) ''
write(*,*) '     N3d     Ndim '
write(*,'(2I8)') N3d, Ndim
write(*,*) ''
write(*,*) '     nup     ndn'
write(*,'(2I8)') nup, ndn
write(*,*) ''
write(*,*) '******************************'
write(*,*) ''

N=N3d
ntol = nup+ndn

!============== The Step 0: preparation ==============

call SetModelParameters

allocate(Hsp_0(Hsize_3dEst))
allocate(H_0(Hsize_3dEst))
allocate(E_all(Hsize_3dEst))
allocate(IndexI_0(SprSize_3dEst))
allocate(IndexJ_0(SprSize_3dEst))
allocate(sparseH_0(SprSize_3dEst))
SprSize_0=SprSize_3dEst

E_0=0.0d0;
H_0=0.0d0;
E_all=0.0d0;
Hsp_0=0
IndexI_0=0
IndexJ_0=0
sparseH_0=0.0d0

!============== The Step 1: Diagonalize 3d system ==============

! 1 stands for (0,0) momentum

! | phi 0 > at k = 0
gk_opt = 1
E_0_opt = 100
do ii=1, nQpoints
   gk = Kmap(ii);
   call GenHsp_kspace(gk, Hsp_0, ksize_0)
   call GenMatrix_kspace(gk, Hsp_0, ksize_0, IndexI_0, IndexJ_0, sparseH_0, SprSize_0)
   !write(*,*) 'Sparse Matrix print'
   !do jj=1,SprSize_0
   !   write(*,*) ' ', IndexI_0(jj), ' ', IndexJ_0(jj), ' ', sparseH_0(jj)
   !enddo
   write(*,*) 'momentum point = ', gk
   write(*,*) 'ksize_0 = ', ksize_0
   write(*,*) 'SprSize_0 = ', SprSize_0
   call MatrixDiagonalization(ksize_0, Hsp_0, SprSize_0, IndexI_0, IndexJ_0, sparseH_0, E_all, H_0)
   do jj=1, ksize_0
      write(*,*) '           E = ', E_all(jj)
   enddo
   !if(E_0 < E_0_opt) then
   !   E_0_opt = E_0
   !   gk_opt = gk
   !endif

enddo

end program







  

subroutine MatrixDiagonalization(Hsize_3d, Hsp_0, SprSize_3d, IndexI_3d, &
                  IndexJ_3d, sparseH_3d, E_0_3d, H_0_3d)
implicit none
INTEGER*8, INTENT(IN) :: Hsize_3d, SprSize_3d
INTEGER*8, DIMENSION(Hsize_3d), INTENT(IN) :: Hsp_0
INTEGER*8, DIMENSION(SprSize_3d), INTENT(IN) :: IndexI_3d, IndexJ_3d
DOUBLE precision, DIMENSION(SprSize_3d), INTENT(IN) :: sparseH_3d
double precision, DIMENSION(Hsize_3d), INTENT(OUT) :: E_0_3d
DOUBLE precision, DIMENSION(Hsize_3d), INTENT(OUT) :: H_0_3d


integer*8 :: info, ii, jj, kk
double precision :: tsum
double precision, allocatable :: WORK(:),Ht(:,:),W_ED(:),Ht0(:,:), Ax(:)

allocate(W_ED(Hsize_3d), Ax(Hsize_3d))
allocate(Ht(Hsize_3d, Hsize_3d), Ht0(Hsize_3d, Hsize_3d))
allocate(WORK(5*Hsize_3d))

Ht=0.0d0
Ht0=0.0d0
!write(*,*) 'start assign'
do ii = 1, SprSize_3d
   Ht(IndexI_3d(ii), IndexJ_3d(ii)) = sparseH_3d(ii)
   Ht0(IndexI_3d(ii), IndexJ_3d(ii)) = sparseH_3d(ii)
enddo

!do ii=1,Hsize_3d
!   do jj=1, Hsize_3d
!      write(*,'(2i10,2F10.2)') ii, jj, Ht(ii, jj), Ht(jj,ii) 
!   enddo
!enddo

!write(*,*) 'assign OK'

If (maxval(transpose(Ht)-Ht).ne.0.0d0.and.minval(transpose(Ht)-Ht).ne.0.0d0) then
write(*,*) "The Matrix is not Hermitian!"
!write(*,*) "The maxminum values of differentce", maxval(transpose(Ht)-Ht)
!write(*,*) "The minimum values of differentce", minval(transpose(Ht)-Ht)
else
write(*,*) "Matrix Hermitian check ok."
end if

info=0
call DSYEV('V','L',Hsize_3d,Ht,Hsize_3d,W_ED,WORK,5*Hsize_3d,info)
write(*,*) info
if (info.eq.0) write(*,*) 'LAPACK routine ok.'


do ii = 1,Hsize_3d
   ! test Ax - Ex = ?
   do kk=1, Hsize_3d
      Ax(kk) = 0
      do jj=1, Hsize_3d
         Ax(kk) = Ax(kk) + Ht0(kk,jj) * Ht(jj,ii)
      enddo 
      Ax(kk) = Ax(kk) - W_Ed(ii) * Ht(kk,ii)    
   enddo
   write(*,*) ii, maxval(Ax), minval(Ax)
enddo
write(*,*) 'Eigenvalue/eigenvector check finish'

do ii = 1, Hsize_3d
   H_0_3d(ii) = Ht(ii,1)
   E_0_3d(ii) = W_ED(ii)
enddo

deallocate(Ht, Ht0, W_ED, Ax, WORK)

end subroutine


subroutine MatrixDiagonalization_complex(Hsize_3d, SprSize_3d, IndexI_3d, &
                  IndexJ_3d, sparseH_3d, E_0_3d, H_0_3d)
implicit none
INTEGER*8, INTENT(IN) :: Hsize_3d, SprSize_3d
INTEGER*8, DIMENSION(SprSize_3d), INTENT(IN) :: IndexI_3d, IndexJ_3d
DOUBLE COMPLEX, DIMENSION(SprSize_3d), INTENT(IN) :: sparseH_3d
double precision, INTENT(OUT) :: E_0_3d
DOUBLE COMPLEX, DIMENSION(Hsize_3d), INTENT(OUT) :: H_0_3d


integer*8 :: LWORK
integer*8 :: info, ii
DOUBLE PRECISION, allocatable ::  RWORK(:)
double complex, allocatable :: WORK(:)
double complex, allocatable :: Ht(:,:)

allocate(Ht(Hsize_3d, Hsize_3d))
Ht=cmplx(0.0d0,0.0d0)

do ii = 1, SprSize_3d
   Ht(IndexI_3d(ii), IndexJ_3d(ii)) = sparseH_3d(ii)
enddo

LWORK=3*Hsize_3d
allocate(WORK(LWORK))
allocate(RWORK(LWORK))
info=0
!call zheev('V', 'U', Hsize_3d, Ht, Hsize_3d, &
!           E_0_3d, WORK, 3*Hsize_3d, RWORK, info)
call ZHEEV('V', 'U', Hsize_3d, Ht, Hsize_3d, &
           E_0_3d, WORK, 3*Hsize_3d, RWORK, info)

if (info.eq.0) then
   write(*,*) 'ZHEEV exit successfully.'
else
   write(*,*) '#info.neq.0# info = ', info
endif

H_0_3d(1:Hsize_3d) = Ht(1:Hsize_3d,1)

deallocate(WORK, RWORK, Ht)

end subroutine




subroutine GenHsp_kspace(kfixed, listqpt, ksize)
use ConstantParas; use NumOfOrbitalAndElectrons
use BettsCluster
implicit none

INTEGER*8, INTENT(IN):: kfixed
INTEGER*8, DIMENSION(Hsize_3dEst), INTENT(OUT):: listqpt
INTEGER*8, INTENT(OUT):: ksize

integer*8 :: iiup,iidn, Hsizet
integer*8 :: listup(65536,2),listdn(65536,2),listversup(65536),listversdn(65536)
integer*8 :: ktempx,ktempy, px, py
integer*8 :: jj, ii, i, kstat, iq
integer*8 :: iktemp1,iktemp2,iktemp
integer*8 :: kstat1,kstat2
integer*8 :: temp,tempi
integer*8 :: iup,idn, jup, jdn
integer*8 :: initstatup, initstatdn
integer*8, external :: sumeverybit, factorial
integer*8 :: ktempx1,ktempx2,ktempy1,ktempy2



iiup=0;iidn=0
do ii=0,2**N-1
        listup(ii+1,1)=0;listup(ii+1,2)=0;listdn(ii+1,1)=0;listdn(ii+1,2)=0;
        listversup(ii+1)=0;listversdn(ii+1)=0
end do

ktempx=0
ktempy=0
iktemp=0
kstat=1
do iup=0,2**N-1
        if(sumeverybit(iup).eq.nup) then                !find one up state
                do i=0,N-1
                   if(BTEST(iup,i)) then
                      ktempx = ktempx + kpoint(i,1)
                      ktempy = ktempy + kpoint(i,2)
                   endif
                enddo
                ktempx = mod(ktempx,kunit)
                ktempy = mod(ktempy,kunit)
                do while (kstat.ne.0)
                   if((kpoint(iktemp,1).eq.ktempx).and.(kpoint(iktemp,2).eq.ktempy)) then
                      kstat=0
                   else
                      iktemp=iktemp+1
                   endif
                   if(iktemp.eq.N) then
                      write(*,*) 'iktemp out of bounds, spin up'
                      write(*,*) iup
                      stop
                   endif
                enddo
                listup(iiup+1,1)=iup;
                listup(iiup+1,2)=iktemp;
                listversup(iup)=iiup;
                iiup=iiup+1;
                ktempx=0;
                ktempy=0;
                iktemp=0;
                kstat=1;

        endif
enddo

!********************************************************************************
! This loop is superfluous at Sz=0
!********************************************************************************

ktempx=0
ktempy=0
iktemp=0
kstat=1

do idn=0,2**N-1
        if(sumeverybit(idn).eq.ndn) then                !find one down state
                do i=0,N-1
                   if(BTEST(idn,i)) then
                      ktempx = ktempx + kpoint(i,1)
                      ktempy = ktempy + kpoint(i,2)
                   endif
                enddo
                ktempx = mod(ktempx,kunit)
                ktempy = mod(ktempy,kunit)
                do while (kstat.ne.0)
                   if(kpoint(iktemp,1).eq.ktempx.and.kpoint(iktemp,2).eq.ktempy) then
                      kstat=0
                   else
                      iktemp=iktemp+1
                   endif
                   if(iktemp.eq.N) then
                      write(*,*) 'iktemp out of bounds, spin down'
                      write(*,*) idn
                      stop
                   endif
                enddo
                listdn(iidn+1,1)=idn;
                listdn(iidn+1,2)=iktemp;
                listversdn(idn)=iidn;
                iidn=iidn+1;
                ktempx=0;
                ktempy=0;
                iktemp=0;
                kstat=1;
        endif
enddo

Hsizet = (factorial(N)/factorial(nup)/factorial(N-nup))*&
(factorial(N)/factorial(ndn)/factorial(N-ndn))

ktempx=0
ktempy=0
iktemp=0
kstat=1
iq=0
listqpt=0
do jj=0,Hsizet-1         !loop B
        jdn=mod(jj,iidn);
        jup=jj/iidn;
        initstatup=listup(jup+1,1);
        initstatdn=listdn(jdn+1,1);
        ktempx=mod(kpoint(listup(jup+1,2),1)+kpoint(listdn(jdn+1,2),1),kunit)
        ktempy=mod(kpoint(listup(jup+1,2),2)+kpoint(listdn(jdn+1,2),2),kunit)
        do while (kstat.ne.0)
           if(kpoint(iktemp,1).eq.ktempx.and.kpoint(iktemp,2).eq.ktempy) then
              kstat=0
           else
              iktemp=iktemp+1
           endif
           if(iktemp.eq.N) then
              write(*,*) 'iktemp out of bounds, column'
              write(*,*) jj
              stop
           endif
        enddo
        if(iktemp.eq.kfixed) then
           iq=iq+1
           listqpt(iq) = initstatup*(2**N)+initstatdn
           write(*,'(I10, B20)') iq, listqpt(iq)
        endif
        iktemp=0
        kstat=1
enddo !loop B

!write(*,*) 'listup&dn set up OK', myid
ksize = iq

write(*,*) 'Inside GenHsp kfixed = ', kfixed
write(*,*) 'px and py = ', kpoint(kfixed, 1), kpoint(kfixed, 2)

end subroutine




subroutine GenMatrix_kspace(kfixed, listqpt, ksize, IndexIt, IndexJt, sparseHt, SprSizet)
use NumOfOrbitalAndElectrons; use ModelParas; use MPIParas; use ConstantParas
use BettsCluster
implicit none

integer*8, INTENT(IN) :: kfixed, ksize
integer*8, DIMENSION(ksize), INTENT(IN) :: listqpt
integer*8, INTENT(INOUT) :: SprSizet
integer*8, DIMENSION(SprSizet), INTENT(OUT) :: IndexIt, IndexJt
DOUBLE PRECISION, DIMENSION(SprSizet), INTENT(OUT) :: sparseHt

integer*8 :: H_index(1:3000)
integer*8 :: jj, kk, l, kkp, qq, ii, pp
integer*8 :: temp, tempi
integer*8 :: initstatup,initstatdn
integer*8 :: afterstatup,afterstatdn,initsign,aftersign
integer*8 :: midstatup,midstatdn, midsign
integer*8 :: iktemp1,iktemp2,iktemp
integer*8 :: kstat1,kstat2
integer*8 :: ktempx,ktempy
integer*8 :: ktempx1,ktempx2,ktempy1,ktempy2
integer*8, external :: sumeverybit, sumbeforebit
integer*8, allocatable :: H(:)
double precision :: H_value(1:3000)


!!!******************Change the parameters here****************
SprSizet = ksize*200

!if(mod(ksize,nprocs).eq.0) then
!   nloc = (ksize / nprocs)
!else
!   nloc = (ksize / nprocs)+1
!endif

allocate(H(ksize))

H=0
H_value=0
H_index=0
tempi=0

do jj= 1, ksize  !loop B

        temp=0
        initstatup=listqpt(jj)/2**N;
        initstatdn=mod(listqpt(jj),2**N);

        do kk=0,N-1  !loop A

           !**************************************************************************
           ! KE
           !**************************************************************************

           if(BTEST(initstatup,kk)) then
              l=jj
              if(H(l).eq.0) then
                 temp=temp+1
                 H(l)=temp
                 H_value(temp)=disp(kk)
                 H_index(temp)=l
              else
                 H_value(H(l))=H_value(H(l))+disp(kk)
              endif
           endif
           if(BTEST(initstatdn,kk)) then
              l=jj
              if(H(l).eq.0) then
                 temp=temp+1
                 H(l)=temp
                 H_value(temp)=disp(kk)
                 H_index(temp)=l
              else
                 H_value(H(l))=H_value(H(l))+disp(kk)
              endif
           endif


           !******************************************************************************
           !               \---
           !                \   +        +
           ! Hubbard U:  U  /  c        c       c      c         (from Richard Scalletar's note)
           !               /--- k'+q,up  k-q,dn  k',up  k,dn
           !               q,k,k'
           !******************************************************************************
           ! https://www.cs.ucdavis.edu/~bai/QUEST/tutorial/hubbard7.pdf

           ktempx1=0
           ktempy1=0
           iktemp1=0
           kstat1=1
           ktempx2=0
           ktempy2=0
           iktemp2=0
           kstat2=1
           !write(*,*) 'kunit = ', kunit
           !write(*,'(2B6)') initstatup, initstatdn
           do kkp=0,N-1
              do qq=0,N-1
                 if(BTEST(initstatdn,kk).and.BTEST(initstatup,kkp)) then
                    ktempx1 = mod(kpoint(kk,1) - kpoint(qq,1)+kunit,kunit)
                    ktempx2 = mod(kpoint(kkp,1) + kpoint(qq,1)+kunit,kunit)
                    ktempy1 = mod(kpoint(kk,2) - kpoint(qq,2)+kunit,kunit)
                    ktempy2 = mod(kpoint(kkp,2) + kpoint(qq,2)+kunit,kunit)
                    do while (kstat1.ne.0)
                       if(kpoint(iktemp1,1).eq.ktempx1.and.kpoint(iktemp1,2).eq.ktempy1) then
                         kstat1=0
                       else
                         iktemp1=iktemp1+1
                       endif
                       if(iktemp1.eq.N) then
                          write(*,*) 'iktemp out of bounds, row k-q'
                          write(*,*) kkp, qq
                          stop
                       endif
                    enddo
                    do while (kstat2.ne.0)
                       if(kpoint(iktemp2,1).eq.ktempx2.and.kpoint(iktemp2,2).eq.ktempy2) then
                         kstat2=0
                       else
                         iktemp2=iktemp2+1
                       endif
                       if(iktemp2.eq.N) then
                          write(*,*) 'iktemp out of bounds, row kp+q'
                          write(*,*) jj
                          stop
                       endif
                    enddo
                    !write(*,*) '   kk,  kp,   qq,  k-q,  kp+q'
                    !write(*,'(5I6)') kk, kkp, qq, iktemp1, iktemp2
                    !set iktemp1 to k-q
                    !set iktemp2 to kp+q
                    afterstatdn=IBCLR(initstatdn,kk)
                    aftersign=((-1)**sumbeforebit(initstatdn,kk))
                    afterstatup=IBCLR(initstatup,kkp)
                    aftersign=aftersign*((-1)**(sumbeforebit(afterstatup,kkp)+sumeverybit(afterstatdn)))
                    if(.not.BTEST(afterstatdn,iktemp1).and..not.BTEST(afterstatup,iktemp2)) then
                       afterstatdn=IBSET(afterstatdn,iktemp1)
                       aftersign=aftersign*((-1)**sumbeforebit(afterstatdn,iktemp1))
                       afterstatup=IBSET(afterstatup,iktemp2)
                       aftersign=aftersign*((-1)**(sumbeforebit(afterstatup,iktemp2)+sumeverybit(afterstatdn)))
                       call BinarySearch(listqpt,ksize,afterstatup*(2**N)+afterstatdn,l)
                       if(l.eq.-1) then
                          write(*,*) 'binary search out of range'
                          !write(*,'(2B6)') afterstatup, afterstatdn
                          stop
                       endif

                       !write(*,'(2I4, 2B10, F10.2)') jj, l, listqpt(jj), afterstatup*(2**N)+afterstatdn, -aftersign*U/dfloat(N)

                       if(H(l).eq.0) then
                          temp=temp+1
                          H(l)=temp
                          H_value(temp)=-aftersign*U/dfloat(N)
                          H_index(temp)=l
                       else
                          H_value(H(l))=H_value(H(l))-aftersign*U/dfloat(N)
                       endif

                    endif
                    iktemp1=0
                    kstat1=1
                    iktemp2=0
                    kstat2=1
                 endif
              enddo
           enddo

        enddo !loop A

        do ii=1,temp
           if(H_value(ii) > 0.00001 .or. H_value(ii) < -0.00001) then
              tempi = tempi + 1
              IndexJt(tempi)  = H_index(ii) !still on the old numbering
              IndexIt(tempi)  = jj
              sparseHt(tempi) = H_value(ii)
           endif
           H(H_index(ii))=0
        enddo

enddo   !loop B

SprSizet = tempi

!write(*,*) 'spr_size in GenMatrix subroutine =',SprSizet

end subroutine





subroutine SetModelParameters
use ModelParas; use NumOfOrbitalAndElectrons; use BettsCluster
use ConstantParas
implicit none

integer*8 :: ii, kk
t=1.00d0; tt=-0.30d0; ttt=0.0d0
miu=0.0d0
U=8.0d0;

!THe hoppping:
phase(:,:)=0.0d0;

selectcase(N3d)
case(4)

 if(Ndim.eq.2) then

site(0,1)=2;site(0,2)=2;site(0,3)=1;site(0,4)=1;
site(1,1)=3;site(1,2)=3;site(1,3)=0;site(1,4)=0;
site(2,1)=0;site(2,2)=0;site(2,3)=3;site(2,4)=3;
site(3,1)=1;site(3,2)=1;site(3,3)=2;site(3,4)=2;

  else if(Ndim.eq.1) then

site(0,1)=1;site(0,2)=3;site(0,3)=0;site(0,4)=0;
site(1,1)=2;site(1,2)=0;site(1,3)=1;site(1,4)=1;
site(2,1)=3;site(2,2)=1;site(2,3)=2;site(2,4)=2;
site(3,1)=0;site(3,2)=2;site(3,3)=3;site(3,4)=3;

  endif



case(8)
!************************************************************
! K-points for the 8-site cluster with (+) integer values
! with a factor of (pi/2) removed
!
!                 |     x     |     x
!                 |      (1,3)|      (3,3)
!                 |           |
!                 x-----------x------
!                 |(0,2)      |(2,2)
!                 |           |
!                 |     x     |     x
!                 |      (1,1)|      (3,1)
!                 |           |
!                 x-----------x------
!                 |(0,0)      |(2,0)
!************************************************************
rpoint(0,1)=0;rpoint(0,2)=0;
rpoint(1,1)=2;rpoint(1,2)=0;
rpoint(2,1)=0;rpoint(2,2)=2;
rpoint(3,1)=2;rpoint(3,2)=2;
rpoint(4,1)=1;rpoint(4,2)=1;
rpoint(5,1)=3;rpoint(5,2)=1;
rpoint(6,1)=1;rpoint(6,2)=3;
rpoint(7,1)=3;rpoint(7,2)=3;

runit=4

! site( ,1) up right
! site( ,2) down left
! site( ,3) left up
! site( ,4) right down

if(Ndim.eq.2) then
site(0,1)=2;site(0,2)=7;site(0,3)=3;site(0,4)=6;
site(1,1)=3;site(1,2)=6;site(1,3)=2;site(1,4)=7;
site(2,1)=5;site(2,2)=0;site(2,3)=4;site(2,4)=1;
site(3,1)=4;site(3,2)=1;site(3,3)=5;site(3,4)=0;
site(4,1)=6;site(4,2)=3;site(4,3)=7;site(4,4)=2;
site(5,1)=7;site(5,2)=2;site(5,3)=6;site(5,4)=3;
site(6,1)=1;site(6,2)=4;site(6,3)=0;site(6,4)=5;
site(7,1)=0;site(7,2)=5;site(7,3)=1;site(7,4)=4;

else if(Ndim.eq.1) then
site(0,1)=1;site(0,2)=7;site(0,3)=0;site(0,4)=0;
site(1,1)=2;site(1,2)=0;site(1,3)=1;site(1,4)=1;
site(2,1)=3;site(2,2)=1;site(2,3)=2;site(2,4)=2;
site(3,1)=4;site(3,2)=2;site(3,3)=3;site(3,4)=3;
site(4,1)=5;site(4,2)=3;site(4,3)=4;site(4,4)=4;
site(5,1)=6;site(5,2)=4;site(5,3)=5;site(5,4)=5;
site(6,1)=7;site(6,2)=5;site(6,3)=6;site(6,4)=6;
site(7,1)=0;site(7,2)=6;site(7,3)=7;site(7,4)=7;
endif


case(10)

  if(Ndim.eq.1) then
site(0,1)=1;site(0,2)=9;site(0,3)=0;site(0,4)=0;
site(1,1)=2;site(1,2)=0;site(1,3)=1;site(1,4)=1;
site(2,1)=3;site(2,2)=1;site(2,3)=2;site(2,4)=2;
site(3,1)=4;site(3,2)=2;site(3,3)=3;site(3,4)=3;
site(4,1)=5;site(4,2)=3;site(4,3)=4;site(4,4)=4;
site(5,1)=6;site(5,2)=4;site(5,3)=5;site(5,4)=5;
site(6,1)=7;site(6,2)=5;site(6,3)=6;site(6,4)=6;
site(7,1)=8;site(7,2)=6;site(7,3)=7;site(7,4)=7;
site(8,1)=9;site(8,2)=7;site(8,3)=8;site(8,4)=8;
site(9,1)=0;site(9,2)=8;site(9,3)=9;site(9,4)=9;
  else
write(*,*) 'input Betts cluster dimension is wrong'
  endif


case(12)

  if(Ndim.eq.1) then
site(0,1)=1;site(0,2)=11;site(0,3)=0;site(0,4)=0;
site(1,1)=2;site(1,2)=0;site(1,3)=1;site(1,4)=1;
site(2,1)=3;site(2,2)=1;site(2,3)=2;site(2,4)=2;
site(3,1)=4;site(3,2)=2;site(3,3)=3;site(3,4)=3;
site(4,1)=5;site(4,2)=3;site(4,3)=4;site(4,4)=4;
site(5,1)=6;site(5,2)=4;site(5,3)=5;site(5,4)=5;
site(6,1)=7;site(6,2)=5;site(6,3)=6;site(6,4)=6;
site(7,1)=8;site(7,2)=6;site(7,3)=7;site(7,4)=7;
site(8,1)=9;site(8,2)=7;site(8,3)=8;site(8,4)=8;
site(9,1)=10;site(9,2)=8;site(9,3)=9;site(9,4)=9;
site(10,1)=11;site(10,2)=9;site(10,3)=10;site(10,4)=10;
site(11,1)=0;site(11,2)=10;site(11,3)=11;site(11,4)=11;
  else
write(*,*) 'input Betts cluster dimension is wrong'
  endif

endselect




selectcase(N3d)
case(4)
if(Ndim.eq.1) then
kpoint(0,1)=0;kpoint(0,2)=0;
kpoint(1,1)=1;kpoint(1,2)=0;
kpoint(2,1)=2;kpoint(2,2)=0;
kpoint(3,1)=3;kpoint(3,2)=0;

kunit=4

nQpoints=3
Kmap(1) = 0
Kmap(2) = 1
Kmap(3) = 2
endif

case(8)
if(Ndim.eq.2) then
!************************************************************
! K-points for the 8-site cluster with (+) integer values
! with a factor of (pi/2) removed
!
!                 |     x     |     x
!                 |      (1,3)|      (3,3)
!                 |           |
!                 x-----------x------
!                 |(0,2)      |(2,2)
!                 |           |
!                 |     x     |     x
!                 |      (1,1)|      (3,1)
!                 |           |
!                 x-----------x------
!                 |(0,0)      |(2,0)
!************************************************************
kpoint(0,1)=0;kpoint(0,2)=0;
kpoint(1,1)=2;kpoint(1,2)=0;
kpoint(2,1)=0;kpoint(2,2)=2;
kpoint(3,1)=2;kpoint(3,2)=2;
kpoint(4,1)=1;kpoint(4,2)=1;
kpoint(5,1)=3;kpoint(5,2)=1;
kpoint(6,1)=1;kpoint(6,2)=3;
kpoint(7,1)=3;kpoint(7,2)=3;

kunit=4

nQpoints=4
Kmap(1) = 0
Kmap(2) = 1
Kmap(3) = 3
Kmap(4) = 4

else if(Ndim.eq.1) then !one dimensional lattice

kpoint(0,1)=0;kpoint(0,2)=0;
kpoint(1,1)=1;kpoint(1,2)=0;
kpoint(2,1)=2;kpoint(2,2)=0;
kpoint(3,1)=3;kpoint(3,2)=0;
kpoint(4,1)=4;kpoint(4,2)=0;
kpoint(5,1)=5;kpoint(5,2)=0;
kpoint(6,1)=6;kpoint(6,2)=0;
kpoint(7,1)=7;kpoint(7,2)=0;

kunit=8

nQpoints=5
Kmap(1) = 0
Kmap(2) = 1
Kmap(3) = 2
Kmap(4) = 3
Kmap(5) = 4

endif

case(10)

if(Ndim.eq.2) then
!************************************************************
! K-points for the 16B-site cluster with (+) integer values
! with a factor of (pi/2) removed
!
!                 p-----x-----x-----x-----x-----p
!                 |                    b  |
!                 |                       |
!                 x     x     p     x     x     x
!                 |  b                    |
!                 |                       |
!                 x     x     x     x     p     x
!                 |              b        |
!                 |                       |
!                 x     p     x     x     x     x
!                 |                       |  b
!                 |                       | 
!                 x     x     x     p     x     x
!                 |        b              |
!                 |                       |
!                 p-----x-----x-----x-----x-----x

!************************************************************

kpoint(0,1)=0;kpoint(0,2)=0;
kpoint(1,1)=3;kpoint(1,2)=1;
kpoint(2,1)=2;kpoint(2,2)=4;
kpoint(3,1)=5;kpoint(3,2)=5;
kpoint(4,1)=6;kpoint(4,2)=2;
kpoint(5,1)=9;kpoint(5,2)=3;
kpoint(6,1)=8;kpoint(6,2)=6;
kpoint(7,1)=1;kpoint(7,2)=7;
kpoint(8,1)=4;kpoint(8,2)=8;
kpoint(9,1)=7;kpoint(9,2)=9;

kunit=10

nQpoints=4
Kmap(1) = 0
Kmap(2) = 1
Kmap(3) = 2
Kmap(4) = 3

else if(Ndim.eq.1) then !one dimensional lattice

kpoint(0,1)=0;kpoint(0,2)=0;
kpoint(1,1)=1;kpoint(1,2)=0;
kpoint(2,1)=2;kpoint(2,2)=0;
kpoint(3,1)=3;kpoint(3,2)=0;
kpoint(4,1)=4;kpoint(4,2)=0;
kpoint(5,1)=5;kpoint(5,2)=0;
kpoint(6,1)=6;kpoint(6,2)=0;
kpoint(7,1)=7;kpoint(7,2)=0;
kpoint(8,1)=8;kpoint(8,2)=0;
kpoint(9,1)=9;kpoint(9,2)=0;

kunit=10

nQpoints=6
Kmap(1) = 0
Kmap(2) = 1
Kmap(3) = 2
Kmap(4) = 3
Kmap(5) = 4
Kmap(6) = 5

endif


case(12)
!************************************************************
! K-points for the 12D-site cluster with (+) integer values
! with a factor of (pi/2) removed
!
!                 x-----x-----x-----x-----x-----x
!                 |                 |      
!                 |                 |      
!                 x     x     x     x     x     x
!                 |                 |      
!                 |                 |      
!                 B     x     x     x     x     x
!                 |        _        |      
!                 |                 |      
!                 x-----x-----x-----B-----x-----x
!                 |                 |      
!                 B                 |      
!                 x     x     x     x     x     x
!                 |                 |        
!                 |        P        |       
!                 x     x     x     x     x     x
!                 |                 B      
!                 |                 |      
!                 B-----x-----x-----x-----x-----x

!************************************************************

kpoint(0,1)=0;kpoint(0,2)=0;
kpoint(1,1)=0;kpoint(1,2)=4;
kpoint(2,1)=3;kpoint(2,2)=3;
kpoint(3,1)=6;kpoint(3,2)=2;
kpoint(4,1)=9;kpoint(4,2)=1;
kpoint(5,1)=0;kpoint(5,2)=8;
kpoint(6,1)=3;kpoint(6,2)=7;
kpoint(7,1)=6;kpoint(7,2)=6;
kpoint(8,1)=9;kpoint(8,2)=5;
kpoint(9,1)=3;kpoint(9,2)=11;
kpoint(10,1)=6;kpoint(10,2)=10;
kpoint(11,1)=9;kpoint(11,2)=9;

kunit=12

nQpoints=5
Kmap(1) = 0
Kmap(2) = 1
Kmap(3) = 2
Kmap(4) = 3
Kmap(5) = 7



case(16)
   selectcase(Nmore)
   case('B')
!************************************************************
! K-points for the 16B-site cluster with (+) integer values
! with a factor of (pi/2) removed
!
!                 x     x     x     x     x
!                 |                       |
!                 |                       |
!                 x-----x-----x-----x-----x-----x
!                 |                       |
!                 |(0,3)(1,3) (2,3) (3,3) |
!                 x     x     x     x     x     x
!                 |                       |
!                 |(0,2)(1,2) (2,2) (3,2) |
!                 x     x     x     x     x     x
!                 |                       |
!                 |(0,1)(1,1) (2,1) (3,1) |
!                 x     x     x     x     x     x
!                 |                       |
!                 |(0,0)(1,0) (2,0) (3,0) |
!                 x-----x-----x-----x-----x-----x
!                 |                       |
!************************************************************
kpoint(0,1)=0;kpoint(0,2)=0;
kpoint(1,1)=1;kpoint(1,2)=0;
kpoint(2,1)=2;kpoint(2,2)=0;
kpoint(3,1)=3;kpoint(3,2)=0;
kpoint(4,1)=0;kpoint(4,2)=1;
kpoint(5,1)=1;kpoint(5,2)=1;
kpoint(6,1)=2;kpoint(6,2)=1;
kpoint(7,1)=3;kpoint(7,2)=1;
kpoint(8,1)=0;kpoint(8,2)=2;
kpoint(9,1)=1;kpoint(9,2)=2;
kpoint(10,1)=2;kpoint(10,2)=2;
kpoint(11,1)=3;kpoint(11,2)=2;
kpoint(12,1)=0;kpoint(12,2)=3;
kpoint(13,1)=1;kpoint(13,2)=3;
kpoint(14,1)=2;kpoint(14,2)=3;
kpoint(15,1)=3;kpoint(15,2)=3;

kunit=4



   case('A')
!************************************************************
! K-points for the 16A-site cluster with (+) integer values
! with a factor of (pi/4) removed
!
!                 |  x     x     x     x     x
!                 |                       |
!                 |                       |
!                 x-----x-----x-----x-----x-----x
!                 |                       |
!                 |  (1,6)(3,6)(5,6)(7,6) |
!                 |  x     x     x     x  |  x
!                 |                       |
!                 |(0,4)(2,4)(4,4)(6,4)   |
!                 x     x     x     x     x     x
!                 |                       |
!                 |(1,2)(3,2) (5,2) (7,2) |
!                 |  x     x     x     x  |  x
!                 |                       |
!                 |(0,0)(2,0) (4,0) (6,0) |
!                 x-----x-----x-----x-----x-----x
!                 |                       |
!************************************************************
kpoint(0,1)=0;   kpoint(0,2)=0;
kpoint(1,1)=2;   kpoint(1,2)=0;
kpoint(2,1)=4;   kpoint(2,2)=0;
kpoint(3,1)=6;   kpoint(3,2)=0;
kpoint(4,1)=1;   kpoint(4,2)=2;
kpoint(5,1)=3;   kpoint(5,2)=2;
kpoint(6,1)=5;   kpoint(6,2)=2;
kpoint(7,1)=7;   kpoint(7,2)=2;
kpoint(8,1)=0;   kpoint(8,2)=4;
kpoint(9,1)=2;   kpoint(9,2)=4;
kpoint(10,1)=4;  kpoint(10,2)=4;
kpoint(11,1)=6;  kpoint(11,2)=4;
kpoint(12,1)=1;  kpoint(12,2)=6;
kpoint(13,1)=3;  kpoint(13,2)=6;
kpoint(14,1)=5;  kpoint(14,2)=6;
kpoint(15,1)=7;  kpoint(15,2)=6;
!**********************************************************************
kunit=8

   endselect
endselect


do kk=0,N3d-1
   disp(kk) = -2*t*(cos(real(kpoint(kk,1))*2.0d0*pi/kunit)+cos(real(kpoint(kk,2))*2.0d0*pi/kunit))&
                -4*tt*(cos(real(kpoint(kk,1))*2.0d0*pi/kunit)*cos(real(kpoint(kk,2))*2.0d0*pi/kunit))&
                  -2*ttt*(cos(2*real(kpoint(kk,1))*2.0d0*pi/kunit)+cos(2*real(kpoint(kk,2))*2.0d0*pi/kunit))&
                    -miu
   lsf(kk) = (cos(real(kpoint(kk,1))*2.0d0*pi/kunit)+cos(real(kpoint(kk,2))*2.0d0*pi/kunit))*0.5&
             -0.3*cos(real(kpoint(kk,1))*2.0d0*pi/kunit)*cos(real(kpoint(kk,2))*2.0d0*pi/kunit); ! lattice structure factor
   write(*,*) 'kk, cos', kk, cos(real(kpoint(kk,1))*2.0d0*pi/kunit)+cos(real(kpoint(kk,2))*2.0d0*pi/kunit)
   write(*,*) 'kk,disp', kk, disp(kk)
enddo

end subroutine

!*****************************************************


subroutine BinarySearch(listqpt,ksize,statupdn,l)
implicit none
INTEGER*8, INTENT(IN) :: ksize
INTEGER*8, DIMENSION(1:ksize), INTENT(IN) :: listqpt
INTEGER*8, INTENT(IN) :: statupdn
INTEGER*8, INTENT(OUT):: l

integer*8:: head, tail, middle

head=1; tail=ksize; middle=(head+tail)/2

do while((listqpt(middle).ne.statupdn).and.(head.le.tail))
   if(statupdn.gt.listqpt(middle)) then
      head = middle+1
   else
      tail = middle-1
   endif
   middle = (head+tail)/2
enddo

if(listqpt(middle).eq.statupdn) then
   l=middle;
else
   l=-1;
endif

end

!*****************************************************

integer*8 function factorial(tempn)
implicit none
integer*8:: tempn,tempn1

factorial=1
tempn1=tempn
do while (tempn1.ne.0)
        factorial=factorial*tempn1;
        tempn1=tempn1-1;
enddo
return
end function

!********************************************************

INTEGER*8 function sumeverybit(tempa)
implicit none
integer*8::tempa,tempa1,tempsum

tempsum=0
tempa1=tempa
do while (tempa1.ne.0)
        tempsum=tempsum+mod(tempa1,2)
        tempa1=ISHFT(tempa1,-1)
enddo
sumeverybit=tempsum
return
end function

!********************************************************

integer*8 function sumbeforebit(temps,ks)
implicit none
integer*8::temps,temps1,is,ks

sumbeforebit=0;

if(ks.eq.0) return

temps1=temps
do is=0,ks-1
        if(BTEST(temps1,is).eqv..true.) then
                sumbeforebit=sumbeforebit+1
        endif
enddo
return
end function

!**********************************************************
integer*8 function ksubtract(kk, kkp)
use BettsCluster; use NumOfOrbitalAndElectrons
implicit none
integer*8 :: kk, kkp
integer*8 :: iktemp, kstat, ktempx, ktempy
integer*8 :: px, py

px = kpoint(kkp,1)
py = kpoint(kkp,2)
      ! code to get kkp
      iktemp=0
      kstat=1
      ktempx = mod(kpoint(kk,1)-px+2*kunit,kunit)
      ktempy = mod(kpoint(kk,2)-py+2*kunit,kunit)
      do while (kstat.ne.0)
         if((kpoint(iktemp,1).eq.ktempx).and.(kpoint(iktemp,2).eq.ktempy)) then
            kstat=0
         else
            iktemp=iktemp+1
         endif
         if(iktemp.eq.N) then
            write(*,*) 'iktemp out of bounds, spin up'
            stop
         endif
      enddo
      ksubtract=iktemp
end function


!**********************************************************
integer*8 function ksum(kk, kkp)
use BettsCluster; use NumOfOrbitalAndElectrons
implicit none
integer*8 :: kk, kkp
integer*8 :: iktemp, kstat, ktempx, ktempy
integer*8 :: px, py

px = kpoint(kkp,1)
py = kpoint(kkp,2)
      ! code to get kkp
      iktemp=0
      kstat=1
      ktempx = mod(kpoint(kk,1)+px+kunit,kunit)
      ktempy = mod(kpoint(kk,2)+py+kunit,kunit)
      do while (kstat.ne.0)
         if((kpoint(iktemp,1).eq.ktempx).and.(kpoint(iktemp,2).eq.ktempy)) then
            kstat=0
         else
            iktemp=iktemp+1
         endif
         if(iktemp.eq.N) then
            write(*,*) 'iktemp out of bounds, spin up'
            stop
         endif
      enddo
      ksum=iktemp
end function


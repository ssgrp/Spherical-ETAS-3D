       module mod_mainfuns
       use mod_kinds, only:dp 
       use mpi

       implicit none
       private
       public::func15,probcalc,bandwcalc,outpmat,bkgdcalc,&
            outprob,outrates

       contains
  
!c***********************************************************************
      subroutine func15(n,b,f,h,ifg)
!c-----------------------------------------------------------------------
!c     likelihood function of modified oomori type
!c     etas point process (with trend: not realized).
!c     the space-time <<anisotropic>> version.
!c \lambda(t,x,y)=\mu+\sum k/(t-t_j+c)^p
!c *\exp{ ( (x-x_j)^2+(y-y_j)^2 )/2/\exp{2\alpham_j} }
!c-----------------------------------------------------------------------
!c      implicit real * 8 (a-h,o-z)
      use mod_kinds, only:dp 
      use mpi
      use mod_basicfuns,only:xlamb,xint
      use mod_eq_data

      use mod_mpi_state
      use mod_med_result
      use mod_eq_data
      
      implicit none

      integer,intent(in):: n
      real(dp),intent(in):: b(n)
      real(dp),intent(out)::f
      real(dp),intent(out)::h(n)
      integer,intent(inout)::ifg
      
      
!c      include 'mpif.h'
!      include 'common1.inc'

      real(dp):: gg(20),gg0(20),gtemp(50)

      real(dp)::ff,df1,df2,df3,df4,df5,df6,df7,df8,df9
      real(dp)::ftemp,df1temp,df2temp,df3temp,df4temp,df5temp
      real(dp)::df6temp,df7temp,df8temp,df9temp
      real(dp)::temp,fv0,fv0temp,fvtemp,h1,h2,h3,h4,h5,h6,h7,h8,h9
      
      integer:: i,ierr,j
!c      dimension b(9),h(9),gtemp(50)

!c------------------Parallel part-----------------------------------      
      
      ff=0d0
      df1=0d0
      df2=0d0
      df3=0d0
      df4=0d0 
      df5=0d0
      df6=0d0
      df7=0d0
      df8=0d0
      df9=0d0
 
      ftemp=0d0
      df1temp=0d0
      df2temp=0d0
      df3temp=0d0
      df4temp=0d0
      df5temp=0d0
      df6temp=0d0
      df7temp=0d0
      df8temp=0d0
      df9temp=0d0

      do i=ista,iend
       if(ind(i) == 1)then
        call xlamb(i,b, temp, gtemp)
        if(temp.gt.1d-25)then 
             ftemp=ftemp+log(temp)
         else 
             ftemp=ftemp-100
        endif
        df1temp=df1temp+gtemp(1)/temp 
        df2temp=df2temp+gtemp(2)/temp 
        df3temp=df3temp+gtemp(3)/temp 
        df4temp=df4temp+gtemp(4)/temp 
        df5temp=df5temp+gtemp(5)/temp 
        df6temp=df6temp+gtemp(6)/temp 
        df7temp=df7temp+gtemp(7)/temp 
        df8temp=df8temp+gtemp(8)/temp 
        df9temp=df9temp+gtemp(9)/temp
       endif
      enddo

!       call mpi_reduce(ftemp,ff,1,mpi_double_precision,mpi_sum,0, &
!         mpi_comm_world,ierr)
!       call mpi_reduce(df1temp,df1,1,mpi_double_precision,mpi_sum,0, &
!         mpi_comm_world,ierr)
!       call mpi_reduce(df2temp,df2,1,mpi_double_precision,mpi_sum,0, &
!         mpi_comm_world,ierr)
!       call mpi_reduce(df3temp,df3,1,mpi_double_precision,mpi_sum,0, &
!         mpi_comm_world,ierr)
!       call mpi_reduce(df4temp,df4,1,mpi_double_precision,mpi_sum,0, &
!         mpi_comm_world,ierr)
!       call mpi_reduce(df5temp,df5,1,mpi_double_precision,mpi_sum,0, &
!         mpi_comm_world,ierr)
!       call mpi_reduce(df6temp,df6,1,mpi_double_precision,mpi_sum,0, &
!         mpi_comm_world,ierr)
!       call mpi_reduce(df7temp,df7,1,mpi_double_precision,mpi_sum,0, &
!         mpi_comm_world,ierr)
!       call mpi_reduce(df8temp,df8,1,mpi_double_precision,mpi_sum,0, &
!         mpi_comm_world,ierr)
!       call mpi_reduce(df9temp,df9,1,mpi_double_precision,mpi_sum,0, &
!         mpi_comm_world,ierr)

       call mpi_allreduce(ftemp,ff,1,mpi_double_precision,mpi_sum, &
         mpi_comm_world,ierr)
       call mpi_allreduce(df1temp,df1,1,mpi_double_precision,mpi_sum, &
         mpi_comm_world,ierr)
       call mpi_allreduce(df2temp,df2,1,mpi_double_precision,mpi_sum, &
         mpi_comm_world,ierr)
       call mpi_allreduce(df3temp,df3,1,mpi_double_precision,mpi_sum, &
         mpi_comm_world,ierr)
       call mpi_allreduce(df4temp,df4,1,mpi_double_precision,mpi_sum, &
         mpi_comm_world,ierr)
       call mpi_allreduce(df5temp,df5,1,mpi_double_precision,mpi_sum, &
         mpi_comm_world,ierr)
       call mpi_allreduce(df6temp,df6,1,mpi_double_precision,mpi_sum, &
         mpi_comm_world,ierr)
       call mpi_allreduce(df7temp,df7,1,mpi_double_precision,mpi_sum, &
         mpi_comm_world,ierr)
       call mpi_allreduce(df8temp,df8,1,mpi_double_precision,mpi_sum, &
         mpi_comm_world,ierr)
       call mpi_allreduce(df9temp,df9,1,mpi_double_precision,mpi_sum, &
         mpi_comm_world,ierr)

!       call mpi_bcast(ff,1,mpi_double_precision,0, &
!         mpi_comm_world,ierr)
!       call mpi_bcast(df1,1,mpi_double_precision,0, &
!         mpi_comm_world,ierr)
!       call mpi_bcast(df2,1,mpi_double_precision,0, &
!         mpi_comm_world,ierr)
!       call mpi_bcast(df3,1,mpi_double_precision,0, &
!         mpi_comm_world,ierr)
!       call mpi_bcast(df4,1,mpi_double_precision,0, &
!         mpi_comm_world,ierr)
!       call mpi_bcast(df5,1,mpi_double_precision,0, &
!         mpi_comm_world,ierr)
!       call mpi_bcast(df6,1,mpi_double_precision,0, &
!         mpi_comm_world,ierr)
!       call mpi_bcast(df7,1,mpi_double_precision,0, &
!         mpi_comm_world,ierr)
!       call mpi_bcast(df8,1,mpi_double_precision,0, &
!         mpi_comm_world,ierr)
!       call mpi_bcast(df9,1,mpi_double_precision,0, &
!         mpi_comm_world,ierr)

       do i=1,n
          gg0(i)=0d0
          gg(i)=0d0
       enddo

       fv0temp=0     

!c------- DIVIDING jobs -------------

       do i=ista2,iend2
       call xint(i,b,fvtemp,gg)
           do j=1,n
              gg0(j)=gg0(j)+gg(j)
           enddo
!c           write(*,*)'in task', myrank, i,'gg3=',gg0(3)
           fv0temp=fv0temp+fvtemp 
!C      write(*,*)'in task', myrank, 'fv0=',fvtemp, tz
       enddo

!       call mpi_reduce(fv0temp,fv0,1,mpi_double_precision,mpi_sum,0, &
!         mpi_comm_world,ierr)
!       call mpi_reduce(gg0(1),h1,1,mpi_double_precision,mpi_sum,0, &
!         mpi_comm_world,ierr)
!       call mpi_reduce(gg0(2),h2,1,mpi_double_precision,mpi_sum,0, &
!         mpi_comm_world,ierr)
!       call mpi_reduce(gg0(3),h3,1,mpi_double_precision,mpi_sum,0, &
!         mpi_comm_world,ierr)
!       call mpi_reduce(gg0(4),h4,1,mpi_double_precision,mpi_sum,0, &
!         mpi_comm_world,ierr)
!       call mpi_reduce(gg0(5),h5,1,mpi_double_precision,mpi_sum,0, &
!         mpi_comm_world,ierr)
!       call mpi_reduce(gg0(6),h6,1,mpi_double_precision,mpi_sum,0, &
!         mpi_comm_world,ierr)
!        call mpi_reduce(gg0(7),h7,1,mpi_double_precision,mpi_sum,0, &
!         mpi_comm_world,ierr)
!        call mpi_reduce(gg0(8),h8,1,mpi_double_precision,mpi_sum,0, &
!         mpi_comm_world,ierr)
!        call mpi_reduce(gg0(9),h9,1,mpi_double_precision,mpi_sum,0, &
!         mpi_comm_world,ierr)

!       if(myrank == 0)then 
!           fv0=fv0+xint0*b(1)*b(1)
!            h1=xint0*b(1)*2
!         endif
         
       call mpi_allreduce(fv0temp,fv0,1,mpi_double_precision,mpi_sum, &
         mpi_comm_world,ierr)
       call mpi_allreduce(gg0(1),h1,1,mpi_double_precision,mpi_sum, &
         mpi_comm_world,ierr)
       call mpi_allreduce(gg0(2),h2,1,mpi_double_precision,mpi_sum, &
         mpi_comm_world,ierr)
       call mpi_allreduce(gg0(3),h3,1,mpi_double_precision,mpi_sum, &
         mpi_comm_world,ierr)
       call mpi_allreduce(gg0(4),h4,1,mpi_double_precision,mpi_sum, &
         mpi_comm_world,ierr)
       call mpi_allreduce(gg0(5),h5,1,mpi_double_precision,mpi_sum, &
         mpi_comm_world,ierr)
       call mpi_allreduce(gg0(6),h6,1,mpi_double_precision,mpi_sum, &
         mpi_comm_world,ierr)
        call mpi_allreduce(gg0(7),h7,1,mpi_double_precision,mpi_sum, &
         mpi_comm_world,ierr)
        call mpi_allreduce(gg0(8),h8,1,mpi_double_precision,mpi_sum, &
         mpi_comm_world,ierr)
        call mpi_allreduce(gg0(9),h9,1,mpi_double_precision,mpi_sum, &
         mpi_comm_world,ierr)

!       if(myrank == 0)then 
        fv0=fv0+xint0*b(1)*b(1)
        h1=xint0*b(1)*2
!       endif

       !c      write(*,*)'in task', myrank, 'ff1=',ff1

!       call mpi_bcast(fv0,1,mpi_double_precision,0, &
!         mpi_comm_world,ierr)
!       call mpi_bcast(h1,1,mpi_double_precision,0, &
!         mpi_comm_world,ierr)
!       call mpi_bcast(h2,1,mpi_double_precision,0, &
!         mpi_comm_world,ierr)
!       call mpi_bcast(h3,1,mpi_double_precision,0, &
!         mpi_comm_world,ierr)
!       call mpi_bcast(h4,1,mpi_double_precision,0, &
!         mpi_comm_world,ierr)
!       call mpi_bcast(h5,1,mpi_double_precision,0, &
!         mpi_comm_world,ierr)
!       call mpi_bcast(h6,1,mpi_double_precision,0, &
!         mpi_comm_world,ierr)
!       call mpi_bcast(h7,1,mpi_double_precision,0, &
!         mpi_comm_world,ierr)
!       call mpi_bcast(h8,1,mpi_double_precision,0, &
!         mpi_comm_world,ierr)
!       call mpi_bcast(h9,1,mpi_double_precision,0, &
!         mpi_comm_world,ierr)

        f=-ff+fv0
        h(1)=-df1+h1
        h(2)=-df2+h2
        h(3)=-df3+h3
        h(4)=-df4+h4
        h(5)=-df5+h5
        h(6)=-df6+h6
        h(7)=-df7+h7
        h(8)=-df8+h8
        h(9)=-df9+h9

        if(myrank == 0)write(*,991)f,ff,fv0, (h(i),i=1,9)
!C        if(myrank == 0) then
!C       print *, 'f=', f
!C        print *, 'ff=', ff
!C        print *, 'fv0=', fv0
!C        print *, 'Gradient=', (h(i),i=1,9)
!C        endif
!c        write(*,992)(b(i),i=1,9)
        if(myrank == 0)write(*,992)(b(i),i=1,9)

      return
 991  format(1x,'function value=', 3f12.4, / ' Gradient= ',/9f12.4)
 992  format(1x,'at b=',/9f9.4/)
!c 993  format(1x,9f12.4)
!c 994  format(1x,i8,f19.7)
      end

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C
!C   This subroutine calculates the variable bandwith for each event
!C         given np and minimum bandwidth xlband
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine bandwcalc(np, xlband)

!c     implicit real*8 (a-h, o-z)
      use mpi
      use mod_kinds, only:dp
      use mod_basicfuns, only:iasign
      use mod_geomsphere,only:arcdist
      use mod_eq_data
      use mod_med_result
      use mod_model
      use mod_mpi_state
      implicit none

      integer::np
      real(dp):: xlband
      
!c      include 'mpif.h'
!      include 'common.inc'

!      real(dp):: dbuffer(100), zbdtemp(nnd)
      real(dp),allocatable:: dbuffer(:), zbdtemp(:)

      integer::i, ierr,j, k, kk
      real(dp):: temp

      allocate(dbuffer(np),zbdtemp(nn))
      
!c      real(dp),external::arcdist

        zbdtemp = 0.0d0
        dbuffer = huge(1.0_dp)
        
      do 20 i=1, nn
!CCC    Dividing jobs according to myrank.
       if(iasign(i, nprocs, myrank) == 1)then
        do 5 j=1,nn
             if(i.ne.j)then
                temp=arcdist(xx(i),yy(i),xx(j),yy(j))
                 do 10 k=1,np
                   if(dbuffer(k).gt.temp)then
                       do kk=np,k+1,-1
                             dbuffer(kk)=dbuffer(kk-1)
                       enddo
                       dbuffer(k)=temp
                       goto 15
                   endif 
 10             continue
             endif
 15          continue
 5        enddo       
          zbdtemp(i)=dbuffer(np)
          if(zbdtemp(i).lt.xlband)zbdtemp(i)=xlband
       endif
 20   enddo 
            
      call mpi_allreduce(zbdtemp, zbandw, nn, mpi_double_precision, &
       mpi_sum, mpi_comm_world, ierr)
      
      if(myrank == 0)then
        open(77, file="bandwidth.dat")
        do i = 1, nn
           write(77,997)xx(i),yy(i),zmg(i),zz(i),zdp(i),zbandw(i)
        enddo
        close(77)
      endif
      
      return
 997  format(1x,6f15.6)
      end

!C-----------------------------------------------------------------------
!C
!C  This subroutine calculates the thinning probabilities
!C     zprob(j)=the probability of j-th event being an independent event
!C
!C-----------------------------------------------------------------------

      subroutine probcalc()
!c      implicit real*8 (a-h, o-z)
      use mpi
      use mod_kinds, only:dp
      use mod_basicfuns, only:xlamb,iasign
      use mod_mpi_state
      use mod_eq_data
      use mod_med_result
      use mod_model
      
      use mod_domain
      
      implicit none
      
!c      include 'mpif.h'
!      include 'common.inc'

      real(dp),allocatable::  gtemp(:), zpbtemp(:)
      
!      gtemp(50), zpbtemp(nnd)
      real(dp)::bk,b2
      integer::i, ierr

      allocate(gtemp(n),zpbtemp(nn))
      
      zpbtemp=0.0_dp

       do i = 1, nn
        if(iasign(i, nprocs, myrank) == 1) then 
           bk=zbkgd(i)*x(1)*x(1)
           call xlamb(i,x, b2, gtemp)
           zpbtemp(i)=bk/b2
        endif
       enddo
      
!C      print *, myrank, 'zpbtemp(1)=', zpbtemp(1)

       call mpi_allreduce(zpbtemp, zprob, nn, mpi_double_precision, &
        mpi_sum, mpi_comm_world, ierr)

!C      if(myrank == 0)then
!C           write(*,*)'Probabilities to be independent events:'
!C           print *, zprob(1)
!C      endif
!c     write(*,*)'finishing background calculation'
     
!c 990  format(8(1x,i5,f7.4))
      end

!C-----------------------------------------------------------------------
!C
!C  This subroutine calculates the thinning probabilities
!C     for the probability matrix
!C  output as (i, j, \rho_ij)
!C-----------------------------------------------------------------------

      subroutine outpmat(fp)
!c      implicit real*8 (a-h, o-z)
!c      include 'mpif.h'

        use mod_kinds, only:dp
        use mpi
        use mod_basicfuns,only:xlamb,iasign
        use mod_geomsphere,only:havad
        use mod_specialfun,only:dbeta
        use mod_params,only:pi
        use mod_model
        use mod_eq_data
        use mod_med_result
        use mod_mpi_state

        use mod_smoothing
        use mod_domain

      implicit none
      
      character*80 fp
      
!      include 'param.inc'
!      include 'common.inc'
!c      parameter (pi=3.1415926535897932384626433832795028841971693993751)


      real(dp):: xmu,a2,c,alfa,p,d,q,eta,gamma
      integer::i,j, ierr
      real(dp)::bk,delt, pr1,pr2,pr3,ssig,bbb1,bbb2,havd
      real(dp)::pr4,s,temprob,tmp,tmp1,tmp2
!      real(dp),external::dbeta

          xmu=x(1)**2
          a2=x(2)**2
          c=x(3)**2
          alfa=x(4)**2
          p=x(5)**2
          d=x(6)**2
          q=x(7)**2
          gamma=x(8)**2
          eta=x(9)**2
      
       do i=1,nn
             if(iasign(i, nprocs, myrank)/=1)lambdas(i)=0.0_dp
       enddo

       call mpi_allreduce(MPI_IN_PLACE,lambdas,nn,mpi_double_precision,mpi_sum,&
            mpi_comm_world, ierr)
       
      if(myrank == min(1,nprocs)) then
          open(39, file=fp)
             
          do i=1,nn
            bk=zbkgd(i)*xmu
!            call xlamb(i,x, b2, gtemp)          
            if(i.gt.1)then
            do j=1, i-1
            
               delt=zz(i)-zz(j)

               pr1=exp(alfa*zmg(j))
               pr2=(p-1)/c*(1d0+delt/c)**(-p)
               ssig=d*exp(gamma*zmg(j))
               havd=havad(xx(i),yy(i),xx(j),yy(j))
               bbb1=ssig+havd
               bbb2=((ssig+1d0)**(1-q)-ssig**(1-q))
               pr3=(1d0-q)/4d0/PI/bbb2*bbb1**(-q)          
   
               tmp1 = zdp(j)/depmax
               tmp2 = 1d0-tmp1

               tmp = zdp(i)/depmax
               if(abs(zdp(i))< 1d-12) tmp = 0.5/depmax
               if(zdp(i).ge.depmax-1d-8) tmp=1d0-0.5d0/depmax
               pr4 = dbeta(tmp, eta*tmp1 + 1d0, eta*tmp2 +1d0)/depmax

               s=a2*pr1*pr2*pr3*pr4

               temprob=s/lambdas(i)

               if(temprob.gt.1e-20) write(39,995)i,j,temprob

            enddo
            endif
            write(39,995)i,i,bk/lambdas(i)

        enddo
        close (39)
        endif
 995    format(1x,i8,i8, f20.16)
      end

!c********************************************************************
!c
!c    function for computing the background rate at the position of
!C    all the events, stored in a common shared vector zbkgd       
!c   
!c********************************************************************   

      subroutine bkgdcalc()
!c     implicit real*8 (a-h,o-z)
      use mpi
      use mod_kinds, only:dp 
      use mod_basicfuns,only:pfisher,pgauss,pcauchy,sfr,dfisher,iasign
      use mod_geomsphere,only:havad
      use mod_fintegral,only:srectint,spolyint
      use mod_specialfun,only:dbeta
      use mod_model
      use mod_eq_data
      use mod_med_result
      use mod_mpi_state

      use mod_smoothing
      use mod_domain

      implicit none
      
!c      include 'mpif.h'
!      include 'common.inc'

!      real(dp),external:: dbeta
!      real(dp),external:: srectint,spolyint
      
      real(dp):: w(6)

      real(dp),allocatable::zbkdtmp(:)
      
      real(dp):: sint,havr0,s,s1
      integer::i,ierr,j
      
!c      external pfisher,pgauss,pcauchy
!c      real*8 w(6),zbkdtmp(nnd)
!c      real*8 sint

      allocate(zbkdtmp(nn))
       
      zbkdtmp=0.0d0
      

      do i = 1, nn
         if(iasign(i, nprocs, myrank) == 1)then
           s=0d0
           do j=1,nn
             havr0=havad(xx(j),yy(j),xx(i),yy(i))
             w(1)=zbandw(j)
             s=s+zprob(j)*dfisher(havr0,w)*dbeta(zdp(i)/depmax, &
              bwh*zdp(j)/depmax+1, bwh*(1-zdp(j)/depmax)+1)/depmax
           enddo
           zbkdtmp(i)=s/(tz-tstart)
         endif
      enddo

      call mpi_allreduce(zbkdtmp, zbkgd, nn, mpi_double_precision, &
           mpi_sum, mpi_comm_world, ierr)

        s=0d0
        do i=1, nn
          if(iasign(i, nprocs, myrank) == 1) then
            w(1)=zbandw(i)
            if(ibndtyp == 1)then
              sint=1d0
            end if 
            if(ibndtyp == 2)then
              sint=srectint(pfisher,w,tx,ty,xx(i),yy(i))
              if(indat(i) == 1) sint=sint+1d0
            end if 
            if(ibndtyp == 3)then
               sint=spolyint(pfisher,w,tx,ty,npoly,xx(i),yy(i))
               if(indat(i) == 1) sint=sint+1d0
            end if 
            s=s+zprob(i)*sint

          endif
        enddo 

         s1=0
         call mpi_allreduce(s,s1,1,mpi_double_precision,mpi_sum, &
          mpi_comm_world, ierr)
!c        write(*,*)myrank,s,s1
          
        xint0=s1

      return 
      end subroutine bkgdcalc



      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   The subroutine outprob outputs the probability of each event
!C      being a background event or not
!C
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        subroutine outprob(file1)
!c     implicit real*8 (a-h,o-z)

!        use mod_kinds, only:dp
        use mod_model
        use mod_eq_data
        use mod_med_result
        use mod_mpi_state

        implicit none
        
        character*80 file1

        integer::i
!c        include 'mpif.h'
!        include 'common.inc'

        if(myrank.eq.min(2,nprocs-1))then
             open(35,file=file1)
             do i=1,nn
                write(35,992)i,zprob(i),zbandw(i),zbkgd(i)
             enddo
             close(35)
        endif
       return 

!c 990  format(6(1x,i5,f7.4))
 992  format(1x,i8,3f20.10)
      end subroutine outprob
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C   The subroutine outback outputs the background rate on a lattice of
!C      mx*my grids
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine outrates(file1)
!c      implicit real*8 (a-h,o-z)

      use mpi
      use mod_kinds, only:dp
      use mod_basicfuns,only:dfisher
      use mod_geomsphere,only:havad
      use mod_specialfun,only:dbeta
      use mod_params,only:pi
      
      use mod_model
      use mod_eq_data
      use mod_med_result
      use mod_mpi_state

      use mod_smoothing
      use mod_domain
      
      implicit none
      character*80 file1

!c      include 'mpif.h'
!c      include 'common.inc'
!      include 'param.inc'
      real(dp):: rtx1, rtx2, rty1, rty2, rtz1, rtz2
      common /range/ rtx1, rtx2, rty1, rty2, rtz1, rtz2

      real(dp):: w(10)

!      real(dp),external::dbeta
      real(dp)::xmu, a2, c, alfa, p, d, q, gamma, eta
      real(dp):: bbb1,bbb2,btemp, delt, havd2,havr0
      real(dp):: pr1, pr2, pr3,pr4, s, s1, ssig, tmp,tmp1,tmp2
      real(dp)::tx1, tx2, ty1, ty2, x0,y0, xlb, z0
      
      integer::i, j, kk,ii
      
          xmu=x(1)**2
          a2=x(2)**2
          c=x(3)**2
          alfa=x(4)**2
          p=x(5)**2
          d=x(6)**2
          q=x(7)**2
          gamma=x(8)**2
          eta=x(9)**2

      if(myrank.eq.0) then
         open(34,file=file1)
         write(34,*)mx, my, mz
  
         tx1=rtx1
         tx2=rtx2
         ty1=rty1
         ty2=rty2

         do kk = 1, mz
            z0 = rtz1 + (kk-1)*(rtz2 - rtz1) / (mz-1)
         do i=1,mx
            x0=tx1+(tx2-tx1)/(mx-1)*(i-1)

            tmp = z0/depmax
            if(z0.le.1e-8)tmp= 0.5/depmax
            if(z0.ge.depmax-1e-8)tmp =1d0-0.5d0/depmax
         
            do j=1,my
               y0=ty1+(ty2-ty1)/(my-1)*(j-1)
               s=0d0     
               s1=0d0
               do ii=1,nn
                  havr0=havad(xx(ii),yy(ii),x0,y0)
                  w(1)=zbandw(ii)
                  btemp = dbeta(tmp, bwh*zdp(ii)/depmax+1, &
                      bwh*(1-zdp(ii)/depmax)+1)/depmax 
     
                  s=s+zprob(ii)*dfisher(havr0,w)*btemp
                  s1=s1+dfisher(havr0,w)*btemp
                  
               enddo
               xlb = 0.0 + xmu * s/(tz-tstart)
               do ii=1, nn
                 delt=tz - zz(ii)
                 pr1=exp(alfa*zmg(ii))
                 pr2=(p-1)/c*(1d0+delt/c)**(-p)
                 ssig=d*exp(gamma*zmg(ii))
                 havd2=havad(xx(ii),yy(ii),x0,y0)
                 bbb1=ssig+havd2
                 bbb2=((ssig+1d0)**(1-q)-ssig**(1-q))
                 pr3=(1d0-q)/4d0/PI/bbb2*bbb1**(-q)

                 tmp1 = zdp(ii)/depmax
                 tmp2 = 1d0-tmp1
                 pr4 = dbeta(tmp,eta*tmp1+1d0,eta*tmp2+1d0)/depmax
                 xlb = xlb + a2*pr1*pr2*pr3*pr4

               enddo           
               write(34,991)s/(tz-tstart),s1/(tz-tstart), xlb
            enddo
         enddo
       enddo
       close(34)
      endif
      return 
991   format(1x,3f20.14)
      
      end subroutine outrates
    

    end module mod_mainfuns

!C      file: etas3ds.f90
!C      associated files: setas8fun.f davidn.f sfrint.f 
!C                        spoly.f setas8out.f
      program etas3ds
!c-----------------------------------------------
!c     3D space-time version of the SETAS model
!c     with parallel computing
!c     Spherical version
!c------------------------------------------------

!c      implicit real*8 (a-h, o-z)
      use mpi
!      use mod_kinds, only:dp 
      use mod_mainfuns,only: bkgdcalc,probcalc,bandwcalc,outrates, &
           outprob,outpmat
      use mod_med_result,only:zprob
      use mod_smoothing
      use mod_mpi_state
      use mod_eq_data
      implicit none
      
!c      include 'mpif.h'
!      include 'common.inc'
      character *80 fn
      integer::i, ioptimise,iterative,ierr
      common /optim/ ioptimise


      call mpi_init(ierr)
      call mpi_comm_size(mpi_comm_world, nprocs, ierr)
      call mpi_comm_rank(mpi_comm_world, myrank, ierr)

      call readata()
      call para(nn)
      call bandwcalc(npp, bwm)

      do i = 1, nn
          zprob(i) = 1
      enddo  

      do iterative = 1, 11
            do i = 1, (12-iterative)
            call bkgdcalc() 
            call probcalc()
            enddo
        if(ioptimise == 1) then 
            call dav()
        endif
      enddo

      fn="rates.dat"
      call outrates(fn)
      fn="probs.dat"
      call outprob(fn)
      fn="pmatr.dat"
      call outpmat(fn)

      call mpi_finalize(ierr)

      end


!c********************************************************************
!c
!c
!c********************************************************************
      subroutine para(nk)

!c      implicit real*8 (a-h,o-z)
        use mod_kinds, only:dp
        use mod_mpi_state
      implicit none
      integer::nk
!      include 'common.inc'

      if(myrank == 0) ista = 1
      if(myrank /= 0) ista = int(sqrt(real(nk,kind=dp)*real(nk,kind=dp)* &
         real(myrank,kind=dp)/real(nprocs,kind=dp)))+1 
      if(myrank /= nprocs-1) iend = int(sqrt(real(nk,kind=dp)*   &
        real(nk,kind=dp)*real(myrank+1,kind=dp)/real(nprocs,kind=dp)))

      if(myrank == nprocs-1) iend = nk

      ista2 = nk*myrank/nprocs+1
      iend2 = nk*(myrank+1)/nprocs 
      return
      end

!c********************************************************************
!c Read data and parameters
!c********************************************************************

      subroutine readata() 
!c     implicit real*8 (a-h, o-z)

      use mod_kinds, only: dp 
      use mpi
      use mod_geomsphere, only: antipl,srecse,spolyse,deg2rad
      use mod_params, only: pi

      use mod_model
      use mod_eq_data
      use mod_domain
      use mod_smoothing
      use mod_med_result

      use mod_mpi_state

      implicit none
      
      real(dp):: rtx1, rtx2, rty1, rty2, rtz1, rtz2
      common /optim/ ioptimise
      common /range/ rtx1, rtx2, rty1, rty2, rtz1, rtz2

      
      character *80 hypodata, fmt
      character *160 ftl
      integer::i, ierr,ioptimise, itemp, nkeep
      real(dp)::xmg0,zmin
      real(dp)::xtemp,ytemp,mgtemp,ttemp,deptemp, prev_zz 
      
      if(myrank.eq.0) then

         write(*,*)'Please input the name of the data file:'
         read(*,*)hypodata
         write(*,*)hypodata
         write(*,*) 'Please input data format:'
         read(*,*)fmt
         write(*,*) 'Input the length of the time intervals:'
         read(*,*)tz
         write(*,*) 'Input number of events:'
         read(*,*)nn
         write(*,*) 'Input starting points of time:'
         read(*,*)  zmin
         write(*,*) 'Input the magnitude threshold and tstart:'
         read(*,*) xmg0, tstart

!c---------define the boundary of the region--------------------------
         write(*,*) 'Input boundary type: 1 no-boundary;'
         write(*,*) '                     2 longitude-latitute range;'
         write(*,*) '                     3 sphere polygon;'
         write(*,*) '                     others: error!'
         read(*,*) ibndtyp 
      endif

      call mpi_bcast(ibndtyp,1,mpi_integer,0,   &
         mpi_comm_world,ierr)

           if(ibndtyp==1)then
              if(myrank==0) write(*,*) 'no-boundary:'
              npoly=0
           endif
      
          if(ibndtyp==2)then
             npoly=2          
             call alloc_domain()
          endif

        if(myrank.eq.0.and.ibndtyp==2)then
           write(*,*) 'Input range (lon1, lon2, lat1, lat2):'
           write(*,*) '(longitude: counterclockwise around north pole.)'
           write(*,*) '(latitude: lat1 > lat2 to exclude rectangle)'
           read(*,*)tx(1), tx(2), ty(1), ty(2)
           print *, tx(1), tx(2), ty(1), ty(2)

           
           tx(1)=deg2rad(tx(1))
           ty(1)=deg2rad(ty(1))
           tx(2)=deg2rad(tx(2))
           ty(2)=deg2rad(ty(2)) 
           write(*,*)(tx(i), ty(i), i=1,2)          
        end if

        
        
         if(ibndtyp.eq.3.and.myrank==0)then
           write(*,*) 'Input the number of vertices for the polygon:'
           read(*,*) npoly
        endif

        call mpi_bcast(npoly,1,mpi_integer,0,        &
         mpi_comm_world,ierr)
        call alloc_domain()
        
        if(ibndtyp.eq.3) then
           if( myrank==0)then         
             write(*,*) 'Input the coordination of the polygon vertices:'
             write(*,*) 'The north pole should not on the boundary.If so,'
             write(*,*) 'add a small polygon to exclude or include it.'
             write(*,*) 'The interior is on the left.'
             do i=1, npoly
                read(*,*)tx(i), ty(i)
                tx(i)=deg2rad(tx(i))
                ty(i)=deg2rad(ty(i))
             enddo
             tx(npoly+1)=tx(1)
             ty(npoly+1)=ty(1)
          end if
          
       endif
       
       call mpi_bcast(tx,npoly+1,mpi_double_precision,0, &
          mpi_comm_world,ierr)
       call mpi_bcast(ty,npoly+1,mpi_double_precision,0, &
         mpi_comm_world,ierr)

         if (ibndtyp.gt.3 .or. ibndtyp.lt.1) then
         stop
         endif

         if(myrank==0) then
         write(*,*)'Input maximum depth:'
         read(*,*) depmax

         write(*,*)'Input numbers of grids for background:'
         read(*,*) mx,my, mz
         write(*,*)'Input ranges of grids:'
         read(*,*) rtx1, rtx2, rty1, rty2, rtz1, rtz2
         print *, rtx1, rtx2, rty1, rty2, rtz1, rtz2
           rtx1=deg2rad(rtx1)
           rtx2=deg2rad(rtx2)
           rty1=deg2rad(rty1)
           rty2=deg2rad(rty2)
        
         open(11, file=hypodata)
         read(11,*)ftl
         write(*,997)ftl

         
         nkeep=0        
         prev_zz = -huge(1.0_dp)

15       read(11,*,end=25)itemp,xtemp,ytemp,mgtemp,ttemp,deptemp
         if(ttemp.gt.tz+zmin)goto 25
          if(ttemp.lt.prev_zz)then
                  write(*,*)'Reverse data:', prev_zz,ttemp
                  stop
          endif
          prev_zz=ttemp
          
          if(mgtemp.gt.xmg0-0.000001.and.ttemp.gt.zmin.and.deptemp.lt.&
               depmax) then
               nkeep=nkeep+1
          endif
          goto 15
25       continue

          nn=nkeep
          write(*,*) nkeep
       endif

       call mpi_bcast(nn,1,mpi_integer,0,           &
            mpi_comm_world,ierr)
       
       call alloc_eq_data()
       call alloc_med_result()


       
       if(myrank==0)then
          rewind(11)
          
         read(11,*)ftl
        
         i=1
         prev_zz = -huge(1.0_dp)

115       read(11,*,end=125)itemp,xtemp,ytemp,mgtemp,ttemp,deptemp
!c          itemp,xx(i),yy(i),zmg(i),zz(i),zdp(i)         
           if(ttemp.gt.tz+zmin)goto 125
           if(ttemp.lt.prev_zz)then
                  write(*,*)'Reverse data:', prev_zz, ttemp
                  stop
          endif
          prev_zz=ttemp

          if(mgtemp.gt.xmg0-0.000001.and.ttemp.gt.zmin.and.deptemp.lt.&
                 depmax) then
               
               xx(i)=deg2rad(xtemp)
               yy(i)=deg2rad(ytemp)
               zz(i)=ttemp
               zmg(i) = mgtemp
               zdp(i) = deptemp
               
               i=i+1 
            endif
            
            goto 115
 125      continue
            nn=i-1
          
!C 
         write(46,*) nn
         nnc = 0

!C     calculate antipoles for each event for Types 2 and 3 boundaries.

        if(ibndtyp.gt.1)then
           do i=1,nn
              ind(i)=1
              indat(i)=1
              call antipl(xx(i),yy(i),xa(i), ya(i))
           enddo
        end if
        
        if(ibndtyp.eq.2)then
           call srecse(tx(1),tx(2),ty(1),ty(2),xx,yy,nn,ind)
           call srecse(tx(1),tx(2),ty(1),ty(2),xa,ya,nn,indat)
        end if
      
         if(ibndtyp.eq.3)then
            call spolyse(tx,ty,npoly,xx,yy,nn,ind)
            call spolyse(tx,ty,npoly,xa,ya,nn,indat)
         end if
         
         do i=1, nn
            if(ind(i).ge.0) ind(i)=1
            if(ind(i).lt.0) ind(i)=0
            if(indat(i).ge.0) indat(i)=1
            if(indat(i).lt.0) indat(i)=0
            write(48, *)i,xx(i),yy(i),zmg(i),zz(i),ind(i),indat(i),&
                zdp(i)
         enddo

         do 10 i=1,nn
            if(zz(i).gt.tz+zmin) ind(i)=0 
            if(zz(i).lt. tstart) ind(i)=0
            zmg(i)=zmg(i)-xmg0
            zz(i)=zz(i)-zmin
            zprob(i)=1.0d0
            write(46,*)i,xx(i),yy(i),zmg(i),zz(i),ind(i),indat(i),zdp(i)

           if(ind(i).eq.1)then
             nnc=nnc+1
             write(47,*)i,xx(i),yy(i),zmg(i),zz(i),ind(i),indat(i), &
               zdp(i)
           endif
   10   continue
        if(myrank.eq.0)print *, nnc, ' of ', nn, ' events selected.'
        tstart=tstart-zmin
        close(46)
        close(47)
        close(48)

      endif

!c---------------------------------------------------------------
!c
!c     Broadcast all the data to other process from process 0
!c
!c---------------------------------------------------------------

      call mpi_bcast(nnc,1,mpi_integer,0,          &
         mpi_comm_world,ierr)
      call mpi_bcast(xx,nn,mpi_double_precision,0, &
         mpi_comm_world,ierr)
      call mpi_bcast(yy,nn,mpi_double_precision,0, &
         mpi_comm_world,ierr)
      call mpi_bcast(zz,nn,mpi_double_precision,0, &
         mpi_comm_world,ierr)
      call mpi_bcast(zmg,nn,mpi_double_precision,0,&
         mpi_comm_world,ierr)
      call mpi_bcast(zdp,nn,mpi_double_precision,0,&
         mpi_comm_world,ierr)
      call mpi_bcast(tz,1,mpi_double_precision,0, &
         mpi_comm_world,ierr)
      call mpi_bcast(tstart,1,mpi_double_precision,0, &
         mpi_comm_world,ierr)
      call mpi_bcast(depmax,1,mpi_double_precision,0, &
         mpi_comm_world,ierr)
      call mpi_bcast(ind,nn,mpi_integer,0,    &
         mpi_comm_world,ierr)
      call mpi_bcast(indat,nn,mpi_integer,0,  &
         mpi_comm_world,ierr)
      call mpi_bcast(zprob,nn,mpi_double_precision,0, &
         mpi_comm_world,ierr)

      call mpi_bcast(mx,1,mpi_integer,0,  &
         mpi_comm_world,ierr)
      call mpi_bcast(my,1,mpi_integer,0,  &
         mpi_comm_world,ierr)
      call mpi_bcast(mz,1,mpi_integer,0,  &
         mpi_comm_world,ierr)
      call mpi_bcast(rtx1,1,mpi_double_precision,0, &
         mpi_comm_world,ierr)
      call mpi_bcast(rtx2,1,mpi_double_precision,0, & 
         mpi_comm_world,ierr)
      call mpi_bcast(rty1,1,mpi_double_precision,0, &
         mpi_comm_world,ierr)
      call mpi_bcast(rty2,1,mpi_double_precision,0, &
         mpi_comm_world,ierr)
      call mpi_bcast(rtz1,1,mpi_double_precision,0, &
         mpi_comm_world,ierr)
      call mpi_bcast(rtz2,1,mpi_double_precision,0, &
         mpi_comm_world,ierr)

!C-------------------------------------------------------
!c    Read and send model parameters
!C------------------------------------------------------

       n=9
       if(myrank == 0)then 
            write(*,*)'Please input initial parameters:'
            write(*,*)'(mu, A, c, alfa, p, d, q, gamma, eta)'  
             read(*,*) (b(i),i=1,n)
            write(6,1020) n
            write(6,1030)  (b(i),i=1,n)
            write(*,*)'Please input n_p'
            write(*,*)' and minimum bandwidth(in degree)'
            write(*,*)' and shape parameter for depth smoothing'
            read(*,*)npp,bwm,bwh
            bwm = bwm /180d0*pi
            write(*,*)'Optimize?(1 for Yes, 0 for No)'
            read(*,*)ioptimise
            print *, npp, bwm, bwh, ioptimise
       endif

       call mpi_bcast(b,n,mpi_double_precision,0,  &
             mpi_comm_world,ierr)
       call mpi_bcast(npp,1,mpi_integer,0,   &
             mpi_comm_world,ierr)
       call mpi_bcast(bwm,1,mpi_double_precision,0, &
             mpi_comm_world,ierr)
       call mpi_bcast(bwh,1,mpi_double_precision,0, &
             mpi_comm_world,ierr)
       call mpi_bcast(ioptimise,1,mpi_integer,0,  &
             mpi_comm_world,ierr)

      do i=1,n
         x(i)=sqrt(b(i))
      enddo

!c 991  format(1x,i6,4f12.4,4f5.2)
!c 993  format(1x,'mx, my= ',i7,i7)
!c 995  format('tx,ty,tz,xmin,ymin,xmg0,zmin,tsta',/8f10.3,
!c     &/' nn =',i6,' nnc=',i6)
 997  format(1x, 'Data set  ',10a8)
 1020 format(1x,3x,'input data'/1h ,5x,'n=',i3,3x,'itr=',i5,3x, &
              'ier=',i3,3x,'eps=',1pe16.7)
 1030 format(1x,'x=',5e14.5)
      end

!c********************************************************************
!c
!c********************************************************************

       subroutine dav()

!c       implicit real * 8 (a-h,o-z)
!       use mod_kinds, only:dp
       use mod_opt,only:davidn
       use mod_mainfuns,only:func15

       use mod_mpi_state
       use mod_model
       
       implicit none
       
!c       include 'common.inc'
       integer::i
 
       do i=1,n
         x(i)=sqrt(b(i))
       enddo
       call davidn(x,n,func15)
     
       do 80 i=1,n
         b(i)=x(i)**2
 80    continue

      if(myrank.eq.0)then 
         write(6,1040) (b(i),i=1,n)
         open(22,file='para')
         write(22,1110)(b(i),i=1,n)
      endif
      return

!c 991  format(1x,f16.8)
!c 999  format(1x, 8f12.8)
!c 1020 format(1x,3x,'input data'/1h ,5x,'n=',i3,3x,'itr=',i5,3x,
!c     2          'ier=',i3,3x,'eps=',1pe16.7)
!c 1030 format(1x,'x=',5e14.5)
 1040 format(1x, /' mle = ',5e12.5/('       ',9e12.5))
 1110 format(1h , 8f20.9)
      end

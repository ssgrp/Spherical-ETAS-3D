C      file: etas3ds.f
C      associated files: setas8fun.f davidn.f sfrint.f 
C                        spoly.f setas8out.f

      programme eats3ds
c-----------------------------------------------
c     3D space-time version of the SETAS model
c     with parallel computing
c     Spherical version
c------------------------------------------------

      implicit real*8 (a-h, o-z)
      include 'mpif.h'
      include 'common.inc'
      character *80 fn
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


c********************************************************************
c
c
c********************************************************************
      subroutine para(nk)

      implicit real*8 (a-h,o-z)
      include 'common.inc'

      if(myrank == 0) ista = 1
      if(myrank /= 0) ista = sqrt(float(nk)*float(nk)*
     &          float(myrank)/float(nprocs))+1 
      if(myrank /= nprocs-1) iend = sqrt(float(nk)*float(nk)*
     &    float(myrank+1)/float(nprocs))

      if(myrank == nprocs-1) iend = nk

      ista2 = nk*myrank/nprocs+1
      iend2 = nk*(myrank+1)/nprocs 
      return
      end

c********************************************************************
c Read data and parameters
c********************************************************************

      subroutine readata() 
      implicit real*8 (a-h, o-z)
      include 'mpif.h'
      include 'common.inc'
      include 'param.inc'
      real*8 rtx1, rtx2, rty1, rty2, rtz1, rtz2
      common /optim/ ioptimise
      common /range/ rtx1, rtx2, rty1, rty2, rtz1, rtz2
      character *80 hypodata, fmt
      character *160 ftl

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

c---------define the boundary of the region--------------------------
         write(*,*) 'Input boundary type: 1 no-boundary;'
         write(*,*) '                     2 longitude-latitute range;'
         write(*,*) '                     3 sphere polygon;'
         write(*,*) '                     others: error!'
         read(*,*) ibndtyp 
         if(ibndtyp.eq.1)then
           write(*,*) 'no-boundary:'
           npoly=0
         endif 
         if(ibndtyp.eq.2)then
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
           npoly=2
         end if 
         if(ibndtyp.eq.3)then
           write(*,*) 'Input the number of vertices for the polygon:'
           read(*,*) npoly
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
         if (ibndtyp.gt.3 .or. ibndtyp.lt.1) then
         stop
         endif

         write(*,*)'Input maximum depth:'
         read(*,*) depmax

         write(*,*)'Input numbers of grids for background:'
         read(*,*)mx,my, mz
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

         i=1
 15      read(11,*,end=25)itemp,xx(i),yy(i),zmg(i),zz(i),zdp(i)
           xx(i)=deg2rad(xx(i))
           yy(i)=deg2rad(yy(i))
           if(zz(i).gt.tz+zmin)goto 25
         if(i.gt.1) then
              if(zz(i).lt.zz(i-1))then
                  write(*,*)'Reverse data:', zz(i-1), zz(i)
                  stop
              endif
         endif

            if(zmg(i).gt.xmg0-0.000001.and.zz(i).gt.zmin.and.zdp(i).lt.
     &           depmax) then
               i=i+1 
            endif
            goto 15
 25      continue
         nn=i-1
C 
         write(46,*) nn
         nnc = 0

C     calculate antipoles for each event for Types 2 and 3 boundaries.

c        if(ibndtyp.gt.1)then
           do i=1,nn
              ind(i)=1
              indat(i)=1
              call antipl(xx(i),yy(i),xa(i), ya(i))
           enddo
c        end if
        
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
            write(48, *)i,xx(i),yy(i),zmg(i),zz(i),ind(i),indat(i),
     &           zdp(i)
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
              write(47,*)i,xx(i),yy(i),zmg(i),zz(i),ind(i),indat(i),
     &          zdp(i)
           endif
   10   continue
        if(myrank.eq.0)print *, nnc, ' of ', nn, ' events selected.'
        tstart=tstart-zmin
        close(46)
        close(47)
        close(48)

      endif

c---------------------------------------------------------------
c
c     Broadcast all the data to other process from process 0
c
c---------------------------------------------------------------

      call mpi_bcast(nn,1,mpi_integer,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(nnc,1,mpi_integer,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(xx,nn,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(yy,nn,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(zz,nn,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(zmg,nn,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(zdp,nn,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(npoly,1,mpi_integer,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(tx,npoly+1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(ty,npoly+1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(tz,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(tstart,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(depmax,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(ind,nn,mpi_integer,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(indat,nn,mpi_integer,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(zprob,nn,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(ibndtyp,1,mpi_integer,0,
     &    mpi_comm_world,ierr)

      call mpi_bcast(mx,1,mpi_integer,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(my,1,mpi_integer,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(mz,1,mpi_integer,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(rtx1,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(rtx2,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(rty1,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(rty2,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(rtz1,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
      call mpi_bcast(rtz2,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)

C-------------------------------------------------------
c    Read and send model parameters
C------------------------------------------------------

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

       call mpi_bcast(b,n,mpi_double_precision,0,
     &        mpi_comm_world,ierr)
       call mpi_bcast(npp,1,mpi_integer,0,
     &        mpi_comm_world,ierr)
       call mpi_bcast(bwm,1,mpi_double_precision,0,
     &        mpi_comm_world,ierr)
       call mpi_bcast(bwh,1,mpi_double_precision,0,
     &        mpi_comm_world,ierr)
       call mpi_bcast(ioptimise,1,mpi_integer,0,
     &        mpi_comm_world,ierr)

      do i=1,n
         x(i)=sqrt(b(i))
      enddo

 991  format(1x,i6,4f12.4,4f5.2)
 993  format(1x,'mx, my= ',i7,i7)
 995  format('tx,ty,tz,xmin,ymin,xmg0,zmin,tsta',/8f10.3,
     &/' nn =',i6,' nnc=',i6)
 997  format(1x, 'Data set  ',10a8)
 1020 format(1x,3x,'input data'/1h ,5x,'n=',i3,3x,'itr=',i5,3x,
     2          'ier=',i3,3x,'eps=',1pe16.7)
 1030 format(1x,'x=',5e14.5)
      end

c********************************************************************
c
c********************************************************************

       subroutine dav()

       implicit real * 8 (a-h,o-z)
       include 'common.inc'
       external func15
 
       do i=1,n
         x(i)=sqrt(b(i))
       enddo
       call davidn(x,n,0,func15)
     
       do 80 i=1,n
         b(i)=x(i)**2
 80    continue

      if(myrank.eq.0)then 
         write(6,1040) (b(i),i=1,n)
         open(22,file='para')
         write(22,1110)(b(i),i=1,n)
      endif
      return

 991  format(1x,f16.8)
 999  format(1x, 8f12.8)
 1020 format(1x,3x,'input data'/1h ,5x,'n=',i3,3x,'itr=',i5,3x,
     2          'ier=',i3,3x,'eps=',1pe16.7)
 1030 format(1x,'x=',5e14.5)
 1040 format(1x, /' mle = ',5e12.5/('       ',9e12.5))
 1110 format(1h , 8f20.9)
      end


c***********************************************************************
      subroutine func15(n,b,f,h,ifg)
c-----------------------------------------------------------------------
c     likelihood function of modified oomori type
c     etas point process (with trend: not realized).
c     the space-time <<anisotropic>> version.
c \lambda(t,x,y)=\mu+\sum k/(t-t_j+c)^p
c *\exp{ ( (x-x_j)^2+(y-y_j)^2 )/2/\exp{2\alpham_j} }
c-----------------------------------------------------------------------
      implicit real * 8 (a-h,o-z)

      include 'mpif.h'
      include 'common1.inc'

      dimension b(9),h(9),gtemp(50)
      real*8 gg(20),gg0(20)

c------------------Parallel part-----------------------------------      
      
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
       if(ind(i).eq.1)then
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

      call mpi_reduce(ftemp,ff,1,mpi_double_precision,mpi_sum,0,
     &    mpi_comm_world,ierr)
       call mpi_reduce(df1temp,df1,1,mpi_double_precision,mpi_sum,0,
     &    mpi_comm_world,ierr)
       call mpi_reduce(df2temp,df2,1,mpi_double_precision,mpi_sum,0,
     &    mpi_comm_world,ierr)
       call mpi_reduce(df3temp,df3,1,mpi_double_precision,mpi_sum,0,
     &    mpi_comm_world,ierr)
       call mpi_reduce(df4temp,df4,1,mpi_double_precision,mpi_sum,0,
     &    mpi_comm_world,ierr)
       call mpi_reduce(df5temp,df5,1,mpi_double_precision,mpi_sum,0,
     &    mpi_comm_world,ierr)
       call mpi_reduce(df6temp,df6,1,mpi_double_precision,mpi_sum,0,
     &    mpi_comm_world,ierr)
       call mpi_reduce(df7temp,df7,1,mpi_double_precision,mpi_sum,0,
     &    mpi_comm_world,ierr)
       call mpi_reduce(df8temp,df8,1,mpi_double_precision,mpi_sum,0,
     &    mpi_comm_world,ierr)
       call mpi_reduce(df9temp,df9,1,mpi_double_precision,mpi_sum,0,
     &    mpi_comm_world,ierr)

       call mpi_bcast(ff,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
       call mpi_bcast(df1,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
       call mpi_bcast(df2,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
       call mpi_bcast(df3,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
       call mpi_bcast(df4,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
       call mpi_bcast(df5,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
       call mpi_bcast(df6,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
       call mpi_bcast(df7,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
       call mpi_bcast(df8,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
       call mpi_bcast(df9,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)

       do i=1,n
          gg0(i)=0d0
          gg(i)=0d0
       enddo

       fv0temp=0     

c------- DIVIDING jobs -------------

       do i=ista2,iend2
       call xint(i,b,fvtemp,gg)
           do j=1,n
              gg0(j)=gg0(j)+gg(j)
           enddo
c           write(*,*)'in task', myrank, i,'gg3=',gg0(3)
           fv0temp=fv0temp+fvtemp 
C      write(*,*)'in task', myrank, 'fv0=',fvtemp, tz
       enddo

       call mpi_reduce(fv0temp,fv0,1,mpi_double_precision,mpi_sum,0,
     &    mpi_comm_world,ierr)
       call mpi_reduce(gg0(1),h1,1,mpi_double_precision,mpi_sum,0,
     &    mpi_comm_world,ierr)
       call mpi_reduce(gg0(2),h2,1,mpi_double_precision,mpi_sum,0,
     &    mpi_comm_world,ierr)
       call mpi_reduce(gg0(3),h3,1,mpi_double_precision,mpi_sum,0,
     &    mpi_comm_world,ierr)
       call mpi_reduce(gg0(4),h4,1,mpi_double_precision,mpi_sum,0,
     &    mpi_comm_world,ierr)
       call mpi_reduce(gg0(5),h5,1,mpi_double_precision,mpi_sum,0,
     &    mpi_comm_world,ierr)
       call mpi_reduce(gg0(6),h6,1,mpi_double_precision,mpi_sum,0,
     &    mpi_comm_world,ierr)
        call mpi_reduce(gg0(7),h7,1,mpi_double_precision,mpi_sum,0,
     &    mpi_comm_world,ierr)
        call mpi_reduce(gg0(8),h8,1,mpi_double_precision,mpi_sum,0,
     &    mpi_comm_world,ierr)
        call mpi_reduce(gg0(9),h9,1,mpi_double_precision,mpi_sum,0,
     &    mpi_comm_world,ierr)

       if(myrank.eq.0)then 
            fv0=fv0+xint0*b(1)*b(1)
            h1=xint0*b(1)*2
       endif
c      write(*,*)'in task', myrank, 'ff1=',ff1

       call mpi_bcast(fv0,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
       call mpi_bcast(h1,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
       call mpi_bcast(h2,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
       call mpi_bcast(h3,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
       call mpi_bcast(h4,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
       call mpi_bcast(h5,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
       call mpi_bcast(h6,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
       call mpi_bcast(h7,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
       call mpi_bcast(h8,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)
       call mpi_bcast(h9,1,mpi_double_precision,0,
     &    mpi_comm_world,ierr)

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

C        if(myrank.eq.0)write(*,991)f,ff,fv0, (h(i),i=1,9)
C        if(myrank.eq.0) then
C       print *, 'f=', f
C        print *, 'ff=', ff
C        print *, 'fv0=', fv0
C        print *, 'Gradient=', (h(i),i=1,9)
C        endif
C        write(*,992)(b(i),i=1,9)
C        if(myrank.eq.0)write(*,992)(b(i),i=1,9)

      return
 991  format(1x,'function value=', 3f12.4, / ' Gradient= ',/9f12.4)
 992  format(1x,'at b=',/9f9.4/)
 993  format(1x,9f12.4)
 994  format(1x,i8,f19.7)
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C   This subroutine calculates the variable bandwith for each event
C         given np and minimum bandwidth xlband
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine bandwcalc(np, xlband)

      implicit real*8 (a-h, o-z)
      include 'mpif.h'
      include 'common.inc'

      real*8 dbuffer(100), zbdtemp(nnd)

      do i=1, nn
        zbdtemp(i) =0.0d0
      enddo 

      do 20 i=1, nn
CCC    Dividing jobs according to myrank.
       if(iasign(i, nprocs, myrank).eq.1)then
        do kkk=1,np
           dbuffer(kkk)=10d39 
        enddo
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
            
       call mpi_allreduce(zbdtemp, zbandw, nn, mpi_double_precision, 
     & mpi_sum, mpi_comm_world, ierr)
      
      if(myrank.eq.0)then
        open(77, file="bandwidth.dat")
        do i = 1, nn
           write(77,997)xx(i),yy(i),zmg(i),zz(i),zdp(i),zbandw(i)
        enddo
        close(77)
      endif
      
      return
 997  format(1x,6f15.6)
      end

C-----------------------------------------------------------------------
C
C  This subroutine calculates the thinning probabilities
C     zprob(j)=the probability of j-th event being an independent event
C
C-----------------------------------------------------------------------

      subroutine probcalc()
      implicit real*8 (a-h, o-z)
      include 'mpif.h'
      include 'common.inc'

      real*8 gtemp(50), zpbtemp(nnd)

      do i=1, nn
          zpbtemp(i)=0.0d0
      enddo

       do i = 1, nn
        if(iasign(i, nprocs, myrank).eq.1) then 
           bk=zbkgd(i)*x(1)*x(1)
           call xlamb(i,x, b2, gtemp)
           zpbtemp(i)=bk/b2
        endif
       enddo
      
C      print *, myrank, 'zpbtemp(1)=', zpbtemp(1)

      call mpi_allreduce(zpbtemp, zprob, nn, mpi_double_precision, 
     &  mpi_sum, mpi_comm_world, ierr)

C      if(myrank.eq.0)then
C           write(*,*)'Probabilities to be independent events:'
C           print *, zprob(1)
C      endif
c     write(*,*)'finishing background calculation'
     
 990  format(8(1x,i5,f7.4))
      end

C-----------------------------------------------------------------------
C
C  This subroutine calculates the thinning probabilities
C     for the probability matrix
C  output as (i, j, \rho_ij)
C-----------------------------------------------------------------------

      subroutine outpmat(fp)
      implicit real*8 (a-h, o-z)
c      include 'mpif.h'
      include 'common.inc'
      parameter (pi=3.1415926535897932384626433832795028841971693993751)

      character*80 fp

      real*8 gtemp(10)

      if(myrank.eq.1) then
          open(39, file=fp)
      
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
            bk=zbkgd(i)*xmu
            call xlamb(i,x, b2, gtemp)
           
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
               if(zdp(i).eq.0) tmp = 0.5/depmax
               if(zdp(i).ge.depmax-1e-8) tmp=1d0-0.5d0/depmax
               pr4 = dbeta(tmp, eta*tmp1 + 1d0, eta*tmp2 +1d0)/depmax

               s=a2*pr1*pr2*pr3*pr4

               temprob=s/b2

               if(temprob.gt.1e-20) write(39,995)i,j,temprob

            enddo
            endif
            write(39,995)i,i,bk/b2

        enddo
        close (39)
        endif
 995    format(1x,i8,i8, f20.16)
      end

c-------------------------------------------------------
      function iasign(im, np, mr)
      implicit real*8 (a-h, o-z)
      
      ia = 0 
      k = mod(im-1, 2*np) 
      if(k.eq.mr.or. k.eq. 2*np-mr-1) ia = 1
      iasign = ia
      
      return 
      end

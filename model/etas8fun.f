c********************************************************************
c
c    function for computing the background rate at the position of
C    all the events, stored in a common shared vector zbkgd       
c   
c********************************************************************   

      subroutine bkgdcalc()
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
      include 'common.inc'

      external pfisher,pgauss,pcauchy
      real*8 w(6),zbkdtmp(nnd)
      real*8 sint

      do i=1, nn
         zbkdtmp(i)=0.0d0
      enddo

      do i = 1, nn
         if(iasign(i, nprocs, myrank).eq.1)then
           s=0d0
           do j=1,nn
             havr0=havad(xx(j),yy(j),xx(i),yy(i))
             w(1)=zbandw(j)
             s=s+zprob(j)*dfisher(havr0,w)*dbeta(zdp(i)/depmax,
     &       bwh*zdp(j)/depmax+1, bwh*(1-zdp(j)/depmax)+1)/depmax
           enddo
           zbkdtmp(i)=s/(tz-tstart)
         endif
      enddo

      call mpi_allreduce(zbkdtmp, zbkgd, nn, mpi_double_precision,
     &        mpi_sum, mpi_comm_world, ierr)

        s=0d0
        do i=1, nn
          if(iasign(i, nprocs, myrank).eq.1) then
            w(1)=zbandw(i)
            if(ibndtyp.eq.1)then
              sint=1d0
            end if 
            if(ibndtyp.eq.2)then
              sint=srectint(pfisher,w,tx,ty,xx(i),yy(i))
              if(indat(i).eq.1) sint=sint+1d0
            end if 
            if(ibndtyp.eq.3)then
               sint=spolyint(pfisher,w,tx,ty,npoly,xx(i),yy(i))
               if(indat(i).eq.1) sint=sint+1d0
            end if 
            s=s+zprob(i)*sint

          endif
        enddo 

C         call mpi_allreduce(s,xint0,1,mpi_double_precision,mpi_sum,
C     &     mpi_comm_world, ierr)
         s1=0
c        write(*,*)myrank,s,s1
         call mpi_reduce(s,s1,1,mpi_double_precision,mpi_sum,0,
     &     mpi_comm_world, ierr)
         call mpi_bcast(s1,1,mpi_double_precision,0,
     &     mpi_comm_world, ierr)
c        write(*,*)myrank,s1
          
        xint0=s1

      return 
      end


C********************************************************************
C   pfisher calculates the integral of the von Mises-Fisher 
C   kernel function
C     with bandwidth w(1) from 0 to r
c   Input:      
C     havr: Haversine of r
c   Output:
c     pfisher: Cumulative probability at r    
C********************************************************************   
       real*8 function pfisher(havr,w)

       implicit real*8 (a-h,o-z)
       real*8 w(6)
       parameter(PI=3.14159265358979323846264338327d0)

       p1fisher_p1=1d0-exp(-havr/2d0/w(1)**2)
       p2fisher_p2=1d0-exp(-0.5d0/w(1)**2)

       pfisher=0.5d0/PI*p1fisher_p1/p2fisher_p2

       return
       end

C********************************************************************



C*******************************************************************
C     sfr(r)=\int f(r)r dr
C    Input: havr haversine of r      
c*******************************************************************
              
        real*8 function sfr(havr,w)
        implicit real*8 (a-h, o-z)
        real*8 w(6),havr
        parameter(PI=3.14159265358979323846264338327d0)

        gamm=w(1)
        d=w(2)
        q=w(3) 
        xmag=w(4)

        ssig=d*exp(gamm*xmag)
        sfr_1=havr+ssig 
        sfr_2=(sfr_1)**(1d0-q)-(ssig)**(1d0-q)
        sfr_3=(ssig+1d0)**(1d0-q)-(ssig)**(1d0-q)

        sfr_23=sfr_2/sfr_3
        sfr=sfr_2/sfr_3/pi/2d0

        return
        end

C********************************************************************
C   dfisher calculates the value of the von Mises-Fisher kernel function
C     at r with bandwidth w(1)
c   Input:      
C     havr: Haversine of r
c   Output:
c     dfisher: density at r    
C
C********************************************************************   
       real*8 function dfisher(havr,w)
       implicit real*8 (a-h,o-z)
       real*8 w(6)
       parameter (PI=3.14159265358979323846264338327)

       d1fisher_d1=exp(-(havr/2d0/w(1)**2))
       d2fisher_d2=1d0-exp(-0.5d0/w(1)**2)

       dfisher=1/8d0/PI/w(1)**2*d1fisher_d1/d2fisher_d2

       return
       end

C*******************************************************************
C       dgamma_sfr(r)=d\int f(r)r dr \over  d\gamma
c*******************************************************************

        real*8 function  dgamma_sfr(havr,param)
        implicit real*8 (a-h, o-z)
        real*8 param(6)
        parameter(PI=3.14159265358979323846264338327d0)

        gamm=param(1)
        d=param(2)
        q=param(3)
        xmag=param(4)

        ssig=d*exp(gamm*xmag)
        sfr_g1=havr+ssig
        if(abs(sfr_g1**(1d0-q)-ssig**(1d0-q)).lt.1d-10) then
           sfr_g2=q/(q-1d0)/ssig     
         else
           sfr_g2=(sfr_g1**(-q)-ssig**(-q))/(sfr_g1**(1d0-q)
     &             -ssig**(1d0-q))
        endif
        sfr_g3=((ssig+1d0)**(-q)-(ssig)**(-q))/((ssig+1d0)**(1d0-q)
     &          -(ssig)**(1d0-q))

        dgamma_sfr=sfr(havr,param)*(1d0-q)*(sfr_g2-sfr_g3)*xmag*ssig

        return
        end

C*******************************************************************
C       dd_sfr(r)=d\int f(r)r dr \over  d\d
c*******************************************************************
              
        real*8 function dd_sfr(havr,param)
        implicit real*8 (a-h, o-z)
        real*8 param(6)
        parameter(PI=3.14159265358979323846264338327d0)

        gamm=param(1)
        d=param(2)
        q=param(3)
        xmag=param(4)

        ssig=d*exp(gamm*xmag)        
        sfr_d1=havr+ssig
        if(abs(sfr_d1**(1d0-q)-ssig**(1d0-q)).lt.1d-10) then
           sfr_d2 = q/(q-1d0)/ssig
         else
          sfr_d2=(sfr_d1**(-q)-ssig**(-q))/(sfr_d1**(1d0-q)
     &          -ssig**(1d0-q))
        endif
       
        sfr_d3=((ssig+1d0)**(-q)-(ssig)**(-q))/((ssig+1d0)**(1-q)
     &         -(ssig)**(1-q))

        dd_sfr=sfr(havr,param)*(1d0-q)*(sfr_d2-sfr_d3)*ssig/d

        return
        end

C*******************************************************************
C       dq_sfr(r)=d\int sf(r)r dr \over  d\q
c*******************************************************************
              
        real*8 function dq_sfr(havr,param)
        implicit real*8 (a-h, o-z)
        real*8 param(6)
        parameter(PI=3.14159265358979323846264338327d0)

        gamm=param(1)
        d=param(2)
        q=param(3)
        xmag=param(4)

        ssig=d*exp(gamm*xmag)
        sfr_q1=havr+ssig

        if(abs(sfr_q1**(1d0-q)-ssig**(1d0-q)).lt.1d-10) then
           sfr_q23 = log(ssig)-1d0/(q-1d0)
         else      
           sfr_q2=sfr_q1**(1d0-q)*log(sfr_q1)
     &              -ssig**(1d0-q)*log(ssig)
           sfr_q3=(sfr_q1)**(1d0-q)-(ssig)**(1d0-q)
           sfr_q23=sfr_q2/sfr_q3
        endif
        sfr_q4=(ssig+1d0)**(1d0-q)*log(ssig+1d0)-
     &           (ssig)**(1d0-q)*log(ssig)
        sfr_q5=(ssig+1d0)**(1d0-q)-(ssig)**(1d0-q)

        dq_sfr=sfr(havr,param) *(sfr_q4/sfr_q5-sfr_q23)

        return
        end


C********************************************************************
C   pgauss calculates the integral  of the gaussian kernel function
C     with bandwidth w(1) from 0 to r
C
C********************************************************************   
       real*8 function pcauchy(r,w)

       implicit real*8 (a-h,o-z)
       real*8 w(6)
       parameter(PI=3.14159265358979323846264338327d0)
           
       pcauchy=1/2.0d0/PI*(1d0-(1+r**2/w(1)**2)**0.5)
      
       return
       end


C********************************************************************
C   dgauss calculates the value of the gaussian kernel function
C     at r with bandwidth w(1)
C
C********************************************************************   
       real*8 function dcauchy(r,w)

       implicit real*8 (a-h,o-z)
       real*8 w(6)
       parameter(PI=3.14159265358979323846264338327d0)

       dcauchy=PI/2.5/w(1)**2/(1+r**2/w(1)**2)**1.5

       return
       end

C********************************************************************
C   pgauss calculates the integral  of the gaussian kernel function
C     with bandwidth w(1) from 0 to r
C
C********************************************************************   
       real*8 function pgauss(r,w)

       implicit real*8 (a-h,o-z)
       real*8 w(6)
       parameter(PI=3.14159265358979323846264338327d0)

       pgauss=1/2d0/PI*(1d0-exp(-r**2/2d0/w(1)**2))
      
       return
       end

C********************************************************************
C   dgauss calculates the value of the gaussian kernel function
C     at r with bandwidth w(1)
C
C********************************************************************   
       real*8 function dgauss(r,w)

       implicit real*8 (a-h,o-z)
       real*8 w(6)
       parameter(PI=3.14159265358979323846264338327d0)

       dgauss=1/2d0/PI/w(1)**2*exp(-r**2/2d0/w(1)**2)
       return
       end

C****************************************************************
C     Function to calculate the conditional function value at
C      each the time of each event
c     Input:
c       i: sequence no of events
c       b(8): square root of the paratemers
c     Output:
c       fv1: conditional function value 
c       g1: gradient of the conditional function
c***************************************************************
      subroutine xlamb(i,b, fv1, g1)
        
      implicit real*8 (a-h,o-z)
      real*8 b(9), g1(50)
      parameter(PI=3.14159265358979323846264338327d0)
      include 'common1.inc'
      xmu=b(1)**2
      a2=b(2)**2
      c=b(3)**2
      alfa=b(4)**2
      p=b(5)**2
      d=b(6)**2
      q=b(7)**2
      gamma=b(8)**2
      eta = b(9)**2

      s=xmu*zbkgd(i)
      sg1=zbkgd(i)
      sg2=0d0
      sg3=0d0
      sg4=0d0
      sg5=0d0
      sg6=0d0
      sg7=0d0
      sg8=0d0
      sg9=0d0

      if(i.gt.1)then
       do j=1, i-1
         delt=zz(i)-zz(j)
         pr1=exp(alfa*zmg(j))
         pr2=(p-1)/c*(1d0+delt/c)**(-p)
         ssig=d*exp(gamma*zmg(j))

         havd2=havad(xx(i),yy(i),xx(j),yy(j))
         bbb1=ssig+havd2
         bbb2=((ssig+1d0)**(1d0-q)-ssig**(1d0-q))
         bbb3=((ssig+1d0)**(-q)-ssig**(-q))
         bbb4=((ssig+1d0)**(1d0-q)*log(ssig+1d0)-ssig**(1-q))*log(ssig)
         pr3=(1d0-q)/4d0/PI/bbb2*bbb1**(-q)
         
         tmp1 = zdp(j)/depmax
         tmp2 = 1d0-tmp1
          
         tmp = zdp(i)/depmax
         if(zdp(i).eq.0)tmp =0.5/depmax
         if(zdp(i).ge.depmax-1e-8)tmp =1d0-0.5d0/depmax

         pr4 = dbeta(tmp,eta*tmp1+1d0,eta*tmp2+1d0)/depmax
       
         s=s+a2*pr1*pr2*pr3*pr4
         sg2=sg2+pr1*pr2*pr3*pr4

         pr2_c=pr2*(-1d0/c-p/(c+delt)+p/c)
         sg3=sg3+a2*pr1*pr2_c*pr3*pr4

         pr1_alfa=pr1*zmg(j)
         sg4=sg4+a2*pr1_alfa*pr2*pr3*pr4

         pr2_p=pr2*(1d0/(p-1d0)-log(1d0+delt/c))
         sg5=sg5+a2*pr1*pr2_p*pr3*pr4

         pr3_d=pr3*((q-1d0)*bbb3/bbb2-q/bbb1)*ssig/d
         sg6=sg6+a2*pr1*pr2*pr3_d*pr4

         pr3_q= pr3*(1.0/(q-1d0)+bbb4/bbb2-log(bbb1))
         sg7=sg7+a2*pr1*pr2*pr3_q*pr4
        
         pr3_gamma=pr3_d*d*zmg(j)
         sg8=sg8+a2*pr1*pr2*pr3_gamma*pr4

         pr4_ln_eta = tmp1*log(tmp)+tmp2*log(1d0-tmp)
     &         - tmp1*psi(eta*tmp1+1d0) - tmp2*psi(eta*tmp2+1d0)
     &         + psi(eta+2d0)
         sg9=sg9+a2*pr1*pr2*pr3*pr4*pr4_ln_eta


       enddo
      endif
      fv1=s
      g1(1)=sg1*2*b(1)
      g1(2)=sg2*2*b(2)
      g1(3)=sg3*2*b(3)
      g1(4)=sg4*2*b(4)
      g1(5)=sg5*2*b(5)
      g1(6)=sg6*2*b(6)
      g1(7)=sg7*2*b(7)
      g1(8)=sg8*2*b(8)
      g1(9)=sg9*2*b(9)

      return
      end


C*************************************************************
C   Function to calculate integral of the contribution of a
C       previous events 
c   Input:
c       i: sequence no of events
c       b(8): square root of the paratemers
c   Output:
c       fv
C       h
C*************************************************************
        subroutine  xint(i,b,fv,h)
        implicit real*8 (a-h,o-z)
        include 'common1.inc'
        real*8 b(9),h(9)
        real*8 w(6)
        external sfr,dgamma_sfr,dd_sfr, dq_sfr

        xmu=b(1)**2
        a2=b(2)**2
        c=b(3)**2
        alfa=b(4)**2
        p=b(5)**2
        d=b(6)**2
        q=b(7)**2
        gamma=b(8)**2

        gi=0d0
        gic=0d0
        gip=0d0

        if(zz(i).gt.tstart)then
           ttemp=tz-zz(i)
           gi=1d0-(1d0+ttemp/c)**(1-p)           
           gic=-(1d0-gi)*(1d0-p)*(1d0/(c+ttemp)-1d0/c) 
           gip=-(1d0-gi)*(log(c)-log(c+ttemp))
        else
           ttemp2=tz-zz(i)
           ttemp1=tstart-zz(i)
           gi2=1-(1d0+ttemp2/c)**(1-p)           
           gic2=-(1d0-gi2)*(1d0-p)*(1d0/(c+ttemp2)-1d0/c) 
           gip2=-(1d0-gi2)*(log(c)-log(c+ttemp2))
           gi1=1-(1d0+ttemp1/c)**(1-p)           
           gic1=-(1d0-gi1)*(1d0-p)*(1d0/(c+ttemp1)-1d0/c) 
           gip1=-(1d0-gi1)*(log(c)-log(c+ttemp1))

           gi=gi2-gi1
           gic=gic2-gic1
           gip=gip2-gip1     
        endif

        w(1)=gamma
        w(2)=d
        w(3)=q
        w(4)=zmg(i)

        if(ibndtyp.eq.1)then
          si=1.0d0
          sid=0.0d0
          siq=0.0d0
          sigamma=0.0d0
        end if 
        if(ibndtyp.eq.2)then
          si=srectint(sfr,w,tx,ty,xx(i),yy(i)) 
          if(indat(i).eq.1)si= si+1d0
          sid=srectint(dd_sfr,w,tx,ty,xx(i),yy(i))
          siq=srectint(dq_sfr,w,tx,ty,xx(i),yy(i))
          sigamma=sid*d*zmg(i)
        end if 
        if(ibndtyp.eq.3)then
          si=spolyint(sfr,w,tx,ty,npoly,xx(i),yy(i))
          if(indat(i).eq.1)si= si+1d0
          sid=spolyint(dd_sfr,w,tx,ty,npoly,xx(i),yy(i)) 
          siq=spolyint(dq_sfr,w,tx,ty,npoly,xx(i),yy(i))
          sigamma=sid*d*zmg(i)
        end if  

        sk=a2*exp(alfa*zmg(i))

        fv=sk*gi*si

        h(1)=0d0
        h(2)=exp(alfa*zmg(i))*gi*si*2d0*b(2)
        h(3)=sk*gic*si*2*b(3)
        h(4)=sk*gi*si*zmg(i)*2*b(4)
        h(5)=sk*gip*si*2*b(5)
        h(6)=sk*gi*sid*2*b(6)
        h(7)=sk*gi*siq*2*b(7)
        h(8)=sk*gi*sigamma*2*b(8)
        h(9)=0d0

        return
        end

C*******************************************************************
C       fr(r)=\int f(r)r dr
c*******************************************************************
              
        real*8 function fr(r,w)

        implicit real*8 (a-h, o-z)

        real*8 w(6)
        parameter(PI=3.14159265358979323846264338327d0)

        gam=w(1)
        d=w(2)
        q=w(3)
        xmag=w(4)
        
        ssig=d*exp(gam*xmag)        
        

        fr=(1d0-(1+r*r/ssig)**(1d0-q))/pi/2d0
        
c        print *,'fr=',fr
        return
        end


C*******************************************************************
C       dgam_fr(r)=d\int f(r)r dr \over  d\alfa
c*******************************************************************
              
        real*8 function dgam_fr(r,param)

        implicit real*8 (a-h, o-z)
        parameter(PI=3.14159265358979323846264338327d0)

        real*8 param(6)

        gam=param(1)
        d=param(2)
        q=param(3)
        xmag=param(4)
        
        ssig=d*exp(gam*xmag)

        dgam_fr=(-q+1)*(1+r*r/ssig)**(-q)*xmag*r*r/ssig/2d0/pi
        
        return
        end

C*******************************************************************
C       dd_fr(r)=d\int f(r)r dr \over  d\d
c*******************************************************************
              
        real*8 function dd_fr(r,param)

        implicit real*8 (a-h, o-z)

        real*8 param(6)
        parameter(PI=3.14159265358979323846264338327d0)

        gam=param(1)
        d=param(2)
        q=param(3)
        xmag=param(4)
        
        ssig=d*exp(gam*xmag)

        dd_fr=(-q+1)*(1+r*r/ssig)**(-q)/d*r*r/ssig/2d0/pi
        
      
        return
        end

C*******************************************************************
C       dq_fr(r)=d\int f(r)r dr \over  dq
c*******************************************************************
              
        real*8 function dq_fr(r,param)

        implicit real*8 (a-h, o-z)
        parameter(PI=3.14159265358979323846264338327d0)

        real*8 param(6)
        
        gam=param(1)
        d=param(2)
        q=param(3)
        xmag=param(4)
       
C        write(*,*) a2, alfa, d, xmag

        ssig=d*exp(gam*xmag)        

        dq_fr=(1d0+r*r/ssig)**(-q+1)*log(1d0+r*r/ssig)/pi/2d0
          
        return
        end

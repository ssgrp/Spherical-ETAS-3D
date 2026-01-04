       module mod_basicfuns
  
       use mod_kinds, only:dp 
       implicit none
       private
       public::pfisher,dfisher,sfr,dgamma_sfr,dd_sfr,dq_sfr,  &
            pcauchy,dcauchy,pgauss,dgauss,fr,dgam_fr,dq_fr,dd_fr, &
            xlamb, xint, &
            iasign
       contains

!C********************************************************************
!C   pfisher calculates the integral of the von Mises-Fisher 
!C   kernel function
!C     with bandwidth w(1) from 0 to r
!c   Input:      
!C     havr: Haversine of r
!c   Output:
!c     pfisher: Cumulative probability at r    
!C********************************************************************   
       real(dp) function pfisher(havr,w)

!c       implicit real*8 (a-h,o-z)
         use mod_kinds, only:dp
         use mod_params, only:pi
       implicit none
       
       real(dp)::w(6),havr

!       include 'param.inc'
       real(dp)::p1fisher_p1,p2fisher_p2
       
!c       parameter(PI=3.14159265358979323846264338327d0)

       p1fisher_p1=1d0-exp(-havr/2d0/w(1)**2)
       p2fisher_p2=1d0-exp(-0.5d0/w(1)**2)

       pfisher=0.5d0/PI*p1fisher_p1/p2fisher_p2

       return
       end function pfisher

!C********************************************************************



!C*******************************************************************
!C     sfr(r)=\int f(r)r dr
!C    Input: havr haversine of r      
!c*******************************************************************
              
        real(dp) function sfr(havr,w)
!c       implicit real*8 (a-h, o-z)
        use mod_kinds, only:dp 
        use mod_params, only:pi
        implicit none

        real(dp):: w(6),havr

        real(dp)::gamm,d,q,xmag, ssig, sfr_1,sfr_2,sfr_3
        real(dp)::sfr_23
        
!c        real*8 w(6),havr
!c        parameter(PI=3.14159265358979323846264338327d0)

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
        end function sfr

!C********************************************************************
!C   dfisher calculates the value of the von Mises-Fisher kernel function
!C     at r with bandwidth w(1)
!c   Input:      
!C     havr: Haversine of r
!c   Output:
!c     dfisher: density at r    
!C
!C********************************************************************   
       real(dp) function dfisher(havr,w)

!c     implicit real*8 (a-h,o-z)
       use mod_kinds, only:dp 
       use mod_params, only: pi
       implicit none
       real(dp)::w(6), havr
!c       parameter (PI=3.14159265358979323846264338327)
        
!        include 'param.inc'
        real(dp):: d1fisher_d1,d2fisher_d2

       d1fisher_d1=exp(-(havr/2d0/w(1)**2))
       d2fisher_d2=1d0-exp(-0.5d0/w(1)**2)

       dfisher=1/8d0/PI/w(1)**2*d1fisher_d1/d2fisher_d2

       return
       end function dfisher

!C*******************************************************************
!C       dgamma_sfr(r)=d\int f(r)r dr \over  d\gamma
!c*******************************************************************

        real(dp) function  dgamma_sfr(havr,param)

          use mod_kinds, only:dp
!          use mod_params, only: pi
        implicit none
      
!c        implicit real*8 (a-h, o-z)
        real(dp)::param(6), havr
        
!c        parameter(PI=3.14159265358979323846264338327d0)
!        include 'param.inc'
        real(dp)::gamm,d,q,xmag,ssig,sfr_g1,sfr_g2,sfr_g3
!c        real(dp),external::sfr

        gamm=param(1)
        d=param(2)
        q=param(3)
        xmag=param(4)

        ssig=d*exp(gamm*xmag)
        sfr_g1=havr+ssig
        if(abs(sfr_g1**(1d0-q)-ssig**(1d0-q)).lt.1d-10) then
           sfr_g2=q/(q-1d0)/ssig     
         else
           sfr_g2=(sfr_g1**(-q)-ssig**(-q))/(sfr_g1**(1d0-q) &
                  -ssig**(1d0-q))
        endif
        sfr_g3=((ssig+1d0)**(-q)-(ssig)**(-q))/((ssig+1d0)**(1d0-q) &
               -(ssig)**(1d0-q))

        dgamma_sfr=sfr(havr,param)*(1d0-q)*(sfr_g2-sfr_g3)*xmag*ssig

        return
        end function dgamma_sfr

!C*******************************************************************
!C       dd_sfr(r)=d\int f(r)r dr \over  d\d
!c*******************************************************************
              
        real(dp) function dd_sfr(havr,param)
!c        implicit real*8 (a-h, o-z)
          use mod_kinds, only:dp
!          use mod_params, only: pi

        implicit none
        real(dp):: param(6),havr
        
!c        parameter(PI=3.14159265358979323846264338327d0)
!        include 'param.inc'
        real(dp)::gamm,d,q,xmag,ssig,sfr_d1,sfr_d2,sfr_d3
!c        real(dp),external::sfr

        gamm=param(1)
        d=param(2)
        q=param(3)
        xmag=param(4)

        ssig=d*exp(gamm*xmag)        
        sfr_d1=havr+ssig
        if(abs(sfr_d1**(1d0-q)-ssig**(1d0-q)).lt.1d-10) then
           sfr_d2 = q/(q-1d0)/ssig
         else
          sfr_d2=(sfr_d1**(-q)-ssig**(-q))/(sfr_d1**(1d0-q) &
             -ssig**(1d0-q))
        endif
       
        sfr_d3=((ssig+1d0)**(-q)-(ssig)**(-q))/((ssig+1d0)**(1-q) &
              -(ssig)**(1-q))

        dd_sfr=sfr(havr,param)*(1d0-q)*(sfr_d2-sfr_d3)*ssig/d

        return
        end function dd_sfr

!C*******************************************************************
!C       dq_sfr(r)=d\int sf(r)r dr \over  d\q
!c*******************************************************************
              
        real(dp) function dq_sfr(havr,param)

!c        implicit real*8 (a-h, o-z)
          use mod_kinds, only:dp
          use mod_params,only:dp
        implicit none
        
        real(dp)::param(6), havr
        
!c        parameter(PI=3.14159265358979323846264338327d0)
!        include 'param.inc'
        real(dp)::gamm,d,q,xmag,ssig,sfr_q1,sfr_q2,sfr_q3
        real(dp)::sfr_q4,sfr_q5,sfr_q23

        gamm=param(1)
        d=param(2)
        q=param(3)
        xmag=param(4)

        ssig=d*exp(gamm*xmag)
        sfr_q1=havr+ssig

        if(abs(sfr_q1**(1d0-q)-ssig**(1d0-q)).lt.1d-10) then
           sfr_q23 = log(ssig)-1d0/(q-1d0)
         else      
           sfr_q2=sfr_q1**(1d0-q)*log(sfr_q1) &
                   -ssig**(1d0-q)*log(ssig)
           sfr_q3=(sfr_q1)**(1d0-q)-(ssig)**(1d0-q)
           sfr_q23=sfr_q2/sfr_q3
        endif
        sfr_q4=(ssig+1d0)**(1d0-q)*log(ssig+1d0)- &
                (ssig)**(1d0-q)*log(ssig)
        sfr_q5=(ssig+1d0)**(1d0-q)-(ssig)**(1d0-q)

        dq_sfr=sfr(havr,param) *(sfr_q4/sfr_q5-sfr_q23)

        return
        end function dq_sfr


!C********************************************************************
!C   pgauss calculates the integral  of the gaussian kernel function
!C     with bandwidth w(1) from 0 to r
!C
!C********************************************************************   
       real(dp) function pcauchy(r,w)

!c       implicit real*8 (a-h,o-z)
       use mod_kinds, only:dp 
       use mod_params, only: pi
       implicit none
       real(dp):: w(6),r
       
!c       parameter(PI=3.14159265358979323846264338327d0)
!        include 'param.inc'
           
       pcauchy=1/2.0d0/PI*(1d0-(1+r**2/w(1)**2)**0.5)
      
       return
       end function pcauchy


!C********************************************************************
!C   dgauss calculates the value of the gaussian kernel function
!C     at r with bandwidth w(1)
!C
!C********************************************************************   
       real(dp) function dcauchy(r,w)

!c       implicit real*8 (a-h,o-z)
        use mod_kinds, only:dp 
        use mod_params, only: pi
        implicit none
        real(dp):: w(6),r
!c       parameter(PI=3.14159265358979323846264338327d0)
!        include 'param.inc'

       dcauchy=PI/2.5/w(1)**2/(1+r**2/w(1)**2)**1.5

       return
       end function dcauchy

!C********************************************************************
!C   pgauss calculates the integral  of the gaussian kernel function
!C     with bandwidth w(1) from 0 to r
!C
!C********************************************************************   
       real(dp) function pgauss(r,w)

!c     implicit real*8 (a-h,o-z)
       use mod_kinds, only:dp 
       use mod_params, only: pi
       implicit none
       
       real(dp)::w(6), r
!c       parameter(PI=3.14159265358979323846264338327d0)
!        include 'param.inc'
       
        pgauss=1/2d0/PI*(1d0-exp(-r**2/2d0/w(1)**2))
      
       return
       end function pgauss

!C********************************************************************
!C   dgauss calculates the value of the gaussian kernel function
!C     at r with bandwidth w(1)
!C
!C********************************************************************   
       real(dp) function dgauss(r,w)

!c     implicit real*8 (a-h,o-z)
        use mod_kinds, only:dp 
          use mod_params, only: pi
       implicit none
       
       real(dp):: w(6),r
       
!        include 'param.inc'
!c       parameter(PI=3.14159265358979323846264338327d0)

       
       dgauss=1/2d0/PI/w(1)**2*exp(-r**2/2d0/w(1)**2)
       return
       end function dgauss

!C****************************************************************
!C     Function to calculate the conditional function value at
!C      each the time of each event
!c     Input:
!c       i: sequence no of events
!c       b(8): square root of the paratemers
!c     Output:
!c       fv1: conditional function value 
!c       g1: gradient of the conditional function
!c***************************************************************
      subroutine xlamb(i,b, fv1, g1)
        
!c     implicit real*8 (a-h,o-z)
        use mod_kinds, only:dp
        use mod_geomsphere, only: havad
        use mod_specialfun,only:psi,dbeta
        use mod_params, only: pi
 
        use mod_eq_data
        use mod_med_result
        use mod_mpi_state
        use mod_domain
         
      implicit none
      real(dp):: b(9), g1(50),fv1
      integer::i
      
!c      parameter(PI=3.14159265358979323846264338327d0)
!      include 'common1.inc'
!      include 'param.inc'

      real(dp)::xmu,a2,c,alfa,p,d,q,gamma,eta
      real(dp)::s,sg1,sg2,sg3,sg4,sg5,sg6,sg7,sg8,sg9
      real(dp)::delt,havd2,bbb1,bbb2,bbb3,bbb4
      real(dp)::pr1,pr1_alfa,pr2,pr2_c,pr2_p,pr3,pr3_d,pr3_q
      real(dp)::pr3_gamma,pr4,pr4_ln_eta,ssig,tmp,tmp1,tmp2
!c      real(dp),external::psi,dbeta
      integer::j
      
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

         pr4_ln_eta = tmp1*log(tmp)+tmp2*log(1d0-tmp) &
              - tmp1*psi(eta*tmp1+1d0) - tmp2*psi(eta*tmp2+1d0) &
              + psi(eta+2d0)
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

      lambdas(i)=fv1
      
      return
      
      end subroutine xlamb


!C*************************************************************
!C   Function to calculate integral of the contribution of a
!C       previous events 
!c   Input:
!c       i: sequence no of events
!c       b(8): square root of the paratemers
!c   Output:
!c       fv
!C       h
!C*************************************************************
        subroutine  xint(i,b,fv,h)
!c     implicit real*8 (a-h,o-z)
        use mod_kinds, only:dp 
        use mod_fintegral,only:srectint,spolyint
 
        use mod_eq_data
        use mod_med_result
        use mod_mpi_state

        use mod_domain
        implicit none
        
        real(dp)::b(9),h(9),fv
        integer:: i
        
!c        real*8 w(6)
!c        real(8),external:: sfr,dgamma_sfr,dd_sfr, dq_sfr

!        include 'common1.inc'
        
!c     real*8 b(9),h(9)
        real(dp)::xmu,a2,c,alfa,p,d,q,gamma
        real(dp)::gi, gic, gip,ttemp, gi1,gi2,gic1,gic2,gip1,gip2,si
        real(dp)::sid,sigamma,siq,sk, ttemp1,ttemp2

        
        real(dp):: w(6)
!c        real(dp),external::sfr,dgamma_sfr,dd_sfr, dq_sfr,srectint
!         real(dp),external::srectint
!         real(dp),external::spolyint
       
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
        end subroutine xint

!C*******************************************************************
!C       fr(r)=\int f(r)r dr
!c*******************************************************************
              
        real(dp) function fr(r,w)

!c        implicit real*8 (a-h, o-z)
         use mod_kinds, only:dp 
          use mod_params, only: pi
         implicit none
        real(dp):: r, w(6)
!        include 'param.inc'
        real(dp)::gam,d,q,xmag,ssig
        
!c        parameter(PI=3.14159265358979323846264338327d0)

        gam=w(1)
        d=w(2)
        q=w(3)
        xmag=w(4)
        
        ssig=d*exp(gam*xmag)        
        

        fr=(1d0-(1+r*r/ssig)**(1d0-q))/pi/2d0
        
!c        print *,'fr=',fr
        return
        end function fr


!C*******************************************************************
!C       dgam_fr(r)=d\int f(r)r dr \over  d\alfa
!c*******************************************************************
              
        real(dp) function dgam_fr(r,param)

        use mod_kinds, only:dp 
        use mod_params, only: pi

        implicit none
        real(dp):: r, param(6)
!        include 'param.inc'
        real(dp)::gam,d,q,xmag,ssig
!c       implicit real*8 (a-h, o-z)
!c        parameter(PI=3.14159265358979323846264338327d0)

!c        real*8 param(6)

        gam=param(1)
        d=param(2)
        q=param(3)
        xmag=param(4)
        
        ssig=d*exp(gam*xmag)

        dgam_fr=(-q+1)*(1+r*r/ssig)**(-q)*xmag*r*r/ssig/2d0/pi
        
        return
        end function dgam_fr

!C*******************************************************************
!C       dd_fr(r)=d\int f(r)r dr \over  d\d
!c*******************************************************************
              
        real(dp) function dd_fr(r,param)

        use mod_kinds, only:dp 
          use mod_params, only: pi

        implicit none
        real(dp)::r, param(6)

!        include 'param.inc'
        real(dp)::gam,d,q,xmag,ssig

!c        implicit real*8 (a-h, o-z)
!c
!c        real*8 param(6)
!c        parameter(PI=3.14159265358979323846264338327d0)

        gam=param(1)
        d=param(2)
        q=param(3)
        xmag=param(4)
        
        ssig=d*exp(gam*xmag)

        dd_fr=(-q+1)*(1+r*r/ssig)**(-q)/d*r*r/ssig/2d0/pi
        
      
        return
        end function dd_fr

!C*******************************************************************
!C       dq_fr(r)=d\int f(r)r dr \over  dq
!c*******************************************************************
              
        real(dp) function dq_fr(r,param)

        use mod_kinds, only:dp 
          use mod_params, only: pi
        implicit none
        real(dp):: r, param(6)
!        include 'param.inc'
        real(dp)::gam,d,q,xmag,ssig

!c        implicit real*8 (a-h, o-z)
!c        parameter(PI=3.14159265358979323846264338327d0)
!c
!c        real*8 param(6)
        
        gam=param(1)
        d=param(2)
        q=param(3)
        xmag=param(4)
       
!C        write(*,*) a2, alfa, d, xmag

        ssig=d*exp(gam*xmag)        

        dq_fr=(1d0+r*r/ssig)**(-q+1)*log(1d0+r*r/ssig)/pi/2d0
          
        return
        end function dq_fr

        
!c-------------------------------------------------------
      integer function iasign(im, np, mr)
!c      implicit real*8 (a-h, o-z)
!      use mod_kinds, only:dp 
      implicit none
      integer::im,np,mr

      integer::ia, k
      
      ia = 0 
      k = mod(im-1, 2*np) 
      if(k == mr.or. k ==  2*np-mr-1) ia = 1
      iasign = ia
      
      return 
    end function iasign

      end module mod_basicfuns
      

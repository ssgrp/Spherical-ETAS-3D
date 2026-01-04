CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   The subroutine outback outputs the background rate on a lattice of
C      mx*my grids
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine outrates(file1)
      implicit real*8 (a-h,o-z)
      character*80 file1
      include 'mpif.h'
      include 'common.inc'
      real*8 rtx1, rtx2, rty1, rty2, rtz1, rtz2
      common /range/ rtx1, rtx2, rty1, rty2, rtz1, rtz2
      include 'param.inc'

      real*8 w(10)

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
                  btemp = dbeta(tmp, bwh*zdp(ii)/depmax+1,
     &                 bwh*(1-zdp(ii)/depmax)+1)/depmax 
     
                  s=s+zprob(ii)*dfisher(havr0,w)*btemp
                  s1=s1+dfisher(havr0,w)*btemp
                  
               enddo
               xlamb = 0.0 + xmu * s/(tz-tstart)
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
                 xlamb = xlamb + a2*pr1*pr2*pr3*pr4

               enddo           
               write(34,991)s/(tz-tstart),s1/(tz-tstart), xlamb
            enddo
         enddo
       enddo
       close(34)
      endif
      return 
 991  format(1x,3f20.14)
      end
    
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   The subroutine outprob outputs the probability of each event
C      being a background event or not
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        subroutine outprob(file1)
        implicit real*8 (a-h,o-z)
        character*80 file1
c        include 'mpif.h'
        include 'common.inc'

        if(myrank.eq.2)then
             open(35,file=file1)
             do i=1,nn
                write(35,992)i,zprob(i),zbandw(i),zbkgd(i)
             enddo
             close(35)
        endif
       return 

 990  format(6(1x,i5,f7.4))
 992  format(1x,i8,3f20.10)
      end

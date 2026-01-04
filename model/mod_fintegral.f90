       module mod_fintegral
       use mod_kinds, only:dp 

       implicit none
       private
       public::sfrint,srectint,spolyint
       
       contains

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C  This function calculates the integral within the triangle
!c
!c  Input:
!c       func,para,x1,y1,x2,y2,xc,yc 
!c  Output:
!c       sfrint
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                             
        real(dp) function sfrint(func,para,x1,y1,x2,y2,xc,yc)
        
!c        implicit real*8 (a-h,o-z)
        use mod_kinds, only:dp
        use mod_geomsphere,only:havad,arcdist,hav,beardirt,dirdif,bisect
        implicit none
        
        real(dp),external:: func
       

        real(dp):: para(6)
        real(dp)::x1,y1,x2,y2,xc,yc
        
        real(dp)::r1,r2,r12,d1,d2,d12,theta,r0,f1,f2,f3,ban1,ban2
        real(dp)::bx,by
!c        real(dp),external:: havad,arcdist,dirdif,bisect,hav
        
!c        real*8 xtemp(3),ytemp(3)

        r1=havad(x1,y1,xc,yc)
        r2=havad(x2,y2,xc,yc)
        r12=havad(x1,y1,x2,y2)

        d1=arcdist(x1,y1,xc,yc)
        d2=arcdist(x2,y2,xc,yc)
        d12=arcdist(x1,y1,x2,y2)

        call beardirt(xc,yc,x1,y1,bx,by,ban1)
        call beardirt(xc,yc,x2,y2,bx,by,ban2)
        theta=dirdif(ban1, ban2, 1)

        if(r1+r2.gt.1d-10)then
           r0= bisect(d1,d2,d12)
        else 
           sfrint=0
           return
        endif

        f1=func(r1,para)
        f2=func(hav(r0),para)
        f3=func(r2,para)

        sfrint=(f1/6d0+f2*2d0/3d0+f3/6d0)*(-theta)

        return
        end function sfrint

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C  This function calculates the integral within the rectangle 
!c
!c  Input:
!c       func,para,tx,ty,xc,yc
!c  Output:
!c       srectint
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        real(dp) function srectint(func,para,tx,ty,xc,yc)

        use mod_kinds, only:dp 
        use mod_geomsphere,only:dirdif
        use mod_params,only:pi
        
        implicit none
!c        implicit real*8 (a-h,o-z)
        real(dp):: para(6),tx(2),ty(2), ttx(5), tty(5), stemp(5)
        
        real(dp),external:: func
        
!        include 'param.inc'

        real(dp):: xc,yc,diflat,diflen, diflon,dlat,dlon,sum
        real(dp):: x1, y1, x2, y2
        integer:: i, j, ndiv


        ttx(1)=tx(1)
        ttx(2)=tx(2)
        ttx(3)=tx(2)
        ttx(4)=tx(1)
        ttx(5)=tx(1)
        
        tty(1)=ty(1)
        tty(2)=ty(1)
        tty(3)=ty(2)
        tty(4)=ty(2)
        tty(5)=ty(1)
 
        sum=0
        do j=1, 4, 2
           stemp(j)=0d0
           diflon = dirdif(ttx(j),ttx(j+1), 0)
           diflen = diflon*cos(tty(j))
           ndiv = int(abs(diflen)/(0.05d0/180d0*pi))+2
           dlon = diflon /(dble(ndiv)-1d0)

           do i=1, ndiv-1
             x1=ttx(j)+dlon*dble(i-1)
             if(x1.gt.pi)x1=x1-pi*2.0d0
             x2=ttx(j)+dlon*dble(i)
             if(x2.gt.pi)x2=x2-pi*2.0d0
             y1=tty(j)
             y2=tty(j)
             stemp(j)=stemp(j)+sfrint(func,para,x1,y1,x2,y2,xc,yc)
           enddo  

        enddo
        
        do j=2,4,2
           stemp(j)=0d0
           diflat = (tty(j+1)-tty(j))
           ndiv = int(abs(diflat)/(0.05d0/180d0*pi)+2)
           dlat = diflat/(dble(ndiv)-1d0)
            
          do i=1, ndiv-1
             y1=tty(j)+dlat*dble(i-1)
             y2=tty(j)+dlat*dble(i)
             x1=ttx(j)
             x2=ttx(j)
             stemp(j)=stemp(j)+sfrint(func,para,x1,y1,x2,y2,xc,yc)
          enddo 

        enddo

        srectint=stemp(1)+stemp(2)+stemp(3)+stemp(4)
        return
        end function srectint


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C  This function calculates the integral within the spherical polygon 
!c  
!c  Input:
!c       func,para,xp,yp,np,xc,yc 
!c  Output:
!c       spolyint
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        real(dp) function spolyint(func,para,xp,yp,np,xc,yc)
       
!c        implicit real*8 (a-h,o-z)
        use mod_kinds, only:dp
        use mod_geomsphere,only:dirdif,havad,hav,beardirt
        use mod_params, only: pi

         implicit none
        
         real(dp):: para(6), xp(np+1), yp(np+1),xc,yc
         integer::np
!c        real*8 fd(10000),ttemp

        real(dp),external:: func

!        include 'param.inc'
        real(dp)::ttemp, alfa, antmp,ba1, ba2, bb, bd, bx
        real(dp)::by,da1,da2,delta,dtmp,ssum, sum1
        integer::i,j,ndiv
!        real(dp),external::simpsum
        real(dp)::s_simpson,fi

        ssum=0d0

        xp(np+1)=xp(1)
        yp(np+1)=yp(1)

        do j=1,np
         call beardirt(xc, yc, xp(j), yp(j), bx, by, da1)
         call beardirt(xc, yc, xp(j+1), yp(j+1), bx, by, da2)
         bd = dirdif(da1,da2,1)      
         if(abs(bd).gt.0d0.and.abs(bd).le.pi) then   
           ndiv= int(abs(bd)/(0.005d0/180d0*pi)) + 2
           if(ndiv.lt.7)ndiv=7
           if(ndiv/2*2.eq. ndiv) ndiv = ndiv +1
           delta = bd/dble(ndiv-1)
           call beardirt(xp(j), yp(j), xc,yc, bx, by, ba1)
           call beardirt(xp(j), yp(j), xp(j+1), yp(j+1), bx, by, ba2)
           alfa = abs(dirdif(ba1,ba2,1))
           bb = havad(xp(j), yp(j), xc, yc)

           s_simpson=0.0_dp
           do i=1, ndiv
              antmp = abs(delta)*dble(i-1)
              dtmp = 2d0*sqrt(bb*(1d0-bb))/((1d0-2d0*bb)*cos(antmp)+  &
                 sin(antmp)/tan(alfa))
              ttemp=hav(atan(dtmp))
              fi=func(ttemp,para)
              
              if(i==1.or.i ==ndiv) then
                 s_simpson=s_simpson + fi
              else if (mod(i,2)==0) then
                 s_simpson=s_simpson +4.0_dp*fi
              else   
                 s_simpson=s_simpson +2.0_dp*fi
              endif  
           enddo
           
           sum1 = (s_simpson/3.0_dp)*(-delta)          
!           sum1 =  simpsum(fd, ndiv)*(-delta)
           ssum = ssum + sum1
        else if (abs(bd).gt.0d0) then
!          never be here since -pi < bd  <= pi
           ssum = ssum + func(1d0,para)*(-bd)
         endif

        enddo
                
        spolyint=ssum
        return
        end function spolyint


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      Simpson integral of an array
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
        real(dp) function simpsum(fd, ndiv)
!c        implicit real*8 (a-h,o-z)
        use mod_kinds, only:dp 
        implicit none
        
        integer::ndiv
        real(dp):: fd(ndiv)

        integer:: i
        real(dp)::ssum

        ssum=fd(1)
        do i=2, ndiv-1, 2
           ssum=ssum+fd(i)*4d0
        enddo
        
        do i=3, ndiv-1, 2
           ssum=ssum+fd(i)*2d0
        enddo

        ssum=ssum+fd(ndiv)
        simpsum=ssum/3d0

        return 
        end  function simpsum

        end module mod_fintegral

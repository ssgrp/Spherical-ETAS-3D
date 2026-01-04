CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  This function calculates the integral within the triangle
c
c  Input:
c       func,para,x1,y1,x2,y2,xc,yc 
c  Output:
c       sfrint
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                             
        real*8 function sfrint(func,para,x1,y1,x2,y2,xc,yc)
        
        implicit real*8 (a-h,o-z)
        external func

        real*8 para(6)  
        real*8 xtemp(3),ytemp(3)

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
        end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  This function calculates the integral within the rectangle 
c
c  Input:
c       func,para,tx,ty,xc,yc
c  Output:
c       srectint
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        real*8 function srectint(func,para,tx,ty,xc,yc)

        implicit real*8 (a-h,o-z)
        real*8 para(6),tx(2),ty(2), ttx(5), tty(5), stemp(5)
        external func
        include 'param.inc'

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
           ndiv = abs(diflen)/(0.05d0/180d0*pi)+2
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
           ndiv = abs(diflat)/(0.05d0/180d0*pi)+2
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
        end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  This function calculates the integral within the spherical polygon 
c  
c  Input:
c       func,para,xp,yp,np,xc,yc 
c  Output:
c       spolyint
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        real*8 function spolyint(func,para,xp,yp,np,xc,yc)
       
        implicit real*8 (a-h,o-z)
        include 'param.inc'
        real*8 para(6), xp(np+1), yp(np+1)
        real*8 fd(10000),ttemp
        external func


        ssum=0d0

        xp(np+1)=xp(1)
        yp(np+1)=yp(1)

        do j=1,np
         call beardirt(xc, yc, xp(j), yp(j), bx, by, da1)
         call beardirt(xc, yc, xp(j+1), yp(j+1), bx, by, da2)
         bd = dirdif(da1,da2,1)      
         if(abs(bd).gt.0d0.and.abs(bd).lt.pi) then   
           ndiv= abs(bd)/(0.005d0/180d0*pi) + 2
           if(ndiv.lt.7)ndiv=7
           if(ndiv/2*2.eq. ndiv) ndiv = ndiv +1
           delta = bd/dble(ndiv-1)
           call beardirt(xp(j), yp(j), xc,yc, bx, by, ba1)
           call beardirt(xp(j), yp(j), xp(j+1), yp(j+1), bx, by, ba2)
           alfa = abs(dirdif(ba1,ba2,1))
           bb = havad(xp(j), yp(j), xc, yc)

           do i=1, ndiv
              antmp = abs(delta)*dble(i-1)
              dtmp = 2d0*sqrt(bb*(1d0-bb))/((1d0-2d0*bb)*cos(antmp)+
     &             sin(antmp)/tan(alfa))
              ttemp=hav(atan(dtmp))
              fd(i)=func(ttemp,para)
           enddo
           sum1 = simpsum(fd, ndiv)*(-delta)
           ssum = ssum + sum1
         else if (abs(bd).gt.0d0) then
            ssum = ssum + func(1d0,para)*(-bd)
         endif

        enddo
                
        spolyint=ssum
        return
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      Simpson integral of an array
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
        real*8 function simpsum(fd, ndiv)
        implicit real*8 (a-h,o-z)
        real*8 fd(ndiv)

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
        end
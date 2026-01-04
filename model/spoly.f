
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  This function is sign function
c  returns the sign of the specified number
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
        integer function sgn(x)
        double precision x
          if (x .gt. 0.0d0) then 
              sgn=1
            elseif (x .lt. 0.0d0) then
              sgn=-1
            else 
              sgn=0
          endif
        return
        end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  This subroutine determin longitude xtar lies in the small arc connecting
c    longitudes xsrc1 and xsrc2 along the equator. 
C  Input: xsrc1, xsrc2, xtar (requires that abs(xsrc1-xsrc2)!= 180 )
c  Output:
c     isbetween: 1 yes, 0, overlapping, -1 no.
c     iecode: 0, when isbetween is 1 or -1; 1, xsrc1 and xtar overlap;
c             2, xsrc2 and xtar overlap; 3, xsrc1, xsrc2, and xtar overlap
c      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine between(xsrc1,xsrc2,xtar, isbetween,iecode)

        implicit real*8 (a-h, o-z)
        real*8 xsrc1,xsrc2, xtar

        isbetween = -1
        iecode=0

        dif = abs(dirdif(xsrc1, xsrc2, 1))
        dif1 =abs(dirdif(xsrc1, xtar, 1))
        dif2 =abs(dirdif(xtar, xsrc2, 1))
        
        if(abs(dif1+dif2-dif).lt.1d-10)then
           isbetween = 1
           if(dif1.lt.1d-10) then
              isbetween =0
              iecode=1
           endif
           if(dif2.lt.1d-10) then
              isbetween =0
              iecode=2
           endif
           if(dif.lt.1d-10) then
              isbetween =0
              iecode=3
           endif
        endif

        return 
        end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  This subroutine is used to obtain crossprod
C  Input: x1,y1,z1,x2,y2,z2
c  Output: x3,y3,z3
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine crossprod(x1,y1,z1,x2,y2,z2, x3,y3,z3)

      implicit real*8 (a-h, o-z)

      x3 = y1 * z2 - y2 * z1
      y3=  z1 * x2 - z2 * x1
      z3 = x1 * y2 - x2 * y1
      return
      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  This function determines whether two points overlap
C  isovlap: 1 overlapping, 0 no
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function isovlap(x1,y1,z1,x2,y2,z2)
      
      implicit real*8 (a-h, o-z)

      isovlap = 0

      if(abs(x2-x1)+abs(y2-y1)+abs(z2-z1).lt.1d-10) isovlap=1

      return
      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  This function determines whether point (xx,yy,zz) is on the arc segment
C  connecting (x1,y1,z1) and (x2,y2,z2)
C  isinarcseg: 1 on the line 
C              0 no 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function isinarcseg(x1,y1,z1,x2,y2,z2,xx,yy,zz)
      implicit real*8 (a-h, o-z)

      include 'param.inc'

      isinarcseg = 0

      call crossprod(x1,y1,z1,x2,y2,z2,xp,yp,zp)

      if(abs(xx*xp+yy*yp+zz*zp).lt.1.0d-10)then
         
         a1 = x1**2 + y1**2 +z1**2
         a2 = x2**2 + y2**2 +z2**2
         ap = xx**2 + yy**2 +zz**2  

         cosa1p = (x1*xx+y1*yy+z1*zz )/sqrt(a1*ap)
         cosa2p = (x2*xx+y2*yy+z2*zz )/sqrt(a2*ap)
         cosa12 = (x1*x2+y1*y2+z1*z2 )/sqrt(a1*a2)

         if(cosa1p.le.-1.0d0)then
            an1p = pi
         elseif (cosa1p.ge.1d0) then
            an1p = 0.d0
         else
            an1p = acos(cosa1p)
         endif
         
         if(cosa2p.le.-1.0d0)then
            an2p = pi
         elseif (cosa2p.ge.1d0) then
            an2p = 0.d0
         else
            an2p = acos(cosa2p)
         endif

         if(cosa12.le.-1.0d0)then
            an12 = pi
         elseif (cosa12.ge.1d0)then
            an12 = 0.d0
         else
            an12 = acos(cosa12)
         endif

         if(abs(an12-an1p-an2p).lt.1.0d-9) then
             isinarcseg =1
         endif
      endif

      return
      end


      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  This subroutine converts latitude and longitude into three-dimensional 
c  Cartesian coordinates
c  Input: xlon, xlat
c  Output: x, y, z
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
        subroutine spc2dc(xlon,xlat, x,y,z)

        implicit real*8 (a-h, o-z)
 
        x = cos(xlon) * cos(xlat)
        y = sin(xlon) * cos(xlat)
        z = sin(xlat)

        return 
        end



ccccCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc 
c   This subroutine evaluate whether the great Arc segment 
c   connecting (plon1, plat1) and (plon2, plat2) intersects
c   with the longitude arc from (xlon1, xlat1) to (xlon2, xlat2)
c 
c    Input:
c       plon1, plat1, plon2, plat2: vertices of the great Arc segment
c                                        
c       xlon1, xlat1, xlon2, xlat2: the other great arc segment
c      
c    Output:
C        
c     iflag1: 1 p1-p2 and x1-x2 intersect; 
c             0 overlapping;
c            -1 no intersection; 
c     iflag2: Only valid when flag1 is 0 or -1, connection with the other segment
c     iflag3: Only valid when flag1 is 0 or -1, overlapping other points
c                  
c     iflag1, flag2(x2 x1 p2 p1)  flag3 (x2,x1,p2,p1):
c      1, 0, 0:   p1-p2 intersects x1-x2 at (vlon,vlat)
c     -1, 0, 0:   p1-p2 and x1-x2 have no intersection
c      
c      0,15, 0:   p1, p2, x1, and x2 overlap at the same location or p1(x1)-p2(x2) or p1(x2)-p2(x1)
c     -1,15, 0:   p1 overlaps p2, x1 overlap x2 at different locations
c        
c     -1, 3, 0:   p1 overlaps p2, but not on x1-x2
c      0, 7, 0:   p1 overlaps p2 at x1
c      0,11, 0:   p1 overlaps p2 at x2
c      0, 3, 3:   p1 overlaps p2 at a point in the x1-x2 segment
c
c     -1,12, 0:   x1 overlaps x2, but not on p1-p2
c      0,13, 0:   x1 overlaps x2 at p1
c      0,14, 0:   x1 overlaps x2 at p2
c      0,12,10:   x1 overlaps x2 at a point in p1-p2
c      
c      0, 5, 0:   p1 overlaps x1, neither p2 is in x1-x2 nor x2 in p1-p2
c      0, 5, 2:   p1 overlaps x1, and p2 is in x1-x2
c      0, 5, 8:   p1 overlaps x1, and x2 is in p1-p2
c      
c      0, 9, 0:   p1 overlaps x2, neither p2 is in x1-x2 nor x1 in p1-p2
c      0, 9, 2:   p1 overlaps x2, and p2 is in x1-x2
c      0, 9, 4:   p1 overlaps x2, and x1 is in p1-p2
c      0,10, 0:   p2 overlaps x2, neither p1 is in x1-x2 nor x1 in p1-p2
c      0,10, 1:   p2 overlaps x2, and p1 is in x1-x2
c      0,10, 4:   p2 overlaps x2, and x1 is in p1-p2
c           
c      0, 6, 0:   p2 overlaps x1, neither p1 is in x1-x2 nor x2 in p1-p2
c      0, 6, 1:   p2 overlaps x1, and p1 is in x1-x2
c      0, 6, 8:   p2 overlaps x1, and x2 is in p1-p2
c
c      0, 0, 1:   p1 is on x1-x2
c      0, 0, 2:   p2 is on x1-x2
c      0, 0, 3:   both p1 and p2 are on x1-x2  
c      0, 0, 4:   x1 is on p1-p2
c      0, 0, 8:   x2 is on p1-p2
c      0, 0,12:   both x1 and x2 are on p1-p2
c      
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        subroutine inter2Arc(plon1, plat1, plon2,plat2,  
     &     xlon1, xlat1, xlon2, xlat2, iflag1,
     &       iflag2, iflag3)

        implicit real*8 (a-h, o-z)
        
        iflag1 = -1
        iflag2 = 0
        iflag3 = 0
        
        call spc2dc(plon1, plat1, px1, py1, pz1)
        call spc2dc(plon2, plat2, px2, py2, pz2)
        call spc2dc(xlon1, xlat1, x1, y1, z1)
        call spc2dc(xlon2, xlat2, x2, y2, z2)
        
cccc Checking overlapping points and setting iflag2
        ivp1p2 = isovlap(px1,py1,pz1,px2, py2, pz2)

        if(ivp1p2.eq.1)iflag2 = ior(iflag2, 3)
        
        ivx1x2 = isovlap(x1,y1,z1,x2,y2,z2)
        if(ivx1x2.eq.1)iflag2 = ior(iflag2, 12)
        
        ivp1x1 = isovlap(px1,py1,pz1,x1,y1,z1)
        if(ivp1x1.eq.1)iflag2 = ior(iflag2, 5)
        
        ivp1x2 = isovlap(px1,py1,pz1,x2,y2,z2)
        if(ivp1x2.eq.1)iflag2 = ior(iflag2,9)
        
        ivp2x1 = isovlap(px2,py2,pz2,x1,y1,z1)
        if(ivp2x1.eq.1)iflag2 = ior(iflag2, 6)
        
        ivp2x2 = isovlap(px2,py2,pz2,x2,y2,z2)
        if(ivp2x2.eq.1)iflag2 = ior(iflag2,10)
        
ccc   Return if there are 2 pairs, a tripple or a quadruple of points overlapping        
        if(iflag2.eq.15) then
          if(ivp1x1.eq.1.or.ivp1x2.eq.1) then
           iflag1 = 0
           iflag3 = 0
           return
          else
           iflag1 = -1
           iflag3 = 0
           return
          endif
        endif
      
        if(iflag2.eq.14 .or. iflag2.eq.13 .or. iflag2.eq.10
     &       .or. iflag2.eq.7) then
           iflag1 = 0
           iflag3 = 0
           return 
        endif
        
ccc   Check whether each point is on the other segement.
        
        insx1 = isinarcseg(px1,py1,pz1,px2,py2,pz2,x1,y1,z1)
        insx2 = isinarcseg(px1,py1,pz1,px2,py2,pz2,x2,y2,z2)
        insp1 = isinarcseg(x1,y1,z1,x2,y2,z2, px1,py1,pz1)
        insp2 = isinarcseg(x1,y1,z1,x2,y2,z2, px2,py2,pz2)

ccc   Setting iflag3 if on the other segment but not overlaping the ends
        
        if(insp1.eq.1.and.ivp1x1.eq.0.and.
     &     ivp1x2.eq.0) iflag3 = ior(iflag3, 1)
      
        if(insp2.eq.1.and.ivp2x1.eq.0.and.
     &     ivp2x2.eq.0) iflag3 = ior(iflag3, 2)
        if(insx1.eq.1.and.ivp1x1.eq.0.and.
     &     ivp2x1.eq.0) iflag3 = ior(iflag3, 4)
        if(insx2.eq.1.and.ivp1x2.eq.0.and.
     &     ivp2x2.eq.0) iflag3 = ior(iflag3, 8)

cc              
        if(iflag2.eq.3 .or.iflag2.eq.12) then
           if(iflag3.eq.0) then
               iflag1 = -1
               return
           else
               iflag1 = 0
               return
           endif
        endif

        if(iflag2.ne.0 .or. iflag3.ne.0) then
           iflag1 = 0
           return
        endif

        call crossprod(px1,py1,pz1, px2,py2, pz2, ppx,ppy,ppz)
        call crossprod(x1,y1,z1,x2,y2,z2,xxx,xxy,xxz)       
        call crossprod(ppx,ppy,ppz, xxx,xxy,xxz, vx,vy,vz)

        is1 = isinarcseg(px1,py1,pz1,px2,py2,pz2,vx,vy,vz)
        is2 = isinarcseg(x1,y1,z1,x2,y2,z2,vx,vy,vz)

        is3 = isinarcseg(px1,py1,pz1,px2,py2,pz2, -vx,-vy,-vz)
        is4 = isinarcseg(x1,y1,z1,x2,y2,z2,-vx,-vy,-vz)

        if((is1.eq.1.and.is2.eq.1).or.(is3.eq.1.and.is4.eq.1))then
           iflag1 = 1
        else
           iflag1 = -1
        endif

        return 
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  This subroutine determines whether the point is inside the rectangle
c   Input:px1, px2, py1, py2, xx, yy, m
C     (lon:-180~0,0~180  lat:-90~0,0~90)
C     px1 left to px2, py1 < py2: Normal inclusion
c     px1 right to px2, py1 > py2: Normal exclusion
c     px1 left to px2, py1 > py2: Double triangles with two poles.
c     px1 right to px2, py1 < py2: stripe opposite to normal inclusion.
c   Output: iflag: -1, outside; 1 inside
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
       subroutine srecse(px1, px2, py1, py2, xx, yy, m, iflag)

       implicit real *8 (a-h, o-z)
      
       real*8 xx(1:m),yy(1:m),aa12(1:m),aa111(1:m),aa112(1:m),aa113(1:m)
       integer sgn
       integer iflag(1:m),iflag11(1:m),iflag12(1:m), m

       do i=1, m
       write(18,*) px1, px2, py1, py2, xx(i), yy(i)   
          iflag(i)=-1
          if(py1.lt.py2)then             
             if(abs(dirdif(px1,xx(i),0)+dirdif(xx(i),px2,0)-dirdif(px1,
     &          px2,0)).lt.1d-12 .and. sgn(yy(i)-py1)*sgn(py2-yy(i))
     &           .ge.0) then
                    iflag(i)=1
             endif
          aa111(i)=dirdif(px1,xx(i),0)
          aa112(i)=dirdif(xx(i),px2,0)
          aa113(i)=dirdif(px1,px2,0)
          aa12(i)=sgn(yy(i)-py1)*sgn(py2-yy(i))
          write(18,*) aa111(i),aa112(i),aa113(i),aa12(i)             
             if(abs(dirdif(px1,xx(i),0)+dirdif(xx(i),px2,0)-dirdif(px1,
     &          px2,0)).lt.1d-12)   then
                iflag11(i)=1
             endif
             if(sgn(yy(i)-py1)*sgn(py2-yy(i)).ge.0)  then
                iflag12(i)=1
             endif
             write(18,*) iflag11(i),iflag12(i)
          endif
          if(py1.ge.py2)then
            if(abs(dirdif(px1,xx(i),0)+dirdif(xx(i),px2,0)-dirdif(px1,
     &           px2,0)).gt.1d-12 .or.
     &          sgn(yy(i)-py1)*sgn(py2-yy(i)).le.0)  then
                    iflag(i)=1
             endif
          endif             
       write(18,*) iflag(i)
       enddo
       return 
       end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  This subroutine dtermines whether the point is inside the polygon
c   Input:px, py, n, xx, yy, m
C
c   output: flag
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
       subroutine spolyse(px, py, n, xx, yy, m, flag)

       implicit real *8 (a-h, o-z)
      
       real*8 px(1:n+1), py(1:n+1)
       real*8 xx(1:m), yy(1:m), plon(1), plat(1)

       integer flag(1:m), npInfo, n, m
       integer pflag(1)

       call findInP(px, py, n, plon(1), plat(1))

       call splnpse(px,py, n, plon, plat, 1, 0, pflag)

       if(pflag(1).eq.-1) npinfo = 1
       if(pflag(1).eq.1) npinfo=0 

       write(*,*) xx(1), yy(1)            
  
       call splnpse(px,py, n, xx, yy, m, npinfo, flag)

       return 
       end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This subroutine determines whether the points (xx,yy) are
c     in a sphere polygon (px, py) or not
c     Input:
c      npinfo: 1, the north pole is inside; 0, the north pole is outside
c     output:
c      flag: -1, outside; 0, on the boundary; 1 inside
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
       subroutine splnpse(px, py, n, xx, yy, m, npinfo, flag)

       implicit real *8 (a-h, o-z)
      
       real*8 px(1:n+1), py(1:n+1), x, y, esp
       real*8 xx(1:m), yy(1:m)

       integer flag(1:m), npInfo, te, n, ip, sgn, is, m
       include 'param.inc'

       esp=1.0d-12

       px(n+1)=px(1)
       py(n+1)=py(1)

       do k=1, m

        x=xx(k)
        y=yy(k)

        ip = npinfo
        te = -1

       
 10     do i = 1, n

           call between(px(i), px(i+1), x, is, iecode)

           d1 = dirdif(x, px(i),1)
           d2 = dirdif(x,px(i+1),1)

           if(is.gt.0)then                  
                  call inter2Arc(px(i),py(i),px(i+1),py(i+1),0.0d0,
     &                pi*0.5d0, x, y, itag, iflag2,iflag3)

                  if(itag.eq.1) ip=ip+1
                  if(itag.eq.0.and.(iand(iflag2,8).eq.8 .or.
     &                             iand(iflag3,8).eq.8))then
                     te = 0
                     goto 910
                  endif
           endif
           if(is.eq.0)then
              if (x.eq.px(i).and.d2.gt.0.0d0)then
                  if(y-py(i).lt.0)then
                     ip =ip+1
                  elseif(y-py(i).eq.0)then
                     te=0
                     goto 910
                  endif
              elseif(x.eq.px(i).and.d2.lt.0.0d0)then
                  if(y-py(i).eq.0)then
                     te=0
                     goto 910
                  endif
              elseif(x.eq.px(i+1).and.d1.gt.0.0d0)then
                  if(y-py(i+1).lt.0) then
                      ip=ip+1
                  endif
                  if(y-py(i+1).eq.0)then
                     te=0
                     goto 910
                  endif
              elseif(x.eq.px(i+1).and.d1.lt.0.0d0)then
                  if(y-py(i).eq.0)then
                     te=0
                     goto 910
                  endif
              elseif(x.eq.px(i+1).and.px(i).eq.x)then
                  if(sgn(y-py(i))*sgn(py(i+1)-y).gt.0.0d0)then
                     te=0
                     goto 910
                  endif
              endif
           endif
            
 19     enddo

 60     if( ip - ip /2 * 2.eq.1) then
           te=1
          else
           te = -1
 69     endif

 910    continue
        flag(k)=te
       enddo

       return 
       end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This function returns the haversin of the arc distance between
C      (lon1, lat1) and (lon2, lat2) in rad. 
C       Input:
c            lon1,lat1, lon2, lat2: longitudes and latitudes in rad
C       Output:
C            havad: haversine of the distance.        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        real*8 function havad(lon1, lat1, lon2, lat2)

        implicit real*8 (a-h, o-z)
        real*8 lon1, lon2, lat1, lat2

        temp = cos(lat1)*cos(lat2)*hav(lon2-lon1)
        temp = temp + hav(lat2-lat1)
        havad=temp

        return
        end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      This function calculates the haversin of radians.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       real*8 function hav(rad)

       implicit real*8 (a-h,o-z)
      
       hav= sin(rad/2d0)**2
       return      
       end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This function calculates the arc distance of the haversin
C       Input:
c            val: haversine
C       Output:
C            ahav
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       real*8 function ahav(val)
       implicit real*8(a-h,o-z)

       if(val.lt.-1d-12.or.val.gt.1d0+1d-12)then
          write(*,*)'ahav() out of range', val
          stop
       endif
       if(val.lt.0d0)val=1d-12
       if(val.gt.1d0)val=1d0-1d-12
       
       temp=sqrt(val)
       ahav=2.0d0*asin(temp)
       return 
       end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       This function returns the arc distance between (lon1, lat1)
C       and (lon2, lat2) in rad. 
C       Input:
c            lon1,lat1, lon2, lat2: longitudes and latitudes in rad
C       Output:
C            arcdist: distance in rad.        
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        real*8 function arcdist(lon1, lat1, lon2, lat2)

        implicit real*8 (a-h, o-z)
        real*8 lon1, lon2, lat1, lat2

        temp = sqrt(havad(lon1, lat1, lon2, lat2))
        arcdist = 2.0*asin(temp)

        return
        end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Subroutine to get bearing starting from 
c      (xlon0, xlat0) to (xlon1, xlat1)
C      bx, by: driection 
C      bangle: bearing angle in radian
c      Input:
C        xlon0,xlat0,xlon1,xlat1: in radian
C      output: 
c        bx, by: direction vector
c        bangle: bearing angle in radian
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       subroutine beardirt(xlon0, xlat0, xlon1, xlat1, bx, by, bangle)

       implicit real*8 (a-h, o-z)
       include 'param.inc'

       dlon=xlon1-xlon0

       by = sin(dlon)*cos(xlat1)
       bx = cos(xlat0)*sin(xlat1)-sin(xlat0)*
     &      cos(xlat1)*cos(dlon)
       bangle= atan2(by,bx)
       bangle=mod((bangle+2.d0*pi),2.d0*pi)

       return
       end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  This subroutine calculates the area of a polygon
c  Input: xlon, xlat, np
C  Output: area
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
       subroutine spa(xlon, xlat, np, area)

       implicit real*8 (a-h, o-z)
       include 'param.inc'

       real*8 xlon(1:np), xlat(1:np)
       real*8 xlon1(10002), xlat1(10002)
 
       eps=1d-20
  
       k = 1
       if(abs(xlon(1)-xlon(np)).gt.eps.or.abs(xlat(1)-xlat(np))
     &     .gt.eps ) then
              xlon1(k)=xlon(1)
              xlat1(k)=xlat(1)
               k=k+1
       endif
       do i=2, np
       if(abs(xlon(i)-xlon(i-1)).gt.eps.or.abs(xlat(i)-xlat(i-1))
     &     .gt.eps ) then
              xlon1(k)=xlon(i)
              xlat1(k)=xlat(i)
              k=k+1
       endif
       enddo
       xlon1(k)=xlon1(1)
       xlat1(k)=xlat1(1)
       xlon1(k+1)=xlon1(2)
       xlat1(k+1)=xlat1(2)
       nt = k-1

       area = 0d0

       do i=2, nt+1

       call beardirt(xlon1(i), xlat1(i), xlon1(i-1), 
     &     xlat1(i-1), bx, by, bangle1)
       call beardirt(xlon1(i), xlat1(i), xlon1(i+1), 
     &     xlat1(i+1), bx, by, bangle2)
       
       temp = dirdif(bangle1,bangle2, 0)

       area = area + temp
       
       enddo

       area = area - (nt-2)*pi 

       return
       end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     This subroutine find an interior point inside a 
C     spherical polygon specified by xlon(), xlat()
C      Input:
C         xlon, xlat: array stored the vertice of the polygon
c                 must in counterclockwise
c         np: number of vertices of the spherical polygon
C      Output: 
c         plon, plat: coordinate of the interior point
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       subroutine findInP(xlon, xlat, np, plon, plat)

       implicit real*8 (a-h, o-z)
       real*8 xlon(1:np), xlat(1:np)
       real*8 xlon1(10002), xlat1(10002)

       eps=1e-20
  
       k = 1
       if(abs(xlon(1)-xlon(np)).gt.eps.or.abs(xlat(1)-xlat(np))
     &     .gt.eps ) then
              xlon1(k)=xlon(1)
              xlat1(k)=xlat(1)
              k=k+1
       endif

       do i=2, np
        if(abs(xlon(i)-xlon(i-1)).gt.eps.or.abs(xlat(i)-xlat(i-1))
     &     .gt.eps ) then
              xlon1(k)=xlon(i)
              xlat1(k)=xlat(i)
              k=k+1
        endif
       enddo
       xlon1(k)=xlon1(1)
       xlat1(k)=xlat1(1)
       nt = k-1

       call beardirt(xlon1(2), xlat1(2), xlon1(1), 
     &     xlat1(1), bx, by, bangle1)
       call beardirt(xlon1(2), xlat1(2), xlon1(3), 
     &     xlat1(3), bx, by, bangle2)

        temp = dirdif(bangle1, bangle2, 0)
        xbangle = dirpa(bangle1, temp/2.0d0, 0) 
        if(xbangle.gt.360d0)xbangle=xbangle-2.0*3.1415926535897d0

        call pwdistbear(xlon1(2), xlat1(2),xbangle, 0.3d0,plon,plat)

        itocheck=3
 
        do 10 while (itocheck.le.nt)
          call inter2Arc(xlon1(2), xlat1(2),plon,plat,  
     &       xlon1(itocheck), xlat1(itocheck), xlon1(itocheck+1),
     &       xlat1(itocheck+1),iff1,iff2,iff3)
             if(iff1.ge.0)then
               tmlon = plon
               tmlat = plat
               call fifpa(xlon1(2), xlat1(2),tmlon,tmlat,0.5d0,
     &                   plon,plat)
               goto 10
             endif
            itocheck=itocheck+1
 10     continue

       return
       end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  This subroutine finds a point along a spherical great arc segment that divides
c  the segment with a given fraction.
c
c   Input:
c       plon1, plat1, plon2, plat2: location of the great arc vertices, in degree.
c       frac: the given fraction, not limited between 0 and 1. 
c   Output:
c       xlon, xlat: coordination of the point found, in degree. 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       subroutine fifpa(plon1, plat1,plon2, plat2,frac, xlon, xlat)

       implicit real*8 (a-h, o-z)

       ad=arcdist(plon1, plat1, plon2, plat2)
       
       a = sin((1d0-frac)*ad)/sin(ad)
       b = sin(frac*ad)/sin(ad)

       x=a*cos(plat1)*cos(plon1)+b*cos(plat2)*cos(plon2)
       y=a*cos(plat1)*sin(plon1)+b*cos(plat2)*sin(plon2)
       z=a*sin(plat1)+b*sin(plat2)

       xlat = atan2(z, sqrt(x**2+y**2))
       xlon = atan2(y,x)

       return
       end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    find a point (plon, plat) from a point (xlon, xlat) with a
c      fixed bearing bangle (in degree) and arc distance arcdis (in degree)
c
c    Input:
c      xlon, xlat: the given reference point, in degree
c      bangle: bearing angle in degree, starting from the north direction anticlockwisely.
c      arcdis: distance in degree
c    Output:
c      plon, plat: coordination of target location in degree        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       subroutine pwdistbear(xlon, xlat, bangle, arcdis, plon, plat)
       implicit real*8 (a-h,o-z)
       include 'param.inc'

       plat = asin(sin(xlat)*cos(arcdis)+cos(xlat)
     &             *sin(arcdis)*cos(bangle))
       plon = xlon+atan2(sin(bangle)*sin(arcdis)*cos(xlat),
     &                cos(arcdis)-sin(xlat)*sin(plat))

       if (plon.gt.pi)  plon=plon-pi*2d0

       return
       end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    This subroutine calculates the Antipodes of a point.
C    Input: xc1, yc1
C    Output: xc2, yc2
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      subroutine antipl(xc1, yc1, xc2, yc2)
      implicit real*8 (a-h,o-z)
      include 'param.inc'

      yc2=-yc1
      
      if (xc1.ge.0d0) then
          xc2=xc1-pi
      else
          xc2=xc1+pi
      endif

      return
      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  This function converts degrees to radians.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      real*8 function deg2rad(deg)
      implicit real*8(a-h,o-z)
      include 'param.inc'
      
      dtr=pi/180d0
 
      deg2rad=deg*dtr

      return
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  This function calculates the difference between two direction.
c   Input:
c       a1, a2: the bearings [0, 2*Pi) or longitudes (-Pi, Pi]
c       ishort: 0, along the natural direction of bearings (clockwise)
c             or longitude (counterclockwise)
c   Output:
c       dirdif: difference between a1 and a2, ranging 0 to 2*Pi if
c              ishort=0 and -Pi to Pi if ishort=1, with a negative value
c              representing short direction is against the natural
c              direction
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      real*8 function dirdif(a1, a2, ishort)
      implicit real*8 (a-h,o-z)
      include 'param.inc'

      a=a2-a1

      if(ishort.eq.0) then
         do while(a.lt.0) 
            a=a+2.0d0*pi
         enddo
         do while(a.gt.2.0d0*pi) 
            a=a-2.0d0*pi
         enddo
      else
        do  while(a.le.-pi)
          a=a+2.0d0*pi
        enddo
        do  while(a.gt.pi)
          a=a-2.0d0*pi
        enddo
      endif
      
      dirdif=a
      return
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  This function calculates the new direction by adding an angle to
c  a known direction
c   Input:
c     a1: the bearing [0, 2*Pi) or longitude (-Pi, Pi]
c     da: the angle to be added, positive along natural direction 
c     itype: 0, bearing (clockwise)
c            1, longitude (counterclockwise)
c   Output:
c     dirpa: new directions by adding da to a1. It ranges between 0 and
c              2*Pi if itype=0 and -Pi to Pi if itype=1.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      real*8 function dirpa(a1, da, itype)
      implicit real*8 (a-h,o-z)
      include 'param.inc'
     
      a = mod (a1 + da,2*pi)

      if(itype.eq.0) then
         if(a.gt.2*pi) a=a-2d0*pi
      else
         if(a.le.-pi)a=a+2d0*pi
         if(a.gt.pi)a=a-2d0*pi
      endif
      
      dirpa=a
      return
      end



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c   This function calculates the arc length of angle bisector
c   
c   Input:
c       r1,r2,r12
c   Output:
c       bisect: distance in rad.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       real*8 function bisect(r1,r2,r12)

       implicit real*8 (a-h, o-z)

       temp1=(cos(r12)-cos(r2)*cos(r1))/sin(r2)/sin(r1)
       if(temp1.gt.1d0)temp1=1d0-1d-10
       if(temp1.lt.-1d0)temp1=-1d00+1d-10

       aa=acos(temp1)

       temp2 =(cos(r2)-cos(r1)*cos(r12))/sin(r1)/sin(r12)
       if(temp2.gt.1d0)temp2=1d0-1d-10
       if(temp2.lt.-1d0)temp2=-1d0+1d-10

       bb=acos(temp2)
       aa2=aa/2d0

       bisect=atan(sin(r1)/(cos(bb)/sin(bb)*sin(aa2)+cos(r1)*cos(aa2)))

       return
       end



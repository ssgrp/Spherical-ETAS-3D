
      integer::NNP,NND
      
      parameter(NNP=5000, NND=99000)

      integer::nprocs, myrank, ista, iend, ista2, iend2, ised
      common/mpi/nprocs, myrank, ista, iend, ista2, iend2, ised

      real(8):: tx(1:NNP),ty(1:NNP),tstart,tz, depmax
      common/domain/tx,ty,tstart,tz, depmax

      real(8)::xx(1:NND),yy(1:NND),zz(1:NND),zmg(1:NND),zdp(1:NND)
      common/data/xx,yy,zz,zmg,zdp

      real(8)::xa(1:NND), ya(1:NND)
      integer::indat(1:NND)
      common/antip/xa, ya,indat

      integer::ibndtyp,npoly,nn,ind(1:NND),nnc,mx,my,mz
      common/data2/ibndtyp,npoly,nn,ind,nnc,mx,my,mz

      real(8)::zprob(1:NND),zbandw(1:NND),zbkgd(1:NND),xint0
      common/back/zprob,zbandw,zbkgd,xint0

      real(8)::bwm,bwh
      integer::npp
      common/smoothing/bwm,bwh,npp
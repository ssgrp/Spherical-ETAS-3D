      parameter(nnp=5000, nnd=99000)

      common/mpi/nprocs, myrank, ista, iend, ista2, iend2, ised
      common/domain/tx(NNP),ty(NNP),tstart,tz, depmax
      common/data/xx(NND),yy(NND),zz(NND),zmg(NND),zdp(NND)
      common/antip/xa(NND), ya(NND),indat(NND)
      common/data2/ibndtyp,npoly,nn,ind(NND),nnc,mx,my,mz
      common/back/zprob(NND),zbandw(NND),zbkgd(NND),xint0
      common/smoothing/npp,bwm,bwh
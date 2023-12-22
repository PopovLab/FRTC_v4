module manager_mod
    !! модуль определяет начальные значения лучей и запускает трассировку
    use kind_module
    implicit none

contains

    subroutine manager(iterat,iw0, ntet, spectr)
        use constants            
        use plasma
        use rt_parameters, only : nr, ipri, iw, nmaxm, pabs0, eps, eps_const            
        use trajectory_module
        use spectrum_mod
        use iterator_mod,only: plost, pnab
        use dispersion_module, only: icall1, icall2, yn3, ivar, izn
        use driver_module !, only: irs, iabsorp
        use trajectory_data
        implicit none
        type (Spectrum) spectr
        type (SpectrumPoint) point
        real(wp) pabs
        !integer iznzap(mpnt),iwzap(mpnt),irszap(mpnt)
        !real(wp) rzap(mpnt),tetzap(mpnt),xmzap(mpnt),yn3zap(mpnt)
        integer ntet, iout, itr,  nnj,  n_it
        integer maxref, iterat, nmax0, ibad, itet, nref
        integer nbad1, nbad2, inz
        integer iw0, ifail, iabsirp, inak0,ib,ie
        integer nmax, i, nb1,nb2
        !integer iznzz, iwzz, irszz
        real(wp) htet, hr, yn, rin, xmin, rstart
        real(wp) powexit, dltpow,  pow1, pgamma, xm
        real(wp) tetin0, tetin, tet

        pabs = spectr%max_power*pabs0/1.d2
        print *, 'pabs =', pabs, spectr%max_power, pabs0
        !lenstor = length
        htet = zero
        hr = 1.d0/dble(nr+1) !sav2008
        if (ntet.ne.1) htet = (tet2-tet1)/(ntet-1)
        irs = 1
        iout = 0
        !mbeg(1) = 1
        itr = 0
        !inak = 0
        nnj = 0
        do n_it = 0,3
            nnj = nnj+nmaxm(n_it+1)
        end do
        maxref = nnj
        if (iterat.lt.3) nmax0=nmaxm(iterat+1)
        if (iterat.ge.3) nmax0=nmaxm(4)
        if (ipri.gt.1) then
            write(*,1001) iterat+1
            write(*,1002)
        end if
        ibad = 0
        !--------------------------------------
        ! begin outer loop on teta
        !--------------------------------------
        do itet = 1,ntet
            nref = 0
            nbad1 = 0
            nbad2 = 0
            icall1 = 0
            icall2 = 0
            tetin = tet1+htet*(itet-1)
            !--------------------------------------
            ! begin inner loop on nz
            !--------------------------------------
            do inz = 1, spectr%size
                itr = itr+1
                current_trajectory => trajectories(itr)
                !ipri          if(ipri.eq.4)  write(23,*)
                point = spectr%data(inz)
                if(iterat.eq.0) then
                    !-----------------------------------------
                    !    find initial radius for a trajectory
                    !    on the 1th iteration
                    !-----------------------------------------
                    call current_trajectory%init(tetin, inz)

                    yn = point%Ntor
                    pow = point%power
                    !yn=ynzm(inz) !sav2008, yn is introduced
                    !pow=pm(inz)
                    irs = 1
                    iw = iw0
                    rin = rini(xmin,tetin,point,hr,ifail)
                    current_trajectory%rin = rin
                    if (ifail.eq.1) then
                        if (ipri.gt.1) write (*,*) 'error: no roots'
                        iabsorp = -1
                        !inak0 = inak
                        go to 10
                    end if

                else
                    if (current_trajectory%mbad.ne.0) then
                        plost = plost+point%power
                        go to 31
                    end if
                    !ib = mbeg(itr)
                    !ie = mend(itr)
                    powexit = point%power
                    dltpow = pabs
                    call dqliter(dltpow,current_trajectory,hr,powexit,iout)
                    if (nmax0.eq.0) then
                        !ib = mbeg(itr)
                        !ie = mend(itr)
                        pow1 = powexit
                        pgamma = 1.d0-pow1/point%power
                        powexit = pow1/pgamma
                        dltpow = powexit-pow1+pabs
                        call dqliter(dltpow,current_trajectory,hr,powexit,iout)
                        powexit = powexit-dltpow+pabs
                        if (powexit.lt.zero) powexit=zero
                        go to 30
                    end if
	                if (iout.eq.0) then
                        go to 30
                    else
                        tetin = current_trajectory%tetzap
                        xmin = current_trajectory%xmzap
                        rin = current_trajectory%rzap
                        yn3 = current_trajectory%yn3zap
                        pow = current_trajectory%powexit
                        irs = current_trajectory%irszap
                        iw =  current_trajectory%iwzap
                        izn = current_trajectory%iznzap
                        !jrad(ie+1) = 1
                        !dland(ie+1) = lfree
                        !inak = lfree-1
                        !call current_trajectory%reset(0)
                        ! продолжение траектории !! ресет не нужен 
                    end if
                end if
                !---------------------------------------
                ! initial parameters for a trajectory
                !---------------------------------------
                xm = xmin
                rstart = rin !sav2008
                tet = tetin
                nmax = nmax0
                iabsorp = 0
                !inak0 = inak
                !-------------------------------------
                ! call ray tracing
                !-------------------------------------
                call traj(xm,tet,rstart,nmax,nb1,nb2,itet,inz, pabs) !sav2009
                eps = eps_const 
                nbad1 = nbad1+nb1
                nbad2 = nbad2+nb2
                current_trajectory%nrefj = current_trajectory%nrefj + nmax
                powexit = pow
                nref = nref+nmax
10              if (iabsorp.lt.0) then
                    !-------------------------------------
                    !    encounted problems
                    !-------------------------------------
                    !if (inak.eq.lenstor-1) then
                    if (current_trajectory%size.eq.max_size-1) then
                        write (*,*) 'fix maximal length'
                        nmax0 = 0
                        do i=1,4
                            nmaxm(i) = 0
                        end do
                        iout = 1
                        goto 20
                    end if
                    if (ipri.gt.1) then
                        tetin0=tet1+htet*(itet-1)
                        write (*,111) tetin0, point%Ntor

111                     format(1x,'traj. with tet0=',f10.5,1x,', Ninput=',f10.5,1x,'failed')
                    end if
                    current_trajectory%mbad = 1 ! плохоая траектория
                    plost= plost+pow
                    !inak = inak0
                    !mend(itr) = inak-1
                    goto 30
                end if
                !---------------------------------------
                ! remember end point of trajectory
                !---------------------------------------
                current_trajectory%rzap   = rzz
                current_trajectory%tetzap = tetzz
                current_trajectory%xmzap  = xmzz
                current_trajectory%yn3zap = yn3
                current_trajectory%iznzap = iznzz
                current_trajectory%iwzap  = iwzz
                current_trajectory%irszap = irszz
                if (iterat.eq.0) then
                    !if (itr.gt.1) mbeg(itr) = mend(itr-1)+2
                    !mend(itr) = inak
                    !jrad(mend(itr)+1) = 0
                    !lfree = mend(itr)+2
                    !inak = lfree-1
                end if
20              continue
                if(iout.ne.0) then
                    print *,'**************'
                    !dcoll(ie+1) = inak
                    !jrad(inak+1) = 0
                    !lfree = inak+2
                end if
                if(current_trajectory%nrefj.gt.maxref.and.pow.gt.pabs) then !forced absorp
                    if(pow.ge.point%power) go to 30 !sav2008
                    !ib = mbeg(itr)
                    !ie = mend(itr)
                    pow1 = pow
                    pgamma = 1.d0-pow1/point%power
                    powexit = pow1/pgamma
                    dltpow = powexit-pow1+pabs
                    call dqliter(dltpow, current_trajectory, hr,powexit,iout)
                    powexit = powexit-dltpow+pabs
                    if(powexit.lt.zero) powexit=zero
                end if
30              continue
                pnab = pnab+powexit
31              continue
            end do
            if(ipri.gt.1) write(*,1003)itet,icall1,icall2,nref,nbad1,nbad2
        end do
1001    format (30x,i4,' iteration')
1002    format (6x,'n',5x,'call2',6x,'call4',6x,'nrefl',4x,'last',5x,'bad2',5x,'bad4')
1003    format (3x,i4,2(1x,i10),2x,i7,2x,i8,2(1x,i7),2(2x,i7))
1004    format(1x,i8)
1005    format(1x,i5)
1006    format (e14.7)
    end    


    real(wp) function rini(xm, tet, point, hr, ifail) !sav2009
        use constants, only : zero
        use rt_parameters, only : inew
        use spectrum_mod, only : SpectrumPoint
        use dispersion_module, only: ivar, yn3
        use dispersion_module, only: disp2_iroot2
        use metrics, only: g22, g33, co, si
        use metrics, only: calculate_metrics
        implicit none

        type(SpectrumPoint), intent(in) :: point
        real(wp), intent(inout)          :: xm
        real(wp), intent(in)             :: tet,  hr
        integer, intent(inout)           :: ifail

        integer :: ntry
        real(wp) :: pa, prt, prm
        real(wp) :: f1,f2

        real(wp),  parameter :: rhostart=1.d0
        integer,   parameter :: ntry_max=5

        ifail = 1
        rini = zero
        ntry = 0
        pa = rhostart
        do while (ntry.lt.ntry_max.and.pa.ge.2d0*hr)
            pa = rhostart-hr*dble(ntry)-1.d-4
            ntry = ntry+1

            ! вычисление g22 и g33
            call calculate_metrics(pa, tet)

            yn3 = point%Ntor*dsqrt(g33)/co 
            xm = point%Npol*dsqrt(g22)/si

            call disp2_iroot2(pa,xm,tet,f1,f2)
            
            if (f1.ge.zero.and.f2.ge.zero) then
                rini = pa
                ifail = 0
                return
            end if
        end do
    end      

    subroutine dqliter(dltpow, traj, h, powexit, iout) !sav2008
        use constants, only: clt, zero
        use rt_parameters, only: itend0, kv
        use iterator_mod, only: vlf, vrt, dflf, dfrt
        use iterator_mod, only: distr
        use decrements, only: pdec1, pdec2, pdec3, pdecv
        use decrements, only: zatukh
        use current, only: dfind
        use plasma, only: vperp
        use iterator_mod, only: psum4
        !use driver_module, only: jrad, iww, length      
        !use driver_module, only: pow, vel, perpn, dland, dcoll, dalf
        use driver_module, only: pow
        use trajectory_data
        implicit none

        type(Trajectory), pointer, intent(in) :: traj
        real(wp), intent(in)   :: dltpow
        real(wp), intent(in)   :: h
        real(wp), intent(out)  :: powexit
        !integer, intent(inout) :: ib, ie
        integer, intent(inout) :: iout

        type(TrajectoryPoint) :: tp
        integer  :: i, iv,  jr, ifast, jchek
        real(wp) :: pdec1z, pdec3z, pintld, pintal
        real(wp) :: v, refr, dek3, argum, valfa
        real(wp) :: df, dfsr, vsr, pcurr, dcv
        real(wp) :: powpr, powd, powcol, powal
        real(wp) :: pil, pic, pia

        pow=powexit
        pdec1=zero
        pdec1z=zero
        pdec3=zero
        pdec3z=zero
        pdecv=zero
        pintld=zero
        pintal=zero
  10    continue
        iout=0
        do i = 1, traj%size
            !-----------------------------------
            ! restore memorized decrements and
            ! integrate power equation
            !------------------------------------
            tp = traj%points(i)
            v=tp%vel
            jr=tp%jrad
            refr=tp%perpn
            ifast=tp%iww
            dek3=zero
            if(itend0.gt.0) then
                argum=clt/(refr*valfa)
                dek3=zatukh(argum,abs(jr),vperp,kv)
            end if
            !!!old variant
            !!!       call raspr(v,abs(jr),iv,df)
            !!!       if(iv.eq.0) iv=1
            !!!!!!!!!!!!!!!!!!!!!!!!!!
            call distr(v,abs(jr),iv,df)
            !!       dfsr=v*df*(vrt-vlf)
            !!       vsr=v*(vrt-vlf)
            dfsr=(vlf*dflf+vrt*dfrt)/2d0*(vrt-vlf) !sav2008
            vsr=(vrt+vlf)*(vrt-vlf)/2d0 !sav2008
            if(jr.lt.0) then !case of turn
                jr=-jr
                !variant        pintld=-dland(i)*df
                !!        pintld=-dland(i)*(dflf+dfrt)/2d0
                pintld=dabs(tp%dland*(dflf+dfrt)/2d0)
                pdec2=dexp(-2d0*tp%dcoll)
                pintal=dabs(tp%dalf*dek3)
                pcurr=pdec2*dexp(-2d0*pintld-2d0*pintal)
                psum4=psum4+pow*(1d0-pcurr)
                dcv=tp%dland/vsr
            else
                pdec2=tp%dcoll
                pdecv=tp%dland
                !!        pdec1=-pdecv*df
                pdec1=dabs(pdecv*df)
                pdec3=dabs(tp%dalf*dek3)
                pintld=(pdec1+pdec1z)/2d0*h
                pintal=(pdec3+pdec3z)/2d0*h
                pdec1z=pdec1
                pdec3z=pdec3
                dcv=pdecv*h/vsr
            end if
            powpr=pow
            if(dltpow.ne.zero) then
                powd=pow*dexp(-2d0*pintld)
                powcol=powd*pdec2
                powal=powcol*dexp(-2d0*pintal)
                pow=powal
            end if
            pil=pintld
            pic=.5d0*dabs(dlog(pdec2))
            pia=pintal
            call dfind(jr,iv,v,powpr,pil,pic,pia,dfsr,dcv, &
                    refr,vlf,vrt,ifast)
            if(pow.lt.dltpow) then
                powexit=pow
                return
            end if
        end do
        jchek= 0 !jrad(ie+1) !!!!!
        !c-------------------------------------------
        !c  check whether trajectory has continuation
        !c---------------------------------------------
        if(jchek.eq.0) then
            iout=1
            powexit=pow
            return
        else
            print *,'dqliter stop'
            stop
            !ib=idnint(dland(ie+1))
            !ie=idnint(dcoll(ie+1))
            goto 10
        end if
    end    

end module manager_mod

      subroutine gibbs(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10
		
	real*8 yy(ni,nd,nk), tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), vv(nl)
      real*8 aalpha(nl), sigmaa
      real*8 lambda(ni), sigmal, api
     
      real*8 sigma_zmu, sigma_xbeta
      real*8 sigma_wbeta, sigma_ubeta
      real*8 a1, b1, a2, b2, a3, b3
      real*8 a4, b4, a5, b5, a6, b6
                
      common /vyy/yy
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vvv/vv
      
      common /vaalpha/aalpha
      common /vsigmaa/sigmaa
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

      common /vsigma_zmu/sigma_zmu
      common /vsigma_xbeta/sigma_xbeta
      common /vsigma_wbeta/sigma_wbeta
      common /vsigma_ubeta/sigma_ubeta
      common /va1/a1
      common /vb1/b1
      common /va2/a2
      common /vb2/b2
      common /va3/a3
      common /vb3/b3
      common /va4/a4
      common /vb4/b4
      common /va5/a5
      common /vb5/b5
      common /va6/a6
      common /vb6/b6

      call gibbs_bij(iseed)          

      call gibbs_taui(iseed)    
     
      call gibbs_rij(iseed)          
      
      
      call gibbs_lambda(iseed)          

      call gibbs_zeta(iseed)     
      

      call gibbs_zmu(iseed)
      
      call gibbs_xbeta(iseed)     


      call gibbs_ubeta(iseed)  
      
      call gibbs_wbeta(iseed)     

      call gibbs_rho(iseed)         

      
      call gibbs_sigmat(iseed)

      call gibbs_sigmar(iseed)

      call gibbs_sigmaz(iseed)     

      call gibbs_sigmal(iseed)          
      
      
      call gibbs_vv(iseed)          
      
	call gibbs_sigmaa(iseed)

      call gibbs_aalpha(iseed)                
      
      end subroutine

                 
	subroutine gibbs_aalpha(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl)
      real*8 aalpha(nl)
      real*8 lambda(ni), sigmal, api

      real*8 a6, b6
      
	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio

      integer ldum

      real*8 zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1
      real*8 bik, bik1, bik2, dd, dd1, dd2
      real*8 temp, rstar, bstar, bstar1, bstar2, bsigma
      real*8 summ, suma, sumb, ttemp
                  
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

      common /va6/a6
      common /vb6/b6
      
	common /vldum/ldum

      external falpha, drnnof, drnunf
      
      do l = 1, nl
          
          ldum = l

          bold = aalpha(l)
      
          nopt = 1
	    reqmin = 1.0d-10 ; konvge = 5 ; kcount = 1000
	    step(1) = 0.2d0 ; start(1) = bold
          call nelmin(falpha, nopt, start, xmin, ynewlo,  
     +                reqmin, step, konvge, kcount, 
     +                icount, numres, ifault)
	    amean = xmin(1)

          aalpha(l) = amean

          summ = 0.d0      
          do i = 1, ni         
          
              zb = 0.d0
              do m = 1, nz
                  zb = zb + zz(i,m)*zmu(m)
              enddo

              wb = 0.d0
              do m = 1, nw
                  wb = wb + ww(i,m)*wbeta(m)
              enddo      
              sigmai = dexp(wb)

              xb = sigmaz*zeta(i)
              do m = 1, nx
                  xb = xb + xx(i,m)*xbeta(m)
              enddo      
              pp = dexp(xb)/(1.d0 + dexp(xb))           
              pp1 = 1.d0 - pp

              ub = sigmal*lambda(i)
              do m = 1, nu
                  ub = ub + uu(i,m)*ubeta(m)
              enddo     

              do k = 1, nk

                  bik = 0.d0 
                  do ll = 1, nl

                      eub = dexp(aalpha(ll) + ub)
                      dd = eub/(1.d0 + eub)                     
                  
                      temp = 2.d0*api*dfloat(ll)*(tik(i,k) + dd)
                      bik = bik + dsqrt(2.d0*delta(ll))*dcos(temp)
                  
                  enddo

                  eub = dexp(aalpha(l) + ub)
                  dd = eub/(1.d0 + eub)                               
                  dd1 = eub/(1.d0 + eub)**2
                  dd2 = eub*(1.d0 - eub)/(1.d0 + eub)**3
              
                  temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                  
                  bik1 = -dsqrt(2.d0*delta(l))*dsin(temp)
     +                    *2.d0*api*dfloat(l)*dd1
                  
                  bik2 = -dsqrt(2.d0*delta(l))
     +                 *(dcos(temp)*(2.d0*api*dfloat(l)*dd1)**2 
     +                   + dsin(temp)*2.d0*api*dfloat(l)*dd2)
                          
                  do j = 1, nd        
            
                      jj = (j - 1)*nk + k      

                      if (jj .eq. 1) then

                          rstar = bij(i,jj) 
     +                          - sigmar*rij(i,jj)
     +                          - zb - sigmat*taui(i)

                          bstar = sigmai*pp*bik
                          bstar1 = sigmai*pp*bik1
                          bstar2 = sigmai*pp*bik2
                      
                          bsigma = (sigmai*pp1)**2
                      
                      else

                          rstar = bij(i,jj) - rho*bij(i,jj-1) 
     +                          - sigmar*(rij(i,jj) 
     +                                    - rho*rij(i,jj-1)) 
     +                          - (1.d0 - rho)
     +                          * (zb + sigmat*taui(i))

                          bstar = (1.d0 - rho)*sigmai*pp*bik
                          bstar1 = (1.d0 - rho)*sigmai*pp*bik1
                          bstar2 = (1.d0 - rho)*sigmai*pp*bik2
                      
                          bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                                        
                      endif
                      
                      summ = summ 
     +                     + (-bstar1**2 + (rstar - bstar)*bstar2)
     +                     / bsigma

                  enddo
            
              enddo

          enddo

          suma = a6; sumb = b6
          do ll = 1, nl
          
              suma = suma + 1.d0/2.d0
              sumb = sumb + aalpha(ll)**2/(2.d0*delta(ll))
          
          enddo
         
          ttemp = aalpha(l)/delta(l)
          
          asigma = summ - suma*(sumb/delta(l) - ttemp**2)/sumb**2
            
	    asigma = -1.d0/asigma*2.d0
          if (asigma .lt. 0.0d0) asigma = -asigma      

          bpdf = -falpha(bold) 
     +         + (bold - amean)**2/(2.d0*asigma)

          do ii = 1, 50
	  
              call rnset(iseed)
		    rv = drnnof()
		    call rnget(iseed)
		    anew = amean + rv*dsqrt(asigma)

		    apdf = -falpha(anew) 
     +             + (anew - amean)**2/(2.d0*asigma)
              
		    ratio = apdf - bpdf 
     
              if (ratio .ge. 0.0d0) then
		        bpdf = apdf
			    bold = anew
	        else
		        call rnset(iseed)
			    u = drnunf()
			    call rnget(iseed)
			    if (dlog(u) .le. ratio) then
			        bpdf = apdf
				    bold = anew
			    endif
		    endif
      		
          enddo
            
          aalpha(l) = bold
              
      enddo

      end subroutine


      real*8 function falpha(astar)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl)
      real*8 aalpha(nl)
      real*8 lambda(ni), sigmal, api

      real*8 a6, b6

      integer ldum
      real*8 astar, zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1, dd, bik
      real*8 temp, rstar, bstar, bsigma
      real*8 summ, suma, sumb, pdf
                  
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

      common /va6/a6
      common /vb6/b6
      
	common /vldum/ldum

      l = ldum
      
      aalpha(l) = astar

      summ = 0.d0      
      do i = 1, ni         
      
          zb = 0.d0
          do m = 1, nz
              zb = zb + zz(i,m)*zmu(m)
          enddo

          wb = 0.d0
          do m = 1, nw
              wb = wb + ww(i,m)*wbeta(m)
          enddo      
          sigmai = dexp(wb)
      
          xb = sigmaz*zeta(i)
          do m = 1, nx
              xb = xb + xx(i,m)*xbeta(m)
          enddo      
          pp = dexp(xb)/(1.d0 + dexp(xb))           
          pp1 = 1.d0 - pp

          ub = sigmal*lambda(i)
          do m = 1, nu
              ub = ub + uu(i,m)*ubeta(m)
          enddo      

          do k = 1, nk

              bik = 0.d0
              do ll = 1, nl
                  eub = dexp(aalpha(ll) + ub)
                  dd = eub/(1.d0 + eub)                      
                  temp = 2.d0*api*dfloat(ll)*(tik(i,k) + dd)
                  bik = bik + dsqrt(2.d0*delta(ll))*dcos(temp)
              enddo

              bstar = sigmai*pp*bik
      
              do j = 1, nd        
            
                  jj = (j - 1)*nk + k      

                  if (jj .eq. 1) then

                      rstar = bij(i,jj) 
     +                      - sigmar*rij(i,jj)  
     +                      - zb - sigmat*taui(i)

                      bstar = sigmai*pp*bik

                      bsigma = (sigmai*pp1)**2
                      
                  else

                      rstar = bij(i,jj) - rho*bij(i,jj-1) 
     +                      - sigmar*(rij(i,jj) - rho*rij(i,jj-1)) 
     +                      - (1.d0 - rho)
     +                      * (zb + sigmat*taui(i))

                      bstar = (1.d0 - rho)*sigmai*pp*bik

                      bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                                        
                  endif

                  summ = summ - (rstar - bstar)**2/(2.d0*bsigma)
                  
              enddo
              
          enddo           

      enddo           
      
      suma = a6; sumb = b6
      do ll = 1, nl
          
          suma = suma + 1.d0/2.d0
          sumb = sumb + aalpha(ll)**2/(2.d0*delta(ll))
          
      enddo
      
      pdf = summ - suma*dlog(sumb)

	falpha = -pdf         	    
     	    
      end function  
                                                           
      
	subroutine gibbs_sigmaa(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

      real*8 aalpha(nl), delta(nl), sigmaa
      real*8 a6, b6

      real*8 shape, scale, rv
      
      common /vaalpha/aalpha
      common /vdelta/delta
      common /vsigmaa/sigmaa

      common /va6/a6
      common /vb6/b6
      
      external drngam
      
      shape = a6; scale = b6
      do l = 1, nl
          
          shape = shape + 1.d0/2.d0
          scale = scale + aalpha(l)**2/(2.d0*delta(l))
          
      enddo

	call rnset(iseed)
	call drngam(1, shape, rv)
	call rnget(iseed)
	sigmaa = scale/rv
      
      end subroutine
      
      
	subroutine gibbs_vv(iseed) 
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10
	integer, parameter :: npoint = 31, ndim = 3

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), vv(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api
      real*8 a5, b5, a6, b6

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio

      real*8 epsilon, astar(npoint), diffa(npoint)
      real*8 yystar(npoint), xxstar(npoint,ndim)
      real*8 xxx(ndim,ndim), xy(ndim), amaxy
      real*8 xhat, detx, xinvy
      
      integer ldum
      real*8 sumv, sumd
                  
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vvv/vv
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

      common /va5/a5
      common /vb5/b5
      common /va6/a6
      common /vb6/b6
      
	common /vldum/ldum

      external fvv, drnnof, drnunf

      epsilon = 0.1d0

      diffa(1) = -1.5d0
      do ii = 2, npoint
          diffa(ii) = diffa(ii-1) + 0.1d0
      enddo
      
      do l = 1, nl-1

          ldum = l

          bold = dlog(vv(l)/(1.d0 - vv(l)))
      
          nopt = 1
	    reqmin = 1.0d-10 ; konvge = 5 ; kcount = 1000
	    step(1) = 0.2d0 ; start(1) = bold
          call nelmin(fvv, nopt, start, xmin, ynewlo,  
     +                reqmin, step, konvge, kcount, 
     +                icount, numres, ifault)
	    amean = xmin(1)
          
          do kk = 1, npoint
              astar(kk) = amean + diffa(kk)*epsilon
          enddo

          do kk = 1, npoint
              yystar(kk) = -fvv(astar(kk))
              xxstar(kk,1) = astar(kk)**2
              xxstar(kk,2) = astar(kk)
              xxstar(kk,3) = 1.d0
          enddo

          amaxy = dabs(yystar(1))
          do kk = 2, npoint
              if (dabs(yystar(kk)) .gt. amaxy) then
                  amaxy = dabs(yystar(kk))
              endif
          enddo

          do kk = 1, npoint
              yystar(kk) = yystar(kk)/amaxy
          enddo
        
          do k1 = 1, ndim
              do k2 = 1, ndim
                  xxx(k1,k2) = 0.d0
                  do kk = 1, npoint
                      xxx(k1,k2) = xxx(k1,k2) 
     +                           + xxstar(kk,k1)
     +                             *xxstar(kk,k2)
                  enddo
              enddo
              xy(k1) = 0.d0
              do k2 = 1, npoint
                  xy(k1) = xy(k1) + xxstar(k2,k1)*yystar(k2)
              enddo
          enddo
        
          detx = xxx(1,1)*xxx(2,2)*xxx(3,3)
     +         + xxx(2,1)*xxx(3,2)*xxx(1,3)
     +         + xxx(3,1)*xxx(1,2)*xxx(2,3)
     +         - xxx(1,1)*xxx(3,2)*xxx(2,3)
     +         - xxx(3,1)*xxx(2,2)*xxx(1,3)
     +         - xxx(2,1)*xxx(1,2)*xxx(3,3)   
                                              
          xinvy = (xxx(2,2)*xxx(3,3) - xxx(2,3)*xxx(3,2))*xy(1) 
     +          + (xxx(1,3)*xxx(3,2) - xxx(1,2)*xxx(3,3))*xy(2)
     +          + (xxx(1,2)*xxx(2,3) - xxx(1,3)*xxx(2,2))*xy(3)     
          
          xhat = xinvy/detx
              
          asigma = -1.d0/(2.d0*amaxy*xhat)*2.d0
          if (asigma .lt. 0.0d0) asigma = -asigma      

          bpdf = -fvv(bold) 
     +         + (bold - amean)**2/(2.d0*asigma)

          do ii = 1, 50
	  
              call rnset(iseed)
	        rv = drnnof()
	        call rnget(iseed)
		    anew = amean + rv*dsqrt(asigma)

		    apdf = -fvv(anew) 
     +             + (anew - amean)**2/(2.d0*asigma)
              
		    ratio = apdf - bpdf 
     
              if (ratio .ge. 0.0d0) then
		        bpdf = apdf
		        bold = anew
	        else
                  call rnset(iseed)
		        u = drnunf()
		        call rnget(iseed)
		        if (dlog(u) .le. ratio) then
		            bpdf = apdf
		            bold = anew
		        endif
		    endif
		
          enddo     
                    
          vv(l) = dexp(bold)/(1.d0 + dexp(bold))
        
      enddo

      sumd = 0.d0
      do l1 = 1, nl
                   
          sumv = dlog(vv(l1))                
          do l2 = 1, l1-1
              sumv = sumv + dlog(1.d0 - vv(l2))
          enddo
          
          delta(l1) = dexp(sumv)

          sumd = sumd + delta(l1)
          
      enddo      
      do ll = 1, nl
          delta(ll) = delta(ll)/sumd
      enddo
      
      end subroutine      

      
      real*8 function fvv(astar)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), vv(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api
      real*8 a5, b5, a6, b6

      integer ldum
      real*8 astar, zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1, dd, bik, sumv, sumd
      real*8 temp, rstar, bstar, bsigma
      real*8 summ, suma, sumb, sumc, pdf
                  
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vvv/vv
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

      common /va5/a5
      common /vb5/b5
      common /va6/a6
      common /vb6/b6
      
	common /vldum/ldum

      l = ldum

      vv(l) = dexp(astar)/(1.d0 + dexp(astar))

      sumd = 0.d0
      do l1 = 1, nl
                   
          sumv = dlog(vv(l1))                
          do l2 = 1, l1-1
              sumv = sumv + dlog(1.d0 - vv(l2))
          enddo
          
          delta(l1) = dexp(sumv)

          sumd = sumd + delta(l1)
          
      enddo      
      do ll = 1, nl
          delta(ll) = delta(ll)/sumd
      enddo
             
      summ = 0.d0      
      do i = 1, ni         
      
          zb = 0.d0
          do m = 1, nz
              zb = zb + zz(i,m)*zmu(m)
          enddo

          wb = 0.d0
          do m = 1, nw
              wb = wb + ww(i,m)*wbeta(m)
          enddo      
          sigmai = dexp(wb)
      
          xb = sigmaz*zeta(i)
          do m = 1, nx
              xb = xb + xx(i,m)*xbeta(m)
          enddo      
          pp = dexp(xb)/(1.d0 + dexp(xb))           
          pp1 = 1.d0 - pp

          ub = sigmal*lambda(i)
          do m = 1, nu
              ub = ub + uu(i,m)*ubeta(m)
          enddo      

          do k = 1, nk

              bik = 0.d0
              do ll = 1, nl
                  eub = dexp(aalpha(ll) + ub)
                  dd = eub/(1.d0 + eub)                      
                  temp = 2.d0*api*dfloat(ll)*(tik(i,k) + dd)
                  bik = bik + dsqrt(2.d0*delta(ll))*dcos(temp)
              enddo
      
              do j = 1, nd        
            
                  jj = (j - 1)*nk + k      

                  if (jj .eq. 1) then

                      rstar = bij(i,jj) 
     +                      - sigmar*rij(i,jj)
     +                      - zb - sigmat*taui(i)

                      bstar = sigmai*pp*bik

                      bsigma = (sigmai*pp1)**2
                      
                  else

                      rstar = bij(i,jj) - rho*bij(i,jj-1) 
     +                      - sigmar*(rij(i,jj) - rho*rij(i,jj-1)) 
     +                      - (1.d0 - rho)
     +                      * (zb + sigmat*taui(i))

                      bstar = (1.d0 - rho)*sigmai*pp*bik

                      bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                                        
                  endif
                            
                  summ = summ - (rstar - bstar)**2/(2.d0*bsigma)
                                    
              enddo
              
          enddo           

      enddo           
      
      suma = a6; sumb = b6
      do ll = 1, nl
          
          suma = suma + 1.d0/2.d0
          sumb = sumb + aalpha(ll)**2/(2.d0*delta(ll))
          
      enddo
      
      sumc = 0.d0
      do ll = 1, nl
          
          sumc = sumc + dlog(delta(ll))/2.d0
          
      enddo
          
      pdf = summ 
     +    + (a5 - 1.d0)*dlog(vv(l)) 
     +    + (b5 - 1.d0)*dlog(1.d0 - vv(l))
     +    - suma*dlog(sumb) - sumc
     +    + astar - 2.d0*dlog(1.d0 + dexp(astar))

	fvv = -pdf         	    
     	    
      end function  
     
  
	subroutine gibbs_sigmal(iseed) 
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api
      real*8 a4, b4

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio

      real*8 zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1
      real*8 bik, bik1, bik2, dd, dd1, dd2
      real*8 temp1, temp2     
      real*8 rstar, bstar, bstar1, bstar2, bsigma
      real*8 star, summ
                        
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

      common /va4/a4
      common /vb4/b4

      external fsigmal, drnnof, drnunf      
     
      bold = dlog(sigmal)
      
      nopt = 1
	reqmin = 1.0d-10 ; konvge = 5 ; kcount = 1000
	step(1) = 0.2d0 ; start(1) = bold
      call nelmin(fsigmal, nopt, start, xmin, 
     +            ynewlo, reqmin, step, konvge, 
     +            kcount, icount, numres, ifault)
	amean = xmin(1)      

      sigmal = dexp(amean)

      summ = 0.d0      
      do i = 1, ni         
            
          zb = 0.d0
          do m = 1, nz
              zb = zb + zz(i,m)*zmu(m)
          enddo

          wb = 0.d0
          do m = 1, nw
              wb = wb + ww(i,m)*wbeta(m)
          enddo      
          sigmai = dexp(wb)

          xb = sigmaz*zeta(i)
          do m = 1, nx
              xb = xb + xx(i,m)*xbeta(m)
          enddo      
          pp = dexp(xb)/(1.d0 + dexp(xb))           
          pp1 = 1.d0 - pp

          ub = sigmal*lambda(i)
          do m = 1, nu
            ub = ub + uu(i,m)*ubeta(m)
          enddo    
          
          temp1 = sigmal*lambda(i)

          do k = 1, nk

              bik = 0.d0; bik1 = 0.d0; bik2 = 0.d0 
              do l = 1, nl
                  
                  eub = dexp(aalpha(l) + ub)
                  
                  dd = eub/(1.d0 + eub)                             
                  dd1 = temp1*eub/(1.d0 + eub)**2
                  dd2 = temp1*eub
     +                * ((1.d0 + temp1)*(1.d0 + eub) 
     +                    - 2.d0*temp1*eub)
     +                / (1.d0 + eub)**3
                  
                  temp2 = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                  
                  bik = bik + dsqrt(2.d0*delta(l))*dcos(temp2)
                  
                  bik1 = bik1 
     +                 - dsqrt(2.d0*delta(l))*dsin(temp2)
     +                   *2.d0*api*dfloat(l)*dd1
                  
                  bik2 = bik2 
     +                 - dsqrt(2.d0*delta(l))
     +                   *(dcos(temp2)
     +                     *(2.d0*api*dfloat(l)*dd1)**2 
     +                     + dsin(temp2)*2.d0*api*dfloat(l)*dd2)
              enddo

              do j = 1, nd        
            
                  jj = (j - 1)*nk + k      

                  if (jj .eq. 1) then

                      rstar = bij(i,jj) 
     +                      - sigmar*rij(i,jj)
     +                      - zb - sigmat*taui(i)

                      bstar = sigmai*pp*bik
                      bstar1 = sigmai*pp*bik1
                      bstar2 = sigmai*pp*bik2
                      
                      bsigma = (sigmai*pp1)**2
                      
                  else

                      rstar = bij(i,jj) - rho*bij(i,jj-1) 
     +                      - sigmar*(rij(i,jj) - rho*rij(i,jj-1)) 
     +                      - (1.d0 - rho)
     +                      * (zb + sigmat*taui(i))

                      bstar = (1.d0 - rho)*sigmai*pp*bik
                      bstar1 = (1.d0 - rho)*sigmai*pp*bik1
                      bstar2 = (1.d0 - rho)*sigmai*pp*bik2
                      
                      bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                                        
                  endif
                  
                  summ = summ 
     +                 + (-bstar1**2 
     +                    + (rstar - bstar)*bstar2)/bsigma

              enddo
            
          enddo

      enddo
      
      asigma = summ - b4*sigmal

	asigma = -1.d0/asigma*2.d0
      if (asigma .lt. 0.0d0) asigma = -asigma      

      bpdf = -fsigmal(bold) 
     +     + (bold - amean)**2/(2.d0*asigma)

      do ii = 1, 50
	  
          call rnset(iseed)
		rv = drnnof()
		call rnget(iseed)
		anew = amean + rv*dsqrt(asigma)

		apdf = -fsigmal(anew) 
     +         + (anew - amean)**2/(2.d0*asigma)
          
		ratio = apdf - bpdf 
     
          if (ratio .ge. 0.0d0) then
		    bpdf = apdf
			bold = anew
	    else
		    call rnset(iseed)
			u = drnunf()
			call rnget(iseed)
			if (dlog(u) .le. ratio) then
			    bpdf = apdf
				bold = anew
			endif
		endif
		
      enddo
            
      sigmal = dexp(bold)         
      
      end subroutine
      

      real*8 function fsigmal(sigma)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api
      real*8 a4, b4

      real*8 sigma, zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1, dd, bik
      real*8 temp, rstar, bstar, bsigma
      real*8 summ, pdf
                  
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vapi/api

      common /va4/a4
      common /vb4/b4
      
      sigmal = dexp(sigma)

      summ = 0.d0      
      do i = 1, ni         

          zb = 0.d0
          do m = 1, nz
              zb = zb + zz(i,m)*zmu(m)
          enddo

          wb = 0.d0
          do m = 1, nw
              wb = wb + ww(i,m)*wbeta(m)
          enddo      
          sigmai = dexp(wb)

          xb = sigmaz*zeta(i)
          do m = 1, nx
              xb = xb + xx(i,m)*xbeta(m)
          enddo      
          pp = dexp(xb)/(1.d0 + dexp(xb))           
          pp1 = 1.d0 - pp

          ub = sigmal*lambda(i)
          do m = 1, nu
              ub = ub + uu(i,m)*ubeta(m)
          enddo      
                                            
          do k = 1, nk

              bik = 0.d0
              do l = 1, nl
                  eub = dexp(aalpha(l) + ub)
                  dd = eub/(1.d0 + eub)
                  temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                  bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
              enddo

              do j = 1, nd        
            
                  jj = (j - 1)*nk + k      

                  if (jj .eq. 1) then

                      rstar = bij(i,jj) 
     +                      - sigmar*rij(i,jj) 
     +                      - zb - sigmat*taui(i)

                      bstar = sigmai*pp*bik

                      bsigma = (sigmai*pp1)**2
                      
                  else

                      rstar = bij(i,jj) - rho*bij(i,jj-1) 
     +                      - sigmar*(rij(i,jj) - rho*rij(i,jj-1)) 
     +                      - (1.d0 - rho)
     +                      * (zb + sigmat*taui(i))

                      bstar = (1.d0 - rho)*sigmai*pp*bik

                      bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                                        
                  endif
                            
                  summ = summ - (rstar - bstar)**2/(2.d0*bsigma)
                                    
              enddo
              
          enddo
          
      enddo           
                             	          
      pdf = summ + a4*sigma - b4*dexp(sigma)

	fsigmal = -pdf         	    
     	    
	end function  

     
	subroutine gibbs_rho(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10
	integer, parameter :: npoint = 31, ndim = 3

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio

      real*8 epsilon, astar(npoint), diffa(npoint)
      real*8 yystar(npoint), xxstar(npoint,ndim)
      real*8 xxx(ndim,ndim), xy(ndim), amaxyy
      real*8 xhat, detxx, xinvy
                        
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

      external frho, drnnof, drnunf      

      epsilon = 0.1d0

      diffa(1) = -1.5d0
      do ii = 2, npoint
        diffa(ii) = diffa(ii-1) + 0.1d0
      enddo
     
      bold = dlog(rho/(1.d0 - rho))
            
      nopt = 1
	reqmin = 1.0d-10 ; konvge = 5 ; kcount = 1000
	step(1) = 0.2d0 ; start(1) = bold
      call nelmin(frho, nopt, start, xmin, 
     +            ynewlo, reqmin, step, konvge, 
     +            kcount, icount, numres, ifault)
	amean = xmin(1)
      
      do kk = 1, npoint
          astar(kk) = amean + diffa(kk)*epsilon
      enddo

      do kk = 1, npoint
          yystar(kk) = -frho(astar(kk))
          xxstar(kk,1) = astar(kk)**2
          xxstar(kk,2) = astar(kk)
          xxstar(kk,3) = 1.d0
      enddo

      amaxyy = dabs(yystar(1))
      do kk = 2, npoint
          if (dabs(yystar(kk)) .gt. amaxyy) then
              amaxyy = dabs(yystar(kk))
          endif
      enddo

      do kk = 1, npoint
          yystar(kk) = yystar(kk)/amaxyy
      enddo
        
      do k1 = 1, ndim
          do k2 = 1, ndim
              xxx(k1,k2) = 0.d0
              do kk = 1, npoint
                  xxx(k1,k2) = xxx(k1,k2) 
     +                       + xxstar(kk,k1)
     +                         *xxstar(kk,k2)
              enddo
          enddo
          xy(k1) = 0.d0
          do k2 = 1, npoint
              xy(k1) = xy(k1) + xxstar(k2,k1)*yystar(k2)
          enddo
      enddo
        
      detxx = xxx(1,1)*xxx(2,2)*xxx(3,3)
     +      + xxx(2,1)*xxx(3,2)*xxx(1,3)
     +      + xxx(3,1)*xxx(1,2)*xxx(2,3)
     +      - xxx(1,1)*xxx(3,2)*xxx(2,3)
     +      - xxx(3,1)*xxx(2,2)*xxx(1,3)
     +      - xxx(2,1)*xxx(1,2)*xxx(3,3)   
                                              
      xinvy = (xxx(2,2)*xxx(3,3) - xxx(2,3)*xxx(3,2))*xy(1) 
     +      + (xxx(1,3)*xxx(3,2) - xxx(1,2)*xxx(3,3))*xy(2)      
     +      + (xxx(1,2)*xxx(2,3) - xxx(1,3)*xxx(2,2))*xy(3)     

      xhat = xinvy/detxx
              
      asigma = -1.d0/(2.d0*amaxyy*xhat)*2.d0
      if (asigma .lt. 0.0d0) asigma = -asigma                

      bpdf = -frho(bold) 
     +     + (bold - amean)**2/(2.d0*asigma)

      do ii = 1, 50
	  
          call rnset(iseed)
		rv = drnnof()
		call rnget(iseed)
		anew = amean + rv*dsqrt(asigma)

		apdf = -frho(anew) 
     +         + (anew - amean)**2/(2.d0*asigma)
          
		ratio = apdf - bpdf 
     
          if (ratio .ge. 0.0d0) then
		    bpdf = apdf
			bold = anew
	    else
		    call rnset(iseed)
			u = drnunf()
			call rnget(iseed)
			if (dlog(u) .le. ratio) then
			    bpdf = apdf
				bold = anew
			endif
		endif
		
      enddo
            
      rho = dexp(bold)/(1.d0 + dexp(bold))
      
      end subroutine
      

      real*8 function frho(rstar)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api

      real*8 rstar, zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1, dd, bik
      real*8 temp, bstar, bsigma
      real*8 summ, pdf
                  
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

      rho = dexp(rstar)/(1.d0 + dexp(rstar))
            
      summ = 0.d0      
      do i = 1, ni         

          zb = 0.d0
          do m = 1, nz
              zb = zb + zz(i,m)*zmu(m)
          enddo

          wb = 0.d0
          do m = 1, nw
              wb = wb + ww(i,m)*wbeta(m)
          enddo      
          sigmai = dexp(wb)
 
          xb = sigmaz*zeta(i)
          do m = 1, nx
              xb = xb + xx(i,m)*xbeta(m)
          enddo      
          pp = dexp(xb)/(1.d0 + dexp(xb)) 
          pp1 = 1.d0 - pp

          ub = sigmal*lambda(i)
          do m = 1, nu
              ub = ub + uu(i,m)*ubeta(m)
          enddo      
                                   
          do k = 1, nk

              bik = 0.d0
              do l = 1, nl
                  eub = dexp(aalpha(l) + ub)
                  dd = eub/(1.d0 + eub)                     
                  temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                  bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
              enddo

              do j = 1, nd        
            
                  jj = (j - 1)*nk + k      

                  if (jj .gt. 1) then

                      bstar = bij(i,jj) - rho*bij(i,jj-1) 
     +                      - sigmar*(rij(i,jj) - rho*rij(i,jj-1)) 
     +                      - (1.d0 - rho)
     +                      * (zb + sigmat*taui(i) 
     +                         + sigmai*pp*bik)

                      bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2

                      summ = summ 
     +                     + dlog(1.d0 - rho**2)/2.d0 
     +                     + bstar**2/(2.d0*bsigma)
                  endif
                                
              enddo
            
          enddo
        
      enddo           
                       	          
      pdf = -summ/2.d0 
     +    + rstar - 2.d0*dlog(1.d0 + dexp(rstar))
      
	frho = -pdf    

	end function 
      

	subroutine gibbs_wbeta(iseed) 
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api
      real*8 sigma_wbeta

	real*8 start(nw), xmin(nw), ynewlo, reqmin, step(nw)
	integer konvge, kcount, icount, numres, ifault, nopt

      real*8 bold(nw), anew(nw)
      real*8 amean(nw), asigma(nw,nw), der2(nw,nw)
	real*8 tol, rsig(nw,nw), hmean(nw)
	real*8 hpdf, bpdf, apdf, ratio, u
	
      real*8 zb, wb, ub, eub
      real*8 xb, pp, pp1, dd, bik
      real*8 temp, bstar1, bstar2, bsigma
                  
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

      common /vsigma_wbeta/sigma_wbeta

      external fwbeta, dlinrg, dchfac, dblinf, drnmvn, drnunf

      do jj = 1, nw
        bold(jj) = wbeta(jj)
      enddo
      
	nopt = nw
	reqmin = 1.0d-10 ; konvge = 5 ; kcount = 1000
	do jj = 1, nw
	  step(jj) = 0.2d0 ; start(jj) = bold(jj)
	enddo	
	call nelmin(fwbeta, nopt, start, xmin, 
     +            ynewlo, reqmin, step, konvge, 
     +            kcount, icount, numres, ifault)       
      do jj = 1, nw
        amean(jj) = xmin(jj)
      enddo

      do j1 = 1, nw
          do j2 = 1, nw
              der2(j1,j2) = 0.d0
          enddo
      enddo

      do j1 = 1, nw
      do j2 = j1, nw

          do i = 1, ni   

              zb = 0.d0
              do m = 1, nz
                  zb = zb + zz(i,m)*zmu(m)
              enddo

              wb = 0.d0
              do m = 1, nw
                  wb = wb + ww(i,m)*amean(m)
              enddo      
            
              xb = sigmaz*zeta(i)
              do m = 1, nx
                  xb = xb + xx(i,m)*xbeta(m)
              enddo
              pp = dexp(xb)/(1.d0 + dexp(xb))           
              pp1 = 1.d0 - pp

              ub = sigmal*lambda(i)
              do m = 1, nu
                  ub = ub + uu(i,m)*ubeta(m)
              enddo      
                             
              do k = 1, nk

                  bik = 0.d0
                  do l = 1, nl
                      eub = dexp(aalpha(l) + ub)
                      dd = eub/(1.d0 + eub)                     
                      temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                      bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
                  enddo

                  do j = 1, nd        
            
                      jj = (j - 1)*nk + k                           

                      if (jj .eq. 1) then

                          bstar1 = bij(i,jj) 
     +                           - sigmar*rij(i,jj) 
     +                           - zb - sigmat*taui(i)
                      
                          bstar2 = pp*bik
                      
                          bsigma = pp1**2
                      
                      else

                          bstar1 = bij(i,jj) - rho*bij(i,jj-1) 
     +                           - sigmar*(rij(i,jj) 
     +                                     - rho*rij(i,jj-1)) 
     +                          - (1.d0 - rho)
     +                          * (zb + sigmat*taui(i)) 
                      
                          bstar2 = (1.d0 - rho)*pp*bik
                      
                          bsigma = (1.d0 - rho**2)*pp1**2
                                        
                      endif
                      
                      der2(j1,j2) = der2(j1,j2) 
     +                            + bstar1*dexp(-wb)
     +                            *( -2.d0*bstar1*dexp(-wb) 
     +                               + bstar2 )
     +                             /bsigma
     +                            *ww(i,j1)*ww(i,j2)
                                
                  enddo
                
              enddo
            
          enddo     

      enddo      
      enddo

      do j1 = 1, nw            
        der2(j1,j1) = der2(j1,j1) - 1.d0/sigma_wbeta        
      enddo
            
      do j1 = 1, nw
          do j2 = j1, nw
              der2(j2,j1) = der2(j1,j2)
          enddo
      enddo                

      do j1 = 1, nw
          do j2 = 1, nw
              der2(j1,j2) = -der2(j1,j2)*0.65d0
          enddo
      enddo                

	call dlinrg(nw, der2, nw, asigma, nw)	
	  	  
	tol = 100.0d0*dmach(4)
	call dchfac(nw, asigma, nw, tol, irank, rsig, nw)

	do jj = 1, nw
	    hmean(jj) = bold(jj) - amean(jj)
	enddo
	hpdf = dblinf(nw, nw, der2, nw, hmean, hmean)
	bpdf = -fwbeta(bold) + hpdf/2.d0

	do ii = 1, 50
	
	    call rnset(iseed)
		call drnmvn(1, nw, rsig, nw, anew, 1)
		call rnget(iseed)		
		do jj = 1, nw
	        anew(jj) = amean(jj) + anew(jj)
		enddo
		
		do jj = 1, nw
		    hmean(jj) = anew(jj) - amean(jj)
		enddo
		hpdf = dblinf(nw, nw, der2, nw, hmean, hmean)
		apdf = -fwbeta(anew) + hpdf/2.d0

		ratio = apdf - bpdf 

		if (ratio .ge. 0.0d0) then
		    bpdf = apdf
			do jj = 1, nw
				bold(jj) = anew(jj)
			enddo
		else
		    call rnset(iseed)
			u = drnunf()
			call rnget(iseed)
			if (dlog(u) .le. ratio) then
			    bpdf = apdf
				do jj = 1, nw
					bold(jj) = anew(jj)
				enddo
			endif
		endif
		
	enddo      

      do jj = 1, nw
        wbeta(jj) = bold(jj)
      enddo
      
      end subroutine
      

      real*8 function fwbeta(wbeta)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api
      real*8 sigma_wbeta

      real*8 zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1, dd, bik
      real*8 temp, bstar, bsigma
      real*8 sum1, sum2, pdf
                  
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

      common /vsigma_wbeta/sigma_wbeta

      sum1 = 0.d0      
      do i = 1, ni         

          zb = 0.d0
          do m = 1, nz
              zb = zb + zz(i,m)*zmu(m)
          enddo

          wb = 0.d0
          do m = 1, nw
              wb = wb + ww(i,m)*wbeta(m)
          enddo      
          sigmai = dexp(wb)

          xb = sigmaz*zeta(i)
          do m = 1, nx
              xb = xb + xx(i,m)*xbeta(m)
          enddo      
          pp = dexp(xb)/(1.d0 + dexp(xb))           
          pp1 = 1.d0 - pp

          ub = sigmal*lambda(i)
          do m = 1, nu
              ub = ub + uu(i,m)*ubeta(m)
          enddo      
                                   
          do k = 1, nk

              bik = 0.d0
              do l = 1, nl
                  eub = dexp(aalpha(l) + ub)
                  dd = eub/(1.d0 + eub)                     
                  temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                  bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
              enddo

              do j = 1, nd        
            
                  jj = (j - 1)*nk + k      

                  if (jj .eq. 1) then

                      bstar = bij(i,jj) 
     +                      - sigmar*rij(i,jj)
     +                      - zb - sigmat*taui(i) - sigmai*pp*bik
                      
                      bsigma = (sigmai*pp1)**2
                      
                  else

                      bstar = bij(i,jj) - rho*bij(i,jj-1) 
     +                      - sigmar*(rij(i,jj) - rho*rij(i,jj-1)) 
     +                      - (1.d0 - rho)
     +                      * (zb + sigmat*taui(i) + sigmai*pp*bik) 
                      
                      bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                                        
                  endif

                  sum1 = sum1 - dlog(sigmai) 
     +                 - bstar**2/(2.d0*bsigma)
                  
              enddo
            
          enddo
        
      enddo           
                       	    
      sum2 = 0.d0
      do jj = 1, nw
        sum2 = sum2 - wbeta(jj)**2/(2.d0*sigma_wbeta)
      enddo
      
      pdf = sum1 + sum2

	fwbeta = -pdf         	    
     	    
	end function  


	subroutine gibbs_ubeta(iseed) 
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api
      real*8 sigma_ubeta

	real*8 start(nu), xmin(nu), ynewlo, reqmin, step(nu)
	integer konvge, kcount, icount, numres, ifault, nopt

      real*8 bold(nu), anew(nu)
      real*8 amean(nu), asigma(nu,nu), der2(nu,nu)
	real*8 tol, rsig(nu,nu), hmean(nu)
	real*8 hpdf, bpdf, apdf, ratio, u

      real*8 zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1
      real*8 bik, bik1, bik2, dd, dd1, dd2
      real*8 temp, rstar, bstar, bstar1, bstar2, bsigma
                  
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

      common /vsigma_ubeta/sigma_ubeta

      external fubeta, dlinrg, dchfac, dblinf, drnmvn, drnunf

      do jj = 1, nu
        bold(jj) = ubeta(jj)
      enddo
      
	nopt = nu
	reqmin = 1.0d-10 ; konvge = 5 ; kcount = 1000
	do jj = 1, nu
	    step(jj) = 0.2d0 ; start(jj) = bold(jj)
	enddo	
	call nelmin(fubeta, nopt, start, xmin, 
     +            ynewlo, reqmin, step, konvge, 
     +            kcount, icount, numres, ifault)       
      do jj = 1, nu
          amean(jj) = xmin(jj)
      enddo

      do j1 = 1, nu
          do j2 = 1, nu
              der2(j1,j2) = 0.d0
          enddo
      enddo

      do j1 = 1, nu
      do j2 = j1, nu

          do i = 1, ni         

              zb = 0.d0
              do m = 1, nz
                  zb = zb + zz(i,m)*zmu(m)
              enddo

              wb = 0.d0
              do m = 1, nw
                  wb = wb + ww(i,m)*wbeta(m)
              enddo      
              sigmai = dexp(wb)
 
              xb = sigmaz*zeta(i)
              do m = 1, nx
                  xb = xb + xx(i,m)*xbeta(m)
              enddo      
              pp = dexp(xb)/(1.d0 + dexp(xb))           
              pp1 = 1.d0 - pp

              ub = sigmal*lambda(i)
              do m = 1, nu
                  ub = ub + uu(i,m)*amean(m)
              enddo     
                                  
              do k = 1, nk

                  bik = 0.d0; bik1 = 0.d0; bik2 = 0.d0 
                  do l = 1, nl

                      eub = dexp(aalpha(l) + ub)
                      
                      dd = eub/(1.d0 + eub)                     
                      dd1 = eub/(1.d0 + eub)**2
                      dd2 = eub*(1.d0 - eub)
     +                    /(1.d0 + eub)**3

                      temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                      bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
                      
                      bik1 = bik1 
     +                     - dsqrt(2.d0*delta(l))*dsin(temp)
     +                       *2.d0*api*dfloat(l)*dd1
                      
                      bik2 = bik2 
     +                - dsqrt(2.d0*delta(l))
     +                *(dcos(temp)*(2.d0*api*dfloat(l)*dd1)**2 
     +                  + dsin(temp)*2.d0*api*dfloat(l)*dd2)
                      
                  enddo

                  do j = 1, nd        
            
                      jj = (j - 1)*nk + k      

                      if (jj .eq. 1) then

                          rstar = bij(i,jj) 
     +                          - sigmar*rij(i,jj) 
     +                          - zb - sigmat*taui(i)
                      
                          bstar = sigmai*pp*bik
                          bstar1 = sigmai*pp*bik1
                          bstar2 = sigmai*pp*bik2
                      
                          bsigma = (sigmai*pp1)**2
                      
                      else

                          rstar = bij(i,jj) - rho*bij(i,jj-1) 
     +                          - sigmar*(rij(i,jj) 
     +                                    - rho*rij(i,jj-1)) 
     +                          - (1.d0 - rho)
     +                          * (zb + sigmat*taui(i)) 
                      
                          bstar = (1.d0 - rho)*sigmai*pp*bik
                          bstar1 = (1.d0 - rho)*sigmai*pp*bik1
                          bstar2 = (1.d0 - rho)*sigmai*pp*bik2
                      
                          bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                                        
                      endif
                      
                      der2(j1,j2) = der2(j1,j2) 
     +                            + ( -bstar1**2 
     +                                + (rstar - bstar)*bstar2 )
     +                              /bsigma
     +                              *uu(i,j1)*uu(i,j2)

                  enddo
                
              enddo

          enddo
                 
      enddo
      enddo

      do j1 = 1, nu            
          der2(j1,j1) = der2(j1,j1) - 1.d0/sigma_ubeta        
      enddo

      do j1 = 1, nu
          do j2 = j1, nu
              der2(j2,j1) = der2(j1,j2)
          enddo
      enddo                

      do j1 = 1, nu
          do j2 = 1, nu
              der2(j1,j2) = -der2(j1,j2)*0.65d0
          enddo
      enddo                

	call dlinrg(nu, der2, nu, asigma, nu)	
	  	  
	tol = 100.0d0*dmach(4)
	call dchfac(nu, asigma, nu, tol, irank, rsig, nu)

	do jj = 1, nu
	    hmean(jj) = bold(jj) - amean(jj)
	enddo
	hpdf = dblinf(nu, nu, der2, nu, hmean, hmean)
	bpdf = -fubeta(bold) + hpdf/2.d0

	do ii = 1, 50
	
          call rnset(iseed)
		call drnmvn(1, nu, rsig, nu, anew, 1)
		call rnget(iseed)		
		do jj = 1, nu
	        anew(jj) = amean(jj) + anew(jj)
		enddo
		
		do jj = 1, nu
		    hmean(jj) = anew(jj) - amean(jj)
		enddo
		hpdf = dblinf(nu, nu, der2, nu, hmean, hmean)
		apdf = -fubeta(anew) + hpdf/2.d0

		ratio = apdf - bpdf 

		if (ratio .ge. 0.0d0) then
		    bpdf = apdf
			do jj = 1, nu
				bold(jj) = anew(jj)
			enddo
		else
		    call rnset(iseed)
			u = drnunf()
			call rnget(iseed)
			if (dlog(u) .le. ratio) then
			    bpdf = apdf
				do jj = 1, nu
					bold(jj) = anew(jj)
				enddo
			endif
		endif
		
	enddo      

      do jj = 1, nu
        ubeta(jj) = bold(jj)
      enddo

      end subroutine
      

      real*8 function fubeta(ubeta)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api
      real*8 sigma_ubeta

      real*8 zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1, dd, bik
      real*8 temp, rstar, bstar, bsigma
      real*8 sum1, sum2, pdf
                  
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

      common /vsigma_ubeta/sigma_ubeta

      sum1 = 0.d0      
      do i = 1, ni         

          zb = 0.d0
          do m = 1, nz
              zb = zb + zz(i,m)*zmu(m)
          enddo

          wb = 0.d0
          do m = 1, nw
              wb = wb + ww(i,m)*wbeta(m)
          enddo      
          sigmai = dexp(wb)

          xb = sigmaz*zeta(i)
          do m = 1, nx
              xb = xb + xx(i,m)*xbeta(m)
          enddo      
          pp = dexp(xb)/(1.d0 + dexp(xb))           
          pp1 = 1.d0 - pp

          ub = sigmal*lambda(i)
          do m = 1, nu
              ub = ub + uu(i,m)*ubeta(m)
          enddo      
                                   
          do k = 1, nk

              bik = 0.d0
              do l = 1, nl
                  eub = dexp(aalpha(l) + ub)
                  dd = eub/(1.d0 + eub)                     
                  temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                  bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
              enddo

              do j = 1, nd        
            
                  jj = (j - 1)*nk + k      

                  if (jj .eq. 1) then

                      rstar = bij(i,jj) 
     +                      - sigmar*rij(i,jj)
     +                      - zb - sigmat*taui(i)
                      
                      bstar = sigmai*pp*bik

                      bsigma = (sigmai*pp1)**2
                      
                  else

                      rstar = bij(i,jj) - rho*bij(i,jj-1) 
     +                      - sigmar*(rij(i,jj) - rho*rij(i,jj-1)) 
     +                      - (1.d0 - rho)
     +                      * (zb + sigmat*taui(i)) 
                      
                      bstar = (1.d0 - rho)*sigmai*pp*bik

                      bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                                        
                  endif
                  
                  sum1 = sum1 - (rstar - bstar)**2/(2.d0*bsigma)

              enddo
            
          enddo
        
      enddo           
                       	    
      sum2 = 0.d0
      do jj = 1, nu
          sum2 = sum2 - ubeta(jj)**2/(2.d0*sigma_ubeta)
      enddo
      
      pdf = sum1 + sum2

	fubeta = -pdf         	    
     	    
	end function  
      
            
	subroutine gibbs_sigmaz(iseed) 
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api
      real*8 a3, b3

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio
      
      real*8 zb, wb, ub, eub, sigmai
      real*8 xb, exb, pp, dd, bik
      real*8 temp1, temp2
      real*8 bstarA, bstarC, bsigma
      real*8 summ
                  
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

      common /va3/a3
      common /vb3/b3

      external fsigmaz, drnnof, drnunf      
     
      bold = dlog(sigmaz)
      
      nopt = 1
	reqmin = 1.0d-10 ; konvge = 5 ; kcount = 1000
	step(1) = 0.2d0 ; start(1) = bold
      call nelmin(fsigmaz, nopt, start, xmin, 
     +            ynewlo, reqmin, step, konvge, 
     +            kcount, icount, numres, ifault)
	amean = xmin(1)      

      sigmaz = dexp(amean)

      summ = 0.d0      
      do i = 1, ni         

          zb = 0.d0
          do m = 1, nz
              zb = zb + zz(i,m)*zmu(m)
          enddo

          wb = 0.d0
          do m = 1, nw
              wb = wb + ww(i,m)*wbeta(m)
          enddo      
          sigmai = dexp(wb)

          xb = sigmaz*zeta(i)
          do m = 1, nx
              xb = xb + xx(i,m)*xbeta(m)
          enddo      
          exb = dexp(xb)
          pp = dexp(xb)/(1.d0 + dexp(xb))           

          ub = sigmal*lambda(i)
          do m = 1, nu
              ub = ub + uu(i,m)*ubeta(m)
          enddo      
        
          temp1 = sigmaz*zeta(i)           
                                   
          do k = 1, nk

              bik = 0.d0
              do l = 1, nl
                  eub = dexp(aalpha(l) + ub)
                  dd = eub/(1.d0 + eub)                     
                  temp2 = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                  bik = bik + dsqrt(2.d0*delta(l))*dcos(temp2)
              enddo

              do j = 1, nd        
            
                  jj = (j - 1)*nk + k      

                  if (jj .eq. 1) then

                      bstarA = bij(i,jj) 
     +                       - sigmar*rij(i,jj) 
     +                       - zb - sigmat*taui(i)

                      bstarC = bstarA
     +                       - sigmai*bik

                      bsigma = sigmai**2
                      
                  else

                      bstarA = bij(i,jj) - rho*bij(i,jj-1) 
     +                       - sigmar*(rij(i,jj) - rho*rij(i,jj-1)) 
     +                       - (1.d0 - rho)
     +                       * (zb + sigmat*taui(i))

                      bstarC = bstarA
     +                       - (1.d0 - rho)*sigmai*bik

                      bsigma = (1.d0 - rho**2)*sigmai**2
                                        
                  endif
                                
                  summ = summ + temp1*exb
     +                 * ( (1.d0 + temp1)*(1.d0 + exb) 
     +                    - exb*temp1 )/(1.d0 + exb)**2
     +                 - bstarC*temp1*exb
     +                 * ( bstarA*(1.d0 + temp1) 
     +                     + bstarC*exb*(1.d0 + 2.d0*temp1) )
     +                 / bsigma

              enddo
            
          enddo
        
      enddo           
                       	          
      asigma = summ - b3*sigmaz

	asigma = -1.d0/asigma*2.d0
      if (asigma .lt. 0.0d0) asigma = -asigma      

      bpdf = -fsigmaz(bold) 
     +     + (bold - amean)**2/(2.d0*asigma)

      do ii = 1, 50
	  
          call rnset(iseed)
		rv = drnnof()
		call rnget(iseed)
		anew = amean + rv*dsqrt(asigma)

		apdf = -fsigmaz(anew) 
     +         + (anew - amean)**2/(2.d0*asigma)
          
		ratio = apdf - bpdf 
     
          if (ratio .ge. 0.0d0) then
		    bpdf = apdf
			bold = anew
	    else
		    call rnset(iseed)
			u = drnunf()
			call rnget(iseed)
			if (dlog(u) .le. ratio) then
			    bpdf = apdf
				bold = anew
			endif
		endif
		
      enddo
            
      sigmaz = dexp(bold)         
      
      end subroutine
      

      real*8 function fsigmaz(sigma)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api
      real*8 a3, b3

      real*8 sigma, zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1, dd, bik
      real*8 temp, bstar, bsigma
      real*8 summ, pdf
                  
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

      common /va3/a3
      common /vb3/b3
      
      sigmaz = dexp(sigma)

      summ = 0.d0      
      do i = 1, ni         

          zb = 0.d0
          do m = 1, nz
              zb = zb + zz(i,m)*zmu(m)
          enddo

          wb = 0.d0
          do m = 1, nw
              wb = wb + ww(i,m)*wbeta(m)
          enddo      
          sigmai = dexp(wb)
 
          xb = sigmaz*zeta(i)
          do m = 1, nx
              xb = xb + xx(i,m)*xbeta(m)
          enddo      
          pp = dexp(xb)/(1.d0 + dexp(xb)) 
          pp1 = 1.d0 - pp

          ub = sigmal*lambda(i)
          do m = 1, nu
              ub = ub + uu(i,m)*ubeta(m)
          enddo      
                                   
          do k = 1, nk

              bik = 0.d0
              do l = 1, nl
                  eub = dexp(aalpha(l) + ub)
                  dd = eub/(1.d0 + eub)                     
                  temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                  bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
              enddo

              do j = 1, nd        
            
                  jj = (j - 1)*nk + k      

                  if (jj .eq. 1) then

                      bstar = bij(i,jj) 
     +                      - sigmar*rij(i,jj)
     +                      - zb - sigmat*taui(i) 
     +                      - sigmai*pp*bik

                      bsigma = (sigmai*pp1)**2
                      
                  else

                      bstar = bij(i,jj) - rho*bij(i,jj-1) 
     +                      - sigmar*(rij(i,jj) - rho*rij(i,jj-1)) 
     +                      - (1.d0 - rho)
     +                      * (zb + sigmat*taui(i) 
     +                         + sigmai*pp*bik)

                      bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                                        
                  endif
              
                  summ = summ 
     +                 - dlog(pp1) 
     +                 - bstar**2/(2.d0*bsigma)
                  
              enddo
            
          enddo
        
      enddo           
                       	          
      pdf = summ + a3*sigma - b3*dexp(sigma)

	fsigmaz = -pdf         	    
     	    
      end function  
         

	subroutine gibbs_xbeta(iseed) 
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api
      real*8 sigma_xbeta

	real*8 start(nx), xmin(nx), ynewlo, reqmin, step(nx)
	integer konvge, kcount, icount, numres, ifault, nopt

      real*8 bold(nx), anew(nx)
      real*8 amean(nx), asigma(nx,nx), der2(nx,nx)
	real*8 tol, rsig(nx,nx), hmean(nx)
	real*8 hpdf, bpdf, apdf, ratio, u

      real*8 zb, wb, ub, eub, sigmai
      real*8 xb, exb, pp, dd, bik
      real*8 temp, bstarA, bstarC, bsigma
                  
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

      common /vsigma_xbeta/sigma_xbeta

      external fxbeta, dlinrg, dchfac, dblinf, drnmvn, drnunf

      do jj = 1, nx
          bold(jj) = xbeta(jj)
      enddo
      
	nopt = nx
	reqmin = 1.0d-10 ; konvge = 5 ; kcount = 1000
	do jj = 1, nx
          step(jj) = 0.2d0 ; start(jj) = bold(jj)
	enddo	
	call nelmin(fxbeta, nopt, start, xmin, 
     +            ynewlo, reqmin, step, konvge, 
     +            kcount, icount, numres, ifault)       
      do jj = 1, nx
          amean(jj) = xmin(jj)
      enddo

      do j1 = 1, nx
          do j2 = 1, nx
              der2(j1,j2) = 0.d0
          enddo
      enddo

      do j1 = 1, nx
      do j2 = j1, nx

          do i = 1, ni   

              zb = 0.d0
              do m = 1, nz
                  zb = zb + zz(i,m)*zmu(m)
              enddo

              wb = 0.d0
              do m = 1, nw
                  wb = wb + ww(i,m)*wbeta(m)
              enddo      
              sigmai = dexp(wb)
            
              xb = sigmaz*zeta(i)
              do m = 1, nx
                  xb = xb + xx(i,m)*amean(m)
              enddo
              exb = dexp(xb)
              pp = dexp(xb)/(1.d0 + dexp(xb))           

              ub = sigmal*lambda(i)
              do m = 1, nu
                  ub = ub + uu(i,m)*ubeta(m)
              enddo      

              do k = 1, nk

                  bik = 0.d0
                  do l = 1, nl
                      eub = dexp(aalpha(l) + ub)
                      dd = eub/(1.d0 + eub)                     
                      temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                      bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
                  enddo

                  do j = 1, nd        
            
                      jj = (j - 1)*nk + k                           
                      
                      if (jj .eq. 1) then

                          bstarA = bij(i,jj) 
     +                           - sigmar*rij(i,jj) 
     +                           - zb - sigmat*taui(i)

                          bstarC = bstarA
     +                           - sigmai*bik

                          bsigma = sigmai**2
                      
                      else

                          bstarA = bij(i,jj) - rho*bij(i,jj-1) 
     +                           - sigmar*(rij(i,jj) 
     +                                     - rho*rij(i,jj-1)) 
     +                           - (1.d0 - rho)
     +                           * (zb + sigmat*taui(i))

                          bstarC = bstarA
     +                           - (1.d0 - rho)*sigmai*bik

                          bsigma = (1.d0 - rho**2)*sigmai**2
                                        
                      endif

                      der2(j1,j2) = der2(j1,j2) 
     +                            + ( exb/(1.d0 + exb)**2
     +                                - bstarC*exb
     +                                  *(bstarA + 2.d0*bstarC*exb)
     +                                  /bsigma )
     +                            * xx(i,j1)*xx(i,j2)
                                
                  enddo
                
              enddo
            
          enddo     

      enddo
      enddo

      do j1 = 1, nx            
          der2(j1,j1) = der2(j1,j1) - 1.d0/sigma_xbeta        
      enddo
            
      do j1 = 1, nx
          do j2 = j1, nx
              der2(j2,j1) = der2(j1,j2)
          enddo
      enddo                

      do j1 = 1, nx
          do j2 = 1, nx
              der2(j1,j2) = -der2(j1,j2)*0.65d0
          enddo
      enddo                

	call dlinrg(nx, der2, nx, asigma, nx)	
	  	  
	tol = 100.0d0*dmach(4)
	call dchfac(nx, asigma, nx, tol, irank, rsig, nx)

	do jj = 1, nx
	    hmean(jj) = bold(jj) - amean(jj)
	enddo
	hpdf = dblinf(nx, nx, der2, nx, hmean, hmean)
	bpdf = -fxbeta(bold) + hpdf/2.d0

	do ii = 1, 50
	
          call rnset(iseed)
		call drnmvn(1, nx, rsig, nx, anew, 1)
		call rnget(iseed)		
		do jj = 1, nx
	        anew(jj) = amean(jj) + anew(jj)
		enddo
		
		do jj = 1, nx
		    hmean(jj) = anew(jj) - amean(jj)
		enddo
		hpdf = dblinf(nx, nx, der2, nx, hmean, hmean)
		apdf = -fxbeta(anew) + hpdf/2.d0

		ratio = apdf - bpdf 

		if (ratio .ge. 0.0d0) then
		    bpdf = apdf
			do jj = 1, nx
				bold(jj) = anew(jj)
			enddo
		else
		    call rnset(iseed)
			u = drnunf()
			call rnget(iseed)
			if (dlog(u) .le. ratio) then
			    bpdf = apdf
				do jj = 1, nx
					bold(jj) = anew(jj)
				enddo
			endif
		endif
		
	enddo      

      do jj = 1, nx
        xbeta(jj) = bold(jj)
      enddo

      end subroutine
      

      real*8 function fxbeta(xbeta)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api
      real*8 sigma_xbeta

      real*8 zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1, dd, bik
      real*8 temp, bstar, bsigma
      real*8 sum1, sum2, pdf
                  
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

      common /vsigma_xbeta/sigma_xbeta

      sum1 = 0.d0      
      do i = 1, ni         

          zb = 0.d0
          do m = 1, nz
              zb = zb + zz(i,m)*zmu(m)
          enddo

          wb = 0.d0
          do m = 1, nw
              wb = wb + ww(i,m)*wbeta(m)
          enddo      
          sigmai = dexp(wb)

          xb = sigmaz*zeta(i)
          do m = 1, nx
              xb = xb + xx(i,m)*xbeta(m)
          enddo      
          pp = dexp(xb)/(1.d0 + dexp(xb)) 
          pp1 = 1.d0 - pp

          ub = sigmal*lambda(i)
          do m = 1, nu
              ub = ub + uu(i,m)*ubeta(m)
          enddo      
                                    
          do k = 1, nk

              bik = 0.d0
              do l = 1, nl
                  eub = dexp(aalpha(l) + ub)
                  dd = eub/(1.d0 + eub)                     
                  temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                  bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
              enddo

              do j = 1, nd        
            
                  jj = (j - 1)*nk + k      

                  if (jj .eq. 1) then

                      bstar = bij(i,jj) 
     +                      - sigmar*rij(i,jj) 
     +                      - zb - sigmat*taui(i) 
     +                      - sigmai*pp*bik

                      bsigma = (sigmai*pp1)**2
                      
                  else

                      bstar = bij(i,jj) - rho*bij(i,jj-1) 
     +                      - sigmar*(rij(i,jj) - rho*rij(i,jj-1)) 
     +                      - (1.d0 - rho)
     +                      * (zb + sigmat*taui(i) 
     +                         + sigmai*pp*bik)

                      bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                                        
                  endif
              
                  sum1 = sum1 - dlog(pp1) - bstar**2/(2.d0*bsigma)
                  
              enddo
            
          enddo
        
      enddo           
                       	    
      sum2 = 0.d0
      do jj = 1, nx
        sum2 = sum2 - xbeta(jj)**2/(2.d0*sigma_xbeta)
      enddo
      
      pdf = sum1 + sum2

	fxbeta = -pdf         	    
     	    
      end function  
      
      
	subroutine gibbs_zmu(iseed)  
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api
      real*8 sigma_zmu

      real*8 wb, ub, eub, sigmai
      real*8 xb, pp, pp1, dd, bik
      real*8 bstar, bsigma
      real*8 temp

      real*8 hmean(nz), hsigma(nz,nz)
      real*8 amean(nz), asigma(nz,nz)
	real*8 tol, rsig(nz,nz), rv(nz)
      
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api
      
      common /vsigma_zmu/sigma_zmu

      external dlinrg, dmach, dchfac, drnmvn

      do j1 = 1, nz
          hmean(j1) = 0.d0
          do j2 = 1, nz
              hsigma(j1,j2) = 0.d0
          enddo
      enddo

      do i = 1, ni                      

          wb = 0.d0
          do m = 1, nw
              wb = wb + ww(i,m)*wbeta(m)
          enddo      
          sigmai = dexp(wb)

          xb = sigmaz*zeta(i)
          do m = 1, nx
              xb = xb + xx(i,m)*xbeta(m)
          enddo      
          pp = dexp(xb)/(1.d0 + dexp(xb))           
          pp1 = 1.d0 - pp

          ub = sigmal*lambda(i)
          do m = 1, nu
              ub = ub + uu(i,m)*ubeta(m)
          enddo      
                      
          do k = 1, nk

              bik = 0.d0
              do l = 1, nl
                  eub = dexp(aalpha(l) + ub)
                  dd = eub/(1.d0 + eub)                     
                  temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                  bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
              enddo

              do j = 1, nd        
            
                  jj = (j - 1)*nk + k      

                  if (jj .eq. 1) then

                      bstar = bij(i,jj) 
     +                      - sigmar*rij(i,jj)
     +                      - sigmat*taui(i) 
     +                      - sigmai*pp*bik

                      bsigma = (sigmai*pp1)**2
                  
                      do j1 = 1, nz
                      
                          hmean(j1) = hmean(j1) 
     +                              + zz(i,j1)*bstar/bsigma
                      
                          do j2 = 1, nz
                          
                              hsigma(j1,j2) = hsigma(j1,j2) 
     +                                      + zz(i,j1)*zz(i,j2)
     +                                        /bsigma
                          
                          enddo
                      
                      enddo           
                      
                  else

                      bstar = bij(i,jj) - rho*bij(i,jj-1) 
     +                      - sigmar*(rij(i,jj) - rho*rij(i,jj-1)) 
     +                      - (1.d0 - rho)*(sigmat*taui(i) 
     +                                      + sigmai*pp*bik)

                      bsigma = (1.d0 + rho)*(sigmai*pp1)**2
                  
                      do j1 = 1, nz
                      
                          hmean(j1) = hmean(j1) 
     +                              + zz(i,j1)*bstar/bsigma
                      
                          do j2 = 1, nz
                          
                              hsigma(j1,j2) = hsigma(j1,j2) 
     +                                      + (1.d0 - rho)
     +                                        *zz(i,j1)*zz(i,j2)
     +                                        /bsigma
                          
                          enddo
                      
                      enddo           
                      
                  endif
                
              enddo
                             
          enddo
        
      enddo     

      do jj = 1, nz
        hsigma(jj,jj) = hsigma(jj,jj) + 1.d0/sigma_zmu
      enddo

	call dlinrg(nz, hsigma, nz, asigma, nz)
	
	do j1 = 1, nz
          amean(j1) = 0.d0
		do j2 = 1, nz
              amean(j1) = amean(j1) + asigma(j1,j2)*hmean(j2)
		enddo
	enddo

	tol = 100.d0*dmach(4)
	call dchfac(nz, asigma, nz, tol, irank, rsig, nz)
	call rnset(iseed)
	call drnmvn(1, nz, rsig, nz, rv, 1)
	call rnget(iseed)      
	
	do jj = 1, nz
          zmu(jj) = amean(jj) + rv(jj)
	enddo      
            
	end subroutine
      
      
	subroutine gibbs_zeta(iseed) 
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio

      integer idum
      real*8 zb, wb, ub, eub, sigmai
      real*8 xb, exb, pp, dd, bik
      real*8 temp, bstarA, bstarC, bsigma, summ
                  
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

	common /vidum/idum

      external fzeta, drnnof, drnunf
      
      do i = 1, ni

          idum = i

          bold = zeta(i)
      
          nopt = 1
	    reqmin = 1.0d-10 ; konvge = 5 ; kcount = 1000
	    step(1) = 0.2d0 ; start(1) = bold
          call nelmin(fzeta, nopt, start, xmin, ynewlo,  
     +                reqmin, step, konvge, kcount, 
     +                icount, numres, ifault)
	    amean = xmin(1)

          zeta(i) = amean

          zb = 0.d0
          do m = 1, nz
              zb = zb + zz(i,m)*zmu(m)
          enddo

          wb = 0.d0
          do m = 1, nw
              wb = wb + ww(i,m)*wbeta(m)
          enddo      
          sigmai = dexp(wb)
                    
          xb = sigmaz*zeta(i)
          do m = 1, nx
              xb = xb + xx(i,m)*xbeta(m)
          enddo
          exb = dexp(xb)
          pp = dexp(xb)/(1.d0 + dexp(xb))           
     
          ub = sigmal*lambda(i)
          do m = 1, nu
              ub = ub + uu(i,m)*ubeta(m)
          enddo      
        
          summ = 0.d0                     
          do k = 1, nk

              bik = 0.d0
              do l = 1, nl
                  eub = dexp(aalpha(l) + ub)
                  dd = eub/(1.d0 + eub)                     
                  temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                  bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
              enddo

              do j = 1, nd        
            
                  jj = (j - 1)*nk + k                           
                  
                  if (jj .eq. 1) then

                      bstarA = bij(i,jj) 
     +                       - sigmar*rij(i,jj)
     +                       - zb - sigmat*taui(i)

                      bstarC = bstarA
     +                       - sigmai*bik

                      bsigma = sigmai**2
                      
                  else

                      bstarA = bij(i,jj) - rho*bij(i,jj-1) 
     +                       - sigmar*(rij(i,jj) - rho*rij(i,jj-1)) 
     +                       - (1.d0 - rho)
     +                       * (zb + sigmat*taui(i))

                      bstarC = bstarA
     +                       - (1.d0 - rho)*sigmai*bik

                      bsigma = (1.d0 - rho**2)*sigmai**2
                                        
                  endif
                  
                  summ = summ 
     +                 + ( exb/(1.d0 + exb)**2
     +                     - bstarC*exb
     +                       *(bstarA + 2.d0*bstarC*exb)
     +                       /bsigma )
     +                 * sigmaz**2
                                
              enddo 
                
          enddo

          asigma = summ - 1.d0
            
	    asigma = -1.d0/asigma*2.d0
          if (asigma .lt. 0.0d0) asigma = -asigma      

          bpdf = -fzeta(bold) 
     +         + (bold - amean)**2/(2.d0*asigma)

          do ii = 1, 50
	  
              call rnset(iseed)
		    rv = drnnof()
		    call rnget(iseed)
		    anew = amean + rv*dsqrt(asigma)

		    apdf = -fzeta(anew) 
     +             + (anew - amean)**2/(2.d0*asigma)
              
		    ratio = apdf - bpdf 
     
              if (ratio .ge. 0.0d0) then
		        bpdf = apdf
			    bold = anew
	        else
		        call rnset(iseed)
			    u = drnunf()
			    call rnget(iseed)
			    if (dlog(u) .le. ratio) then
			        bpdf = apdf
				    bold = anew
			    endif
		    endif
      		
          enddo
            
          zeta(i) = bold
        
      enddo

      end subroutine


      real*8 function fzeta(zstar)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api

      integer idum
      real*8 zstar, zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1, dd, bik
      real*8 temp, bstar, bsigma
      real*8 summ, pdf
                  
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

	common /vidum/idum

      i = idum
      
      zeta(i) = zstar

      zb = 0.d0
      do m = 1, nz
          zb = zb + zz(i,m)*zmu(m)
      enddo

      wb = 0.d0
      do m = 1, nw
          wb = wb + ww(i,m)*wbeta(m)
      enddo      
      sigmai = dexp(wb)
      
      xb = sigmaz*zeta(i)
      do m = 1, nx
          xb = xb + xx(i,m)*xbeta(m)
      enddo      
      pp = dexp(xb)/(1.d0 + dexp(xb)) 
      pp1 = 1.d0 - pp

      ub = sigmal*lambda(i)
      do m = 1, nu
          ub = ub + uu(i,m)*ubeta(m)
      enddo      

      summ = 0.d0                                   
      do k = 1, nk

          bik = 0.d0
          do l = 1, nl
              eub = dexp(aalpha(l) + ub)
              dd = eub/(1.d0 + eub)                     
              temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
              bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
          enddo

          do j = 1, nd        
            
              jj = (j - 1)*nk + k      

              if (jj .eq. 1) then

                  bstar = bij(i,jj) 
     +                  - sigmar*rij(i,jj)
     +                  - zb - sigmat*taui(i) 
     +                  - sigmai*pp*bik

                  bsigma = (sigmai*pp1)**2
                      
              else

                  bstar = bij(i,jj) - rho*bij(i,jj-1) 
     +                  - sigmar*(rij(i,jj) - rho*rij(i,jj-1)) 
     +                  - (1.d0 - rho)
     +                  * (zb + sigmat*taui(i) 
     +                     + sigmai*pp*bik)

                  bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                                        
              endif
              
              summ = summ - dlog(pp1) - bstar**2/(2.d0*bsigma)

          enddo
        
      enddo           
                       	          
      pdf = summ - zeta(i)**2/2.d0

	fzeta = -pdf         	    
     	    
	end function  

            
	subroutine gibbs_lambda(iseed) 
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio

      integer idum
      real*8 zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1
      real*8 bik, bik1, bik2, dd, dd1, dd2
      real*8 temp, rstar, bstar, bstar1, bstar2
      real*8 summ, bsigma
                  
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

	common /vidum/idum

      external flambda, drnnof, drnunf
      
      do i = 1, ni

          idum = i

          bold = lambda(i)
      
          nopt = 1
	    reqmin = 1.0d-10 ; konvge = 5 ; kcount = 1000
	    step(1) = 0.2d0 ; start(1) = bold
          call nelmin(flambda, nopt, start, xmin, ynewlo,  
     +                reqmin, step, konvge, kcount, 
     +                icount, numres, ifault)
	    amean = xmin(1)

          lambda(i) = amean
        
          zb = 0.d0
          do m = 1, nz
              zb = zb + zz(i,m)*zmu(m)
          enddo

          wb = 0.d0
          do m = 1, nw
              wb = wb + ww(i,m)*wbeta(m)
          enddo      
          sigmai = dexp(wb)

          xb = sigmaz*zeta(i)
          do m = 1, nx
              xb = xb + xx(i,m)*xbeta(m)
          enddo      
          pp = dexp(xb)/(1.d0 + dexp(xb))           
          pp1 = 1.d0 - pp

          ub = sigmal*lambda(i)
          do m = 1, nu
              ub = ub + uu(i,m)*ubeta(m)
          enddo     

          summ = 0.d0
          do k = 1, nk

              bik = 0.d0; bik1 = 0.d0; bik2 = 0.d0 
              do l = 1, nl

                  eub = dexp(aalpha(l) + ub)
          
                  dd = eub/(1.d0 + eub)                     
                  dd1 = eub/(1.d0 + eub)**2
                  dd2 = eub*(1.d0 - eub)
     +                /(1.d0 + eub)**3
                  
                  temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                  bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
                  
                  bik1 = bik1 
     +                 - dsqrt(2.d0*delta(l))*dsin(temp)
     +                   *2.d0*api*dfloat(l)*dd1
                  
                  bik2 = bik2 
     +                 - dsqrt(2.d0*delta(l))
     +                 *(dcos(temp)*(2.d0*api*dfloat(l)*dd1)**2 
     +                   + dsin(temp)*2.d0*api*dfloat(l)*dd2)
              enddo

              do j = 1, nd        
            
                  jj = (j - 1)*nk + k      
                  
                  if (jj .eq. 1) then

                      rstar = bij(i,jj) 
     +                      - sigmar*rij(i,jj)
     +                      - zb - sigmat*taui(i)

                      bstar = sigmai*pp*bik
                      bstar1 = sigmai*pp*bik1
                      bstar2 = sigmai*pp*bik2
                      
                      bsigma = (sigmai*pp1)**2
                      
                  else

                      rstar = bij(i,jj) - rho*bij(i,jj-1) 
     +                      - sigmar*(rij(i,jj) - rho*rij(i,jj-1)) 
     +                      - (1.d0 - rho)
     +                      * (zb + sigmat*taui(i))

                      bstar = (1.d0 - rho)*sigmai*pp*bik
                      bstar1 = (1.d0 - rho)*sigmai*pp*bik1
                      bstar2 = (1.d0 - rho)*sigmai*pp*bik2
                      
                      bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                                        
                  endif
                            
                  summ = summ 
     +                 + (-bstar1**2 
     +                     + (rstar - bstar)*bstar2)/bsigma
     +                    *sigmal**2
                  
              enddo
            
          enddo

          asigma = summ - 1.d0
            
	    asigma = -1.d0/asigma*2.d0
          if (asigma .lt. 0.0d0) asigma = -asigma      

          bpdf = -flambda(bold) 
     +         + (bold - amean)**2/(2.d0*asigma)

          do ii = 1, 50
	  
              call rnset(iseed)
		    rv = drnnof()
		    call rnget(iseed)
		    anew = amean + rv*dsqrt(asigma)

		    apdf = -flambda(anew) 
     +             + (anew - amean)**2/(2.d0*asigma)
              
		    ratio = apdf - bpdf 
     
              if (ratio .ge. 0.0d0) then
		        bpdf = apdf
			    bold = anew
	        else
		        call rnset(iseed)
			    u = drnunf()
			    call rnget(iseed)
			    if (dlog(u) .le. ratio) then
			        bpdf = apdf
				    bold = anew
			    endif
		    endif
      		
          enddo
            
          lambda(i) = bold
                  
      enddo
      
      end subroutine


      real*8 function flambda(astar)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api

      integer idum
      real*8 astar, zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1, dd, bik
      real*8 rstar, bstar, bsigma
      real*8 temp, summ, pdf
                  
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

	common /vidum/idum

      i = idum
      
      lambda(i) = astar

      zb = 0.d0
      do m = 1, nz
          zb = zb + zz(i,m)*zmu(m)
      enddo

      wb = 0.d0
      do m = 1, nw
          wb = wb + ww(i,m)*wbeta(m)
      enddo      
      sigmai = dexp(wb)
      
      xb = sigmaz*zeta(i)
      do m = 1, nx
          xb = xb + xx(i,m)*xbeta(m)
      enddo      
      pp = dexp(xb)/(1.d0 + dexp(xb))           
      pp1 = 1.d0 - pp

      ub = sigmal*lambda(i)
      do m = 1, nu
          ub = ub + uu(i,m)*ubeta(m)
      enddo      

      summ = 0.d0                                   
      do k = 1, nk

          bik = 0.d0
          do l = 1, nl
              eub = dexp(aalpha(l) + ub)
              dd = eub/(1.d0 + eub)                                   
              temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
              bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
          enddo
      
          do j = 1, nd        
            
              jj = (j - 1)*nk + k      

              if (jj .eq. 1) then

                  rstar = bij(i,jj) 
     +                  - sigmar*rij(i,jj) 
     +                  - zb - sigmat*taui(i)

                  bstar = sigmai*pp*bik

                  bsigma = (sigmai*pp1)**2
                      
              else

                  rstar = bij(i,jj) - rho*bij(i,jj-1) 
     +                  - sigmar*(rij(i,jj) - rho*rij(i,jj-1)) 
     +                  - (1.d0 - rho)
     +                  * (zb + sigmat*taui(i))

                  bstar = (1.d0 - rho)*sigmai*pp*bik

                  bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                                        
              endif
                            
              summ = summ - (rstar - bstar)**2/(2.d0*bsigma)

          enddo
          
      enddo           
                       	          
      pdf = summ - lambda(i)**2/2.d0

	flambda = -pdf         	    
     	    
      end function  
                  
      
	subroutine gibbs_sigmar(iseed)  
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api
      real*8 a1, b1

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio

      real*8 zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1, dd, bik
      real*8 bstar, rstar, bsigma
      real*8 temp, summ

      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

      common /va1/a1
      common /vb1/b1

      external fsigmar, drnnof, drnunf      

      bold = dlog(sigmar)
      
      nopt = 1
	reqmin = 1.0d-10 ; konvge = 5 ; kcount = 1000
	step(1) = 0.2d0 ; start(1) = bold
      call nelmin(fsigmar, nopt, start, xmin, 
     +            ynewlo, reqmin, step, konvge, 
     +            kcount, icount, numres, ifault)
	amean = xmin(1)      

      sigmar = dexp(amean)

      summ = 0.d0
      do i = 1, ni    

          zb = 0.d0
          do m = 1, nz
              zb = zb + zz(i,m)*zmu(m)
          enddo

          wb = 0.d0
          do m = 1, nw
              wb = wb + ww(i,m)*wbeta(m)
          enddo      
          sigmai = dexp(wb)

          xb = sigmaz*zeta(i)
          do m = 1, nx
              xb = xb + xx(i,m)*xbeta(m)
          enddo
          pp = dexp(xb)/(1.d0 + dexp(xb))           
          pp1 = 1.d0 - pp                      

          ub = sigmal*lambda(i)
          do m = 1, nu
              ub = ub + uu(i,m)*ubeta(m)
          enddo      
      
          do k = 1, nk

              bik = 0.d0
              do l = 1, nl
                  eub = dexp(aalpha(l) + ub)
                  dd = eub/(1.d0 + eub)                     
                  temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                  bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
              enddo

              do j = 1, nd        
            
                  jj = (j - 1)*nk + k      
                  
                  if (jj .eq. 1) then

                      bstar = bij(i,jj) 
     +                      - zb 
     +                      - sigmat*taui(i) 
     +                      - sigmai*pp*bik
                     
                      rstar = rij(i,jj) 
                      
                      bsigma = (sigmai*pp1)**2
                      
                  else

                      bstar = bij(i,jj) - rho*bij(i,jj-1) 
     +                      - (1.d0 - rho)
     +                      * (zb + sigmat*taui(i) 
     +                         + sigmai*pp*bik)
                     
                      rstar = rij(i,jj) - rho*rij(i,jj-1)
                      
                      bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                                        
                  endif
                                    
                  summ = summ + rstar*sigmar
     +                          *(bstar - 2.d0*rstar*sigmar)
     +                          /bsigma
                  
              enddo
        
          enddo
                
      enddo     

      asigma = summ - b1*sigmar

	asigma = -1.d0/asigma*2.d0
      if (asigma .lt. 0.0d0) asigma = -asigma      

      bpdf = -fsigmar(bold) 
     +     + (bold - amean)**2/(2.d0*asigma)

      do ii = 1, 50
	  
          call rnset(iseed)
		rv = drnnof()
		call rnget(iseed)
		anew = amean + rv*dsqrt(asigma)

		apdf = -fsigmar(anew) 
     +         + (anew - amean)**2/(2.d0*asigma)
          
		ratio = apdf - bpdf 
     
          if (ratio .ge. 0.0d0) then
		    bpdf = apdf
			bold = anew
	    else
		    call rnset(iseed)
			u = drnunf()
			call rnget(iseed)
			if (dlog(u) .le. ratio) then
			    bpdf = apdf
				bold = anew
			endif
		endif
		
      enddo
            
      sigmar = dexp(bold)                
          
      end subroutine


      real*8 function fsigmar(sigma)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api
      real*8 a1, b1
     
      real*8 sigma, zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1, dd, bik
      real*8 bstar, bmean
      real*8 temp, summ, pdf

      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

      common /va1/a1
      common /vb1/b1
      
      sigmar = dexp(sigma)

      summ = 0.d0
      do i = 1, ni    

          zb = 0.d0
          do m = 1, nz
              zb = zb + zz(i,m)*zmu(m)
          enddo

          wb = 0.d0
          do m = 1, nw
              wb = wb + ww(i,m)*wbeta(m)
          enddo      
          sigmai = dexp(wb)

          xb = sigmaz*zeta(i)
          do m = 1, nx
              xb = xb + xx(i,m)*xbeta(m)
          enddo
          pp = dexp(xb)/(1.d0 + dexp(xb))           
          pp1 = 1.d0 - pp                      

          ub = sigmal*lambda(i)
          do m = 1, nu
              ub = ub + uu(i,m)*ubeta(m)
          enddo      
      
          do k = 1, nk

              bik = 0.d0
              do l = 1, nl
                  eub = dexp(aalpha(l) + ub)
                  dd = eub/(1.d0 + eub)                     
                  temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                  bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
              enddo

              do j = 1, nd        
            
                  jj = (j - 1)*nk + k      
                                        
                  if (jj .eq. 1) then

                      bstar = bij(i,jj) 
     +                      - sigmar*rij(i,jj)
     +                      - zb 
     +                      - sigmat*taui(i) 
     +                      - sigmai*pp*bik

                      bsigma = (sigmai*pp1)**2
                      
                  else

                      bstar = bij(i,jj) - rho*bij(i,jj-1) 
     +                      - sigmar*(rij(i,jj) - rho*rij(i,jj-1))
     +                      - (1.d0 - rho)
     +                      * (zb + sigmat*taui(i) 
     +                         + sigmai*pp*bik)

                      bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                                        
                  endif
                                    
                  summ = summ - bstar**2/(2.d0*bsigma)
                  
              enddo
        
          enddo
                
      enddo     

      pdf = summ + a1*sigma - b1*dexp(sigma)

      fsigmar = -pdf
      
      end function  
                  
      
	subroutine gibbs_rij(iseed)   
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api
      
      real*8 zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1, dd, bik
      real*8 temp, bmean1, bmean2, bsigma

      real*8 hmean, hsigma, amean, asigma, rv

      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api
      
      external drnnof

      do i = 1, ni    

          zb = 0.d0
          do m = 1, nz
              zb = zb + zz(i,m)*zmu(m)
          enddo

          wb = 0.d0
          do m = 1, nw
              wb = wb + ww(i,m)*wbeta(m)
          enddo      
          sigmai = dexp(wb)

          xb = sigmaz*zeta(i)
          do m = 1, nx
              xb = xb + xx(i,m)*xbeta(m)
          enddo
          pp = dexp(xb)/(1.d0 + dexp(xb))           
          pp1 = 1.d0 - pp                      

          ub = sigmal*lambda(i)
          do m = 1, nu
              ub = ub + uu(i,m)*ubeta(m)
          enddo      
      
          do k = 1, nk

              bik = 0.d0
              do l = 1, nl
                  eub = dexp(aalpha(l) + ub)
                  dd = eub/(1.d0 + eub)                     
                  temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                  bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
              enddo

              do j = 1, nd        
            
                  jj = (j - 1)*nk + k      

                  if (jj .eq. 1) then

                      bmean1 = bij(i,jj)  
     +                       - zb - sigmat*taui(i) - sigmai*pp*bik

                      bmean2 = bij(i,jj+1) - rho*bij(i,jj) 
     +                       - (1.d0 - rho)*(zb + sigmat*taui(i) 
     +                                       + sigmai*pp*bik)

                      bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                  
                      hmean = sigmar*(bmean1/(sigmai*pp1)**2 
     +                                - rho*(bmean2 - sigmar
     +                                       *rij(i,jj+1))/bsigma)
     +                                            
                      hsigma = sigmar**2*(1.d0/(sigmai*pp1)**2 
     +                                    + rho**2/bsigma)
     +                       + 1.d0

                  else if (jj .eq. ndk) then

                      bmean1 = bij(i,jj) - rho*bij(i,jj-1) 
     +                       - (1.d0 - rho)*(zb + sigmat*taui(i) 
     +                                       + sigmai*pp*bik)

                      bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                  
                      hmean = sigmar*(bmean1 + sigmar*rho
     +                                         *rij(i,jj-1))
     +                      / bsigma
                      
                      hsigma = sigmar**2/bsigma
     +                       + 1.d0
                      
                  else

                      bmean1 = bij(i,jj) - rho*bij(i,jj-1) 
     +                       - (1.d0 - rho)*(zb + sigmat*taui(i) 
     +                                       + sigmai*pp*bik)

                      bmean2 = bij(i,jj+1) - rho*bij(i,jj) 
     +                       - (1.d0 - rho)*(zb + sigmat*taui(i) 
     +                                       + sigmai*pp*bik)

                      bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                  
                      hmean = sigmar*(bmean1 + sigmar*rho
     +                                         *rij(i,jj-1)
     +                                - rho*(bmean2 - sigmar
     +                                       *rij(i,jj+1)))
     +                      / bsigma
                      
                      hsigma = (1.d0 + rho**2)*sigmar**2/bsigma
     +                       + 1.d0
                      
                  endif

                  amean = hmean/hsigma   
                  asigma = 1.d0/hsigma

                  call rnset(iseed)
		        rv = drnnof()
		        call rnget(iseed)

		        rij(i,jj) = amean + rv*dsqrt(asigma)
                  
              enddo
        
          enddo
                          
      enddo     
             
      end subroutine
            
      
	subroutine gibbs_sigmat(iseed)  
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api
      real*8 a2, b2

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio

      real*8 zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1, dd, bik
      real*8 bstar, bsigma
      real*8 temp, summ

      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

      common /va2/a2
      common /vb2/b2

      external fsigmat, drnnof, drnunf      

      bold = dlog(sigmat)
      
      nopt = 1
	reqmin = 1.0d-10 ; konvge = 5 ; kcount = 1000
	step(1) = 0.2d0 ; start(1) = bold
      call nelmin(fsigmat, nopt, start, xmin, 
     +            ynewlo, reqmin, step, konvge, 
     +            kcount, icount, numres, ifault)
	amean = xmin(1)      

      sigmat = dexp(amean)

      summ = 0.d0
      do i = 1, ni    

          zb = 0.d0
          do m = 1, nz
              zb = zb + zz(i,m)*zmu(m)
          enddo

          wb = 0.d0
          do m = 1, nw
              wb = wb + ww(i,m)*wbeta(m)
          enddo      
          sigmai = dexp(wb)

          xb = sigmaz*zeta(i)
          do m = 1, nx
              xb = xb + xx(i,m)*xbeta(m)
          enddo
          pp = dexp(xb)/(1.d0 + dexp(xb))           
          pp1 = 1.d0 - pp                      

          ub = sigmal*lambda(i)
          do m = 1, nu
              ub = ub + uu(i,m)*ubeta(m)
          enddo      
      
          do k = 1, nk

              bik = 0.d0
              do l = 1, nl
                  eub = dexp(aalpha(l) + ub)
                  dd = eub/(1.d0 + eub)                     
                  temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                  bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
              enddo

              do j = 1, nd        
            
                  jj = (j - 1)*nk + k      
                  
                  if (jj .eq. 1) then

                      bstar = bij(i,jj)  
     +                      - sigmar*rij(i,jj)
     +                      - zb - 2.d0*sigmat*taui(i) 
     +                      - sigmai*pp*bik

                      bsigma = (sigmai*pp1)**2
                      
                  else

                      bstar = bij(i,jj) - rho*bij(i,jj-1) 
     +                      - sigmar*(rij(i,jj) - rho*rij(i,jj-1))
     +                      - (1.d0 - rho)
     +                      * (zb + 2.d0*sigmat*taui(i) 
     +                         + sigmai*pp*bik)
                     
                      bsigma = (1.d0 + rho)*(sigmai*pp1)**2
                                        
                  endif
                                    
                  summ = summ + sigmat*taui(i)*bstar/bsigma
                  
              enddo
        
          enddo
                
      enddo     

      asigma = summ - b2*sigmat

	asigma = -1.d0/asigma*2.d0
      if (asigma .lt. 0.0d0) asigma = -asigma      

      bpdf = -fsigmat(bold) 
     +     + (bold - amean)**2/(2.d0*asigma)

      do ii = 1, 50
	  
          call rnset(iseed)
		rv = drnnof()
		call rnget(iseed)
		anew = amean + rv*dsqrt(asigma)

		apdf = -fsigmat(anew) 
     +         + (anew - amean)**2/(2.d0*asigma)
          
		ratio = apdf - bpdf 
     
          if (ratio .ge. 0.0d0) then
		    bpdf = apdf
			bold = anew
	    else
		    call rnset(iseed)
			u = drnunf()
			call rnget(iseed)
			if (dlog(u) .le. ratio) then
			    bpdf = apdf
				bold = anew
			endif
		endif
		
      enddo
            
      sigmat = dexp(bold)                
          
      end subroutine


      real*8 function fsigmat(sigma)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api
      real*8 a2, b2
     
      real*8 sigma, zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1, dd, bik
      real*8 bstar, bsigma
      real*8 temp, summ, pdf

      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

      common /va2/a2
      common /vb2/b2
      
      sigmat = dexp(sigma)

      summ = 0.d0
      do i = 1, ni    

          zb = 0.d0
          do m = 1, nz
              zb = zb + zz(i,m)*zmu(m)
          enddo

          wb = 0.d0
          do m = 1, nw
              wb = wb + ww(i,m)*wbeta(m)
          enddo      
          sigmai = dexp(wb)

          xb = sigmaz*zeta(i)
          do m = 1, nx
              xb = xb + xx(i,m)*xbeta(m)
          enddo
          pp = dexp(xb)/(1.d0 + dexp(xb))           
          pp1 = 1.d0 - pp                      

          ub = sigmal*lambda(i)
          do m = 1, nu
              ub = ub + uu(i,m)*ubeta(m)
          enddo      
      
          do k = 1, nk

              bik = 0.d0
              do l = 1, nl
                  eub = dexp(aalpha(l) + ub)
                  dd = eub/(1.d0 + eub)                     
                  temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                  bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
              enddo

              do j = 1, nd        
            
                  jj = (j - 1)*nk + k      
                                        
                  if (jj .eq. 1) then

                      bstar = bij(i,jj) 
     +                      - sigmar*rij(i,jj)
     +                      - zb 
     +                      - sigmat*taui(i) 
     +                      - sigmai*pp*bik

                      bsigma = (sigmai*pp1)**2
                      
                  else

                      bstar = bij(i,jj) - rho*bij(i,jj-1) 
     +                      - sigmar*(rij(i,jj) - rho*rij(i,jj-1))
     +                      - (1.d0 - rho)
     +                      * (zb + sigmat*taui(i) 
     +                         + sigmai*pp*bik)

                      bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                                        
                  endif
                                    
                  summ = summ - bstar**2/(2.d0*bsigma)
                  
              enddo
        
          enddo
                
      enddo     

      pdf = summ + a2*sigma - b2*dexp(sigma)

      fsigmat = -pdf
      
      end function  
      

	subroutine gibbs_taui(iseed)   
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api
      
      real*8 zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1, dd, bik, bstar
      real*8 temp, bsigma

      real*8 hmean, hsigma, amean, asigma, rv

      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api
      
      external drnnof

      do i = 1, ni    

          zb = 0.d0
          do m = 1, nz
              zb = zb + zz(i,m)*zmu(m)
          enddo

          wb = 0.d0
          do m = 1, nw
              wb = wb + ww(i,m)*wbeta(m)
          enddo      
          sigmai = dexp(wb)

          xb = sigmaz*zeta(i)
          do m = 1, nx
              xb = xb + xx(i,m)*xbeta(m)
          enddo
          pp = dexp(xb)/(1.d0 + dexp(xb))           
          pp1 = 1.d0 - pp                      

          ub = sigmal*lambda(i)
          do m = 1, nu
              ub = ub + uu(i,m)*ubeta(m)
          enddo      
      
          hmean = 0.d0
          hsigma = 1.d0                      
          do k = 1, nk

              bik = 0.d0
              do l = 1, nl
                  eub = dexp(aalpha(l) + ub)
                  dd = eub/(1.d0 + eub)                     
                  temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
                  bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
              enddo

              do j = 1, nd        
            
                  jj = (j - 1)*nk + k      

                  if (jj .eq. 1) then

                      bstar = bij(i,jj)  
     +                      - sigmar*rij(i,jj)
     +                      - zb - sigmai*pp*bik

                      bsigma = (sigmai*pp1)**2
                  
                      hmean = hmean + sigmat*bstar/bsigma
                      hsigma = hsigma + sigmat**2/bsigma
                      
                  else

                      bstar = bij(i,jj) - rho*bij(i,jj-1) 
     +                      - sigmar*(rij(i,jj) - rho*rij(i,jj-1))
     +                      - (1.d0 - rho)*(zb + sigmai*pp*bik)

                      bsigma = (1.d0 + rho)*(sigmai*pp1)**2
                  
                      hmean = hmean + sigmat*bstar/bsigma
                      hsigma = hsigma 
     +                       + (1.d0 - rho)*sigmat**2/bsigma
                      
                  endif
                
              enddo
        
          enddo
        
          amean = hmean/hsigma   
          asigma = 1.d0/hsigma

          call rnset(iseed)
		rv = drnnof()
		call rnget(iseed)

		taui(i) = amean + rv*dsqrt(asigma)
                  
      enddo     
             
      end subroutine
      

	subroutine gibbs_bij(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 yy(ni,nd,nk), tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz    
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

	real*8 bold, anew, amean, asigma
	real*8 rv, u, bpdf, apdf, ratio

      integer idum, jdum, kdum, jjdum
      real*8 wb, sigmai, xb, pp, pp1, bsigma
                  
      common /vyy/yy
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta      
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

	common /vidum/idum
	common /vjdum/jdum
	common /vkdum/kdum
	common /vjjdum/jjdum

      external fbij, drnnof, drnunf      

      do i = 1, ni

          idum = i

          wb = 0.d0
          do m = 1, nw
              wb = wb + ww(i,m)*wbeta(m)
          enddo      
          sigmai = dexp(wb)          

          xb = sigmaz*zeta(i)
          do m = 1, nx
              xb = xb + xx(i,m)*xbeta(m)
          enddo
          pp = dexp(xb)/(1.d0 + dexp(xb))           
          pp1 = 1.d0 - pp
          
          do j = 1, nd

              jdum = j
        
              do k = 1, nk

                  kdum = k
            
                  jj = (j - 1)*nk + k

                  jjdum = jj


                  bold = bij(i,jj)
      
                  nopt = 1
	            reqmin = 1.0d-10 ; konvge = 5 ; kcount = 1000
	            step(1) = 0.2d0 ; start(1) = bold
                  call nelmin(fbij, nopt, start, xmin, ynewlo,  
     +                        reqmin, step, konvge, kcount, 
     +                        icount, numres, ifault)
	            amean = xmin(1)

                  if (jj .eq. 1) then

                      bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                      
                      asigma = -dexp(amean) 
     +                       - 1.d0/(sigmai*pp1)**2
     +                       - rho**2/bsigma
 
                  else if (jj .eq. ndk) then

                      bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                      
                      asigma = -dexp(amean) 
     +                       - 1.d0/bsigma
                
                  else
                      
                      bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
                      
                      asigma = -dexp(amean) 
     +                       - 1.d0/bsigma 
     +                       - rho**2/bsigma
                                                    
                  endif 
                                                
	            asigma = -1.d0/asigma*2.d0
                  if (asigma .lt. 0.0d0) asigma = -asigma      

                  bpdf = -fbij(bold) 
     +                 + (bold - amean)**2/(2.d0*asigma)

                  do ii = 1, 50
	  
                      call rnset(iseed)
		            rv = drnnof()
		            call rnget(iseed)
		            anew = amean + rv*dsqrt(asigma)

		            apdf = -fbij(anew) 
     +                     + (anew - amean)**2/(2.d0*asigma)
                      
                      ratio = apdf - bpdf 
     
                      if (ratio .ge. 0.0d0) then
		                bpdf = apdf
			            bold = anew
                          iiflag = 1
	                else
		                call rnset(iseed)
			            u = drnunf()
			            call rnget(iseed)
			            if (dlog(u) .le. ratio) then
			                bpdf = apdf
				            bold = anew
                              iiflag = 1
			            endif
		            endif
      		
                  enddo

                  bij(i,jj) = bold
                                             
              enddo
              
          enddo
          
      enddo     
                       
      end subroutine	


	subroutine gibbs_bij_mle(iseed)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 yy(ni,nd,nk), tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz    
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api

	real*8 start(1), xmin(1), ynewlo, reqmin, step(1)
	integer konvge, kcount, icount, numres, ifault, nopt

      integer idum, jdum, kdum, jjdum
                  
      common /vyy/yy
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta      
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

	common /vidum/idum
	common /vjdum/jdum
	common /vkdum/kdum
	common /vjjdum/jjdum

      external fbij      

      do i = 1, ni

          idum = i
          
          do j = 1, nd

              jdum = j
        
              do k = 1, nk

                  kdum = k
            
                  jj = (j - 1)*nk + k

                  jjdum = jj

                  nopt = 1
	            reqmin = 1.0d-10 ; konvge = 5 ; kcount = 1000
	            step(1) = 0.2d0 ; start(1) = bij(i,jj)
                  call nelmin(fbij, nopt, start, xmin, ynewlo,  
     +                        reqmin, step, konvge, kcount, 
     +                        icount, numres, ifault)
	            bij(i,jj) = xmin(1)
                                            
              enddo
              
          enddo
          
      enddo     
                       
      end subroutine	
      
      
      real*8 function fbij(bstar)
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10

	real*8 yy(ni,nd,nk), tik(ni,nk)
	real*8 zz(ni,nz), xx(ni,nx)
	real*8 ww(ni,nw), uu(ni,nu)

      real*8 taui(ni)
      real*8 rij(ni,ndk), bij(ni,ndk)
      real*8 zmu(nz), xbeta(nx)
      real*8 wbeta(nw), ubeta(nu)
      real*8 sigmar, sigmat, rho 
      real*8 zeta(ni), sigmaz    
      real*8 delta(nl), aalpha(nl)
      real*8 lambda(ni), sigmal, api

      integer idum, jdum, kdum, jjdum
      real*8 bstar, zb, wb, ub, eub, sigmai
      real*8 xb, pp, pp1, dd, bik
      real*8 bmean1, bmean2, bsigma
      real*8 temp, pdf

      common /vyy/yy
      common /vtik/tik

      common /vzz/zz
      common /vxx/xx
      common /vww/ww
      common /vuu/uu
      
      common /vtaui/taui
      common /vrij/rij
      common /vbij/bij

      common /vzmu/zmu
      common /vxbeta/xbeta
      common /vwbeta/wbeta
      common /vubeta/ubeta

      common /vsigmar/sigmar
      common /vsigmat/sigmat
      common /vrho/rho

      common /vzeta/zeta
      common /vsigmaz/sigmaz

      common /vdelta/delta      
      common /vaalpha/aalpha
      
      common /vlambda/lambda
      common /vsigmal/sigmal
      common /vapi/api

	common /vidum/idum
	common /vjdum/jdum
	common /vkdum/kdum
	common /vjjdum/jjdum

      i = idum
      j = jdum
      k = kdum
      jj = jjdum

      bij(i,jj) = bstar            

      zb = 0.d0
      do m = 1, nz
          zb = zb + zz(i,m)*zmu(m)
      enddo

      wb = 0.d0
      do m = 1, nw
          wb = wb + ww(i,m)*wbeta(m)
      enddo      
      sigmai = dexp(wb)
                  
      xb = sigmaz*zeta(i)
      do m = 1, nx
          xb = xb + xx(i,m)*xbeta(m)
      enddo
      pp = dexp(xb)/(1.d0 + dexp(xb))           
      pp1 = 1.d0 - pp

      ub = sigmal*lambda(i)
      do m = 1, nu
          ub = ub + uu(i,m)*ubeta(m)
      enddo      

      bik = 0.d0 
      do l = 1, nl
          eub = dexp(aalpha(l) + ub)
          dd = eub/(1.d0 + eub)                     
          temp = 2.d0*api*dfloat(l)*(tik(i,k) + dd)
          bik = bik + dsqrt(2.d0*delta(l))*dcos(temp)
      enddo
      
      if (jj .eq. 1) then

          bmean1 = sigmar*rij(i,jj)
     +           + zb + sigmat*taui(i) + sigmai*pp*bik

          bmean2 = sigmar*(rij(i,jj+1) - rho*rij(i,jj))
     +           + (1.d0 - rho)*(zb + sigmat*taui(i) 
     +                           + sigmai*pp*bik)
          
          bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
      
          pdf = yy(i,j,k)*bij(i,jj) - dexp(bij(i,jj))
     +        - (bij(i,jj) - bmean1)**2
     +        / (2.d0*(sigmai*pp1)**2)
     +        - (bij(i,jj+1) - bmean2 - rho*bij(i,jj))**2
     +        / (2.d0*bsigma)

      else if (jj .eq. ndk) then

          bmean1 = sigmar*(rij(i,jj) - rho*rij(i,jj-1))
     +           + (1.d0 - rho)*(zb + sigmat*taui(i) 
     +                           + sigmai*pp*bik)
          
          bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
      
          pdf = yy(i,j,k)*bij(i,jj) - dexp(bij(i,jj))
     +        - (bij(i,jj) - bmean1 - rho*bij(i,jj-1))**2
     +        / (2.d0*bsigma)
                
      else
                
          bmean1 = sigmar*(rij(i,jj) - rho*rij(i,jj-1))
     +           + (1.d0 - rho)*(zb + sigmat*taui(i) 
     +                           + sigmai*pp*bik)

          bmean2 = sigmar*(rij(i,jj+1) - rho*rij(i,jj))
     +           + (1.d0 - rho)*(zb + sigmat*taui(i) 
     +                           + sigmai*pp*bik)
          
          bsigma = (1.d0 - rho**2)*(sigmai*pp1)**2
      
          pdf = yy(i,j,k)*bij(i,jj) - dexp(bij(i,jj))
     +        - (bij(i,jj) - bmean1 - rho*bij(i,jj-1))**2
     +        / (2.d0*bsigma)
     +        - (bij(i,jj+1) - bmean2 - rho*bij(i,jj))**2
     +        / (2.d0*bsigma)
                                    
      endif 
                                              
	fbij = -pdf    

	end function
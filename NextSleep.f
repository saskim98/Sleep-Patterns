	program main
      implicit real*8 (a-h, o-z)
	integer, parameter :: ni = 89, nd = 7, nk = 96
	integer, parameter :: nz = 6, nx = 6
	integer, parameter :: nw = 6, nu = 5
	integer, parameter :: ndk = nd*nk, nl = 10
      integer, parameter :: niter = 100000
	
	real*8 yold, bmi, race1, race2, family, parent
	
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
            
    	open(unit =  5, file = 'NextSleepData89.txt') 
    	open(unit =  6, file = 'NextSleep_Initial.txt') 

    	open(unit = 11, file = 'NextSleep_Output1.txt') 
    	open(unit = 12, file = 'NextSleep_Output2.txt') 
    	open(unit = 13, file = 'NextSleep_Output3.txt') 
    	open(unit = 14, file = 'NextSleep_Output4.txt') 
    	open(unit = 15, file = 'NextSleep_Output5.txt') 
      
      iseed = 9999999

      do i = 1, ni
        do j = 1, nd
            do k = 1, nk
                yy(i,j,k) = -999.d0
            enddo
        enddo
      enddo

	do ii = 1, 59808
		
	  read(5,*) i, j, k, 
     +            yold, bmi, race1, race2, family, parent
  
        bmi = (bmi - 0.16854d0)/0.37647d0
        race1 = (race1 - 0.25843d0)/0.44025d0
        race2 = (race2 - 0.37079d0)/0.48575d0
        family = (family - 5.04494d0)/1.56605d0
        parent = (parent - 3.70787d0)/1.77868d0
              		  
	  yy(i,j,k) = yold
	  tik(i,k) = dfloat(k)/dfloat(nk)
        
        xx(i,1) = 1.d0 
        xx(i,2) = bmi 
        xx(i,3) = race1
        xx(i,4) = race2
        xx(i,5) = family
        xx(i,6) = parent

        zz(i,1) = 1.d0 
        zz(i,2) = bmi 
        zz(i,3) = race1
        zz(i,4) = race2
        zz(i,5) = family
        zz(i,6) = parent

        ww(i,1) = 1.d0
        ww(i,2) = bmi 
        ww(i,3) = race1
        ww(i,4) = race2
        ww(i,5) = family
        ww(i,6) = parent

        uu(i,1) = bmi 
        uu(i,2) = race1
        uu(i,3) = race2
        uu(i,4) = family
        uu(i,5) = parent
        	  	  	  	         
      enddo  
      
      sigma_zmu   = 100.d0
      sigma_xbeta = 100.d0
      sigma_wbeta = 100.d0
      sigma_ubeta = 100.d0
      
      a1 = 1.0d0 ; b1 = 0.1d0
      a2 = 1.0d0 ; b2 = 0.1d0
      a3 = 1.0d0 ; b3 = 0.1d0
      a4 = 1.0d0 ; b4 = 0.1d0
      a5 = 1.0d0 ; b5 = 1.0d0      
      a6 = 1.0d0 ; b6 = 1.0d0      

      api = dconst('PI')

      read(6,*) zmu
      read(6,*) xbeta
      read(6,*) wbeta 
      read(6,*) ubeta
      read(6,*) sigmar
      read(6,*) sigmat
      read(6,*) rho 
      read(6,*) sigmaz
      read(6,*) sigmal
      read(6,*) sigmaa
      read(6,*) aalpha
      read(6,*) delta
      read(6,*) vv
      
      vv(nl) = 1.d0
      
      sumd = 0.d0
      do l1 = 1, nl
                   
          sumv = dlog(vv(l1))                
          do l2 = 1, l1-1
              sumv = sumv + dlog(1.d0 - vv(l2))
          enddo
          
          delta(l1) = dexp(sumv)
          sumd = sumd + delta(l1)
          
      enddo      
      do l = 1, nl
          delta(l) = delta(l)/sumd
      enddo                                                 
                                           
      do i = 1, ni
        taui(i) = 0.d0
        zeta(i) = 0.d0
        lambda(i) = 0.d0
        do j = 1, nd        
            do k = 1, nk            
                jj = (j - 1)*nk + k                
                bij(i,jj) = 0.d0
                rij(i,jj) = 0.d0
            enddo
        enddo
      enddo

      call rnset(iseed)            
      
      call gibbs_bij_mle(iseed)          
      call gibbs_taui(iseed)     
      call gibbs_rij(iseed)     
      call gibbs_lambda(iseed)          
      call gibbs_zeta(iseed)        
      do ir = 1, 10
          call gibbs_bij(iseed)          
          call gibbs_taui(iseed)     
          call gibbs_rij(iseed)     
          call gibbs_lambda(iseed)          
          call gibbs_zeta(iseed)        
      enddo
                                          
      do ir = 1, niter
              
          call gibbs(iseed)
        
          write(11,1) ir, zmu, xbeta, wbeta, ubeta
          write(12,1) ir, sigmar, sigmat, rho, 
     +                    sigmaz, sigmal, sigmaa
          write(13,2) ir, aalpha
          write(14,2) ir, delta
          write(15,2) ir, vv
                                   
      enddo    	         
                            
      call rnget(iseed)     
            
    1 format(i6, 1000f20.10)      
    2 format(i6, 1000f30.20)      
    3 format(2i6, 1000f20.10)      
    4 format(4i6, 1000f20.10)      
    5 format(1000f30.10)      

      end program
                
      include 'NextSleep_Gibbs.f'
      include 'optim1.f'

      

! -----------------------------------------------------------------
!  ckthrm.f
!  Thermochemistry Library
!  Version: 2.0
!  Last Update: Nov 3, 2019
!
!  Programmer: Americo Barbosa da Cunha Junior
!              americo.cunhajr@gmail.com
! -----------------------------------------------------------------
!  Copyright (c) 2010-2019, Americo Barbosa da Cunha Junior
!  All rights reserved.
! -----------------------------------------------------------------
!  This is the implementation file for for CKTHRM module,
!  computational library with routines to do the interface with
!  Chemkin-II code.
! -----------------------------------------------------------------




!------------------------------------------------------------
!     ck_wrk_len
!
!     This routine computes work vectors dimensions
!
!     Parameters:
!     leniwk - integer   work vector length
!     lenrwk - double    work vector length
!     lencwk - character work vector length
!------------------------------------------------------------

      subroutine ck_wrk_len(leniwk, lenrwk, lencwk)
      implicit none
      integer :: iflag = 0, lout = 6, link = 25
      integer :: leniwk, lenrwk, lencwk
      
!     opening link file
      open(unit=link, file='chem.bin', form='unformatted', status='old')
     
!     computing work vectors minimum dimensions
      call cklen(link,lout,leniwk,lenrwk,lencwk,iflag)
      if (iflag .gt. 0) then
            stop
      endif
      
!     closing link file
      close(link)
      
      return
      end subroutine
!------------------------------------------------------------





!------------------------------------------------------------
!     ck_wrk_init
!
!     This subroutine starts CHEMKIN-II work vectors
!
!     Parameters:
!     leniwk - integer   work vector length
!     lenrwk - double    work vector length
!     lencwk - character work vector length
!     iwrk   - integer   work vector
!     rwrk   - double    work vector
!     cwrk   - character work vector
!     mm     - number of elements
!     kk     - number of species
!     ii     - number of reactions
!------------------------------------------------------------

      subroutine ck_wrk_init(leniwk, lenrwk, lencwk, iwrk, rwrk, cwrk)
      implicit none
      integer             :: lin = 5, lout = 6, link = 25
      integer             :: iflag = 0
      integer             :: leniwk, lenrwk, lencwk
      integer             :: mm, kk, ii
      integer             :: iwrk(leniwk)
      double precision    :: rwrk(lenrwk)
      character(len = 16) :: cwrk(lencwk)
     
!     opening link file
      open(unit=link, file='chem.bin', form='unformatted', status='old')
     
!     initializing work vectors
      call ckinit(leniwk,lenrwk,lencwk,link,lout,iwrk,rwrk,cwrk,iflag)
      if (iflag .gt. 0) then
            stop
      endif
     
!     closing link file
      close(link)
     
      return
      end subroutine
!------------------------------------------------------------




!------------------------------------------------------------
!     ck_wrk_indx
!
!     This subroutine starts CHEMKIN-II and sets chemical
!     kinetics mechanism parameters.
!
!     Parameters:
!     iwrk   - integer   work vector
!     rwrk   - double    work vector
!     mm     - number of elements
!     kk     - number of species
!     ii     - number of reactions
!------------------------------------------------------------

      subroutine ck_wrk_indx(iwrk, rwrk, mm, kk, ii)
      implicit none
      integer             :: nfit = 0
      integer             :: mm, kk, ii
      integer             :: iwrk(*)
      double precision    :: rwrk(*)
     
!     reaction mechanism size
      call ckindx(iwrk, rwrk, mm, kk, ii, nfit)
     
      return
      end subroutine
!------------------------------------------------------------




!------------------------------------------------------------
!     ck_composition
!
!     This subroutine sets composition vector with values.
!
!     Parameters:
!     kk     - number of species
!     iwrk   - integer   work vector
!     rwrk   - double    work vector
!     cwrk   - character work vector
!     phi    - composition vector
!------------------------------------------------------------

      subroutine ck_composition(kk, iwrk, rwrk, cwrk, phi)
      implicit none
      integer             :: lin = 5, lout = 6
      integer             :: i, ilen, knum, nval
      integer             :: iwrk(*), kk
      double precision    :: rwrk(*), xx(kk), phi(kk+2)
      double precision    :: T, p, RU, RUC, patm, val
      character(len = 16) :: cwrk(*), ksym(kk)
      character(len = 80) :: line = ' '
      logical             :: ierr, kerr = .false.
     
!     character strings of species names
      call cksyms(cwrk, lout, ksym, ierr)
      if (ierr) then
            kerr = .true.
      endif
     
!     initializing molar fraction vector
      do i = 1, KK
        xx(i) = 0.0
      end do
      
!     initializing composition
      do i = 1, KK+2
        phi(i) = 0.0
      end do
      
!     input system temperature
      write (lout, '(/A)') ' system temperature (K):'
      read  (lin,    *) T
      write (lout,7105) T
      
!     input pressure
      write (lout, '(/A)') ' system pressure (atm):'
      read  (lin,    *) p
      write (lout,7105) p
     
!     input initial non-zero chemical species and they moles number
      write (lout, '(/A)') ' chemical species and their moles number:'
      write (lout, '(/A)') ' --- type END to finish the input ---'
     
   40 continue
     
      write (lout, '(/A)') ' next specie and it number of moles:'
      read  (lin,  '(A)', end = 45) line
      write (lout, '(1X,A)') line
      
      ilen = index(line, '!')
      if (ilen .eq. 1) then
            go to 40
      endif
      
      ilen = ilen - 1
      if (ilen .le. 0) then
            ilen = len(line)
      endif
      
      if (index(line(:ilen), 'END') .eq. 0) then
         if (line(:ilen) .ne. ' ') then
            call cksnum(line(:ilen),1,lout,ksym,kk,knum,nval,val,ierr)
            if (ierr) then
               write (lout,*) ' error reading moles...'
               kerr = .true.
            else
               xx(knum) = val
            endif
         endif
         go to 40
      endif
      
   45 continue
     
      if (kerr) then
            stop
      endif
     
!     normalizing molar fraction vector
      call ck_normalize(kk, xx(1))
      
!     converting from molar to mass fraction
      call ckxty(xx(1), iwrk, rwrk, phi(3))
     
!     computing system enthalpy in CGS unit
      call ckhbms(T,phi(3),iwrk,rwrk,phi(1))
     
!     computing atmospheric pressure in CGS
      call ckrp(iwrk,rwrk,RU,RUC,patm)

!     converting pressure to CGS units
      phi(2) = p*patm
      
      
!     input/output formats
 7105 format (12E11.3)
     
      return
      end subroutine
!------------------------------------------------------------




!------------------------------------------------------------
!     ck_species_heading
!
!     This subroutine prints species names on the screen.
!
!     Parameters:
!     kk     - number of species
!     cwrk   - character work vector
!------------------------------------------------------------

      subroutine ck_species_heading(kk,cwrk,ksym)
      implicit none
      integer             :: lin = 5, lout = 6
      integer             :: kk, i
      character(len = 16) :: cwrk(*), ksym(kk)
      logical             :: ierr
     
!     character strings of species names
      call cksyms(cwrk, lout, ksym, ierr)
     
!      write (lout, 7125) (ksym(i)(:11), i=1, kk)
      
!     input/output formats
 7125 format (100(A11))
     
      return
      end subroutine
!------------------------------------------------------------




!------------------------------------------------------------
!     ck_normalize
!
!     This subroutine normalize a vector
!
!     Parameters:
!     n - vector dimension
!     v - vector
!------------------------------------------------------------

      subroutine ck_normalize(nn, vv)
      implicit none
      integer          :: i, nn
      double precision :: vv(nn), tot = 0.d0
     
!     normalizing the mole fractions
      do i = 1, nn
         tot = tot + vv(i)
      end do
      do i = 1, nn
         vv(i) = vv(i) / tot
      end do
      
      return
      end subroutine
!------------------------------------------------------------




!------------------------------------------------------------
!     ck_h2t_newt
!
!     This function returns the mixture temperature
!     in Kelvin given it specific enthalpy and mass
!     fraction of species. (using Newton method)
!
!     Parameters:
!     kk    - number of species
!     phi   - composition vector
!     Tmax  - temperature upper bound
!     Tmin  - temperature lower bound
!     tol   - error tolerance
!     iwrk  - integer work vector
!     rwrk  - double  work vector
!     Tsafe - temperature
!------------------------------------------------------------

      subroutine ck_h2t_newt(kk, phi, tol, iwrk, rwrk, Tguest, Tsafe)
      implicit none
      integer          :: it, kk, iwrk(*), maxit = 100
      double precision :: rwrk(*), phi(kk+2)
      double precision :: tol, Tguest, Tsafe
      double precision :: h0, xn, x0, fx0, dfx0, dx0
      
!     mean enthalpy of the mixture in mass units
      h0 = phi(1)
      
!     initializing Newton parameters
      x0  = Tguest
      dx0 = 1.0
      
      do while ( it .lt. maxit .and. abs(dx0) .gt. tol )
!       computing the function at x0
        call ckhbms(x0, phi(3), iwrk, rwrk, fx0)
        fx0 = fx0 - h0
!       computing the function derivative at x0
        call ckcpbs(x0, phi(3), iwrk, rwrk, dfx0)
!       treatment to avoid division by zero
        if ( abs(dfx0) .le. tol ) then
            dfx0 = 1.0
        end if
!       computing the Newton quotient
        dx0 = fx0/dfx0
!       computing the root aproximation xn
        xn = x0 - dx0
!       updating iteration counter
        it = it + 1
!       updating the initial guest for the root
        x0 = xn
      end do
      
      if ( abs(dx0) .le. tol ) then
        Tsafe = xn
      else
        Tsafe = 0.0
!        print *,'ck_h2t_newt: maximum # of iterations'
!        stop
      end if
      
      return
      end subroutine
!------------------------------------------------------------




!------------------------------------------------------------
!     ck_h2t_bisec
!
!     This function returns the mixture temperature
!     in Kelvin given it specific enthalpy and mass
!     fraction of species. (using Bisection method)
!
!     Parameters:
!     kk    - number of species
!     phi   - composition vector
!     Tmax  - temperature upper bound
!     Tmin  - temperature lower bound
!     tol   - error tolerance
!     iwrk  - integer work vector
!     rwrk  - double  work vector
!     Tsafe - temperature
!------------------------------------------------------------

      subroutine ck_h2t_bisec(kk,phi,Tmax,Tmin,tol,iwrk,rwrk,Tsafe)
      implicit none
      integer          :: kk, iwrk(*), it, maxit = 30
      double precision :: rwrk(*), phi(kk+2)
      double precision :: Tmin, Tmax, tol, Tsafe
      double precision :: h0, fmin, fmax
      double precision :: xmin, xmax, xmid, fxmid
      
      xmin = Tmin
      xmax = Tmax
      
!     mean enthalpy of the mixture in mass units
      h0 = phi(1)
!     mean enthalpy of the mixture at Tmin in mass units
      call ckhbms(xmin, phi(3), iwrk, rwrk, fmin)
      fmin = fmin - h0
!     mean enthalpy of the mixture at Tmax in mass units
      call ckhbms(xmax, phi(3), iwrk, rwrk, fmax)
      fmax = fmax - h0
      
!     checking if is there a root in the interval [Tmim,Tmax]
      if( fmin*fmax .gt. tol ) then
          Tsafe = 0.0
!         print *, 'ck_h2t_bisec: root must be bracketed in [Tmim,Tmax]'
!         stop
          return
      end if
      
!     setting bisection parameters
      it    = 1
      
      do while ( it .lt. maxit )
!       computing the inital guest x0 for the root
        xmid = xmin + 0.5*(xmax - xmin)
!       compute the function at x0
        call ckhbms(xmid, phi(3), iwrk, rwrk, fxmid)
        fxmid = fxmid - h0
!       checking the convergence
        if ( abs(fxmid) .le. tol ) then
            Tsafe = xmid
            return
        end if
!       setting the interval boundary
        if ( fmin*fxmid .ge. tol ) then
            xmin = xmid
        else
            xmax = xmid
        end if
!       updating iteration counter
        it = it + 1
      end do
      
      Tsafe = xmid
      
      return
      end subroutine
!------------------------------------------------------------




!------------------------------------------------------------
!     ck_h2t_safe
!
!     This function returns the mixture temperature
!     in Kelvin given it specific enthalpy and mass
!     fraction of species. (using Newton method)
!
!     Parameters:
!     kk    - number of species
!     phi   - composition vector
!     Tmax  - temperature upper bound
!     Tmin  - temperature lower bound
!     tol   - error tolerance
!     iwrk  - integer work vector
!     rwrk  - double  work vector
!------------------------------------------------------------

!      subroutine ck_h2t_safe(kk, phi, Tmax, Tmin, tol, iwrk, rwrk, Tsafe)
!      implicit none
!      integer          :: i, kk, iwrk(*), maxit = 100
!      double precision :: rwrk(*), phi(kk+2)
!      double precision :: x1, x2, tol, Tmin, Tmax, Tsafe, h0
!      double precision :: dfun, dx, dxold, fun, funh, funl, temp, xh, xl
!      
!      x1 = Tmin
!      x2 = Tmax
!     
!     mean enthalpy of the mixture in mass units     
!      h0 = phi(1)
!     mean enthalpy of mixture for low temp. in mass units
!      call ckhbms(x1, phi(3), iwrk, rwrk, funl)
!      funl = funl - h0
!     mean enthalpy of mixture for high temp. in mass units
!      call ckhbms(x2, phi(3), iwrk, rwrk, funh)
!      funh = funh - h0
!      
!      
!     check if is there a root in the interval [Tmim,Tmax]
!      if( funl*funh .gt. tol ) then
!         print *, 'root must be bracketed in [Tmim,Tmax]'
!         stop
!      end if
!     
!      if( abs(funl) .lt. tol ) then
!        Tsafe = x1
!        return
!      end if
!      
!      if( abs(funh) .lt. tol ) then
!        Tsafe = x2
!        return
!      end if
!     
!      if( funl .lt. tol ) then
!        xl = x1
!        xh = x2
!      else
!        xh = x1
!        xl = x2
!      end if
!      
!      Tsafe = 0.5*(x1 + x2)
!      dxold = abs(x2 - x1)
!      dx    = dxold
!      
!      call ckhbms(Tsafe, phi(3), iwrk, rwrk, fun)
!      fun = fun - h0
!      call ckcpbs(Tsafe, phi(3), iwrk, rwrk, dfun)
!      
!      do i = 1, maxit
!        if( ((Tsafe-xh)*dfun-fun)*((Tsafe-xl)*dfun-fun) .gt. 0.0
!        2                  .or. abs(2.0*fun) .gt. abs(dxold*dfun) ) then
!          dxold = dx
!          dx    = 0.5*(xh - xl)
!          Tsafe = xl + dx
!          
!          if(xl .eq. Tsafe) then
!            return
!          end if
!        else
!          dxold = dx
!          dx    = fun/dfun
!          temp  = Tsafe
!          Tsafe = Tsafe - dx
!          
!          if(temp .eq. Tsafe) then
!            return
!          end if
!        end if
!        
!        if( abs(dx) .lt. tol ) then
!            return
!        end if
!        
!        call ckhbms(Tsafe, phi(3), iwrk, rwrk, fun)
!        fun = fun - h0
!        call ckcpbs(Tsafe, phi(3), iwrk, rwrk, dfun)
!        
!        if( fun .lt. tol ) then
!          xl = Tsafe
!        else
!          xh = Tsafe
!        end if
!      end do
!      
!      print *,'ck_h2t_safe exceeding maximum iterations'
!      stop
!      
!      return
!      end subroutine
!------------------------------------------------------------

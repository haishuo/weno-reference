c
c program: weno.f
c
c purpose: to solve the two dimensional Euler
c          equations of compressible gas dynamics in a rectangle
c          using WENO method.
c
c reference: Efficient implementation of weighted ENO schemes,
c            by Guang-Shan Jiang and Chi-Wang Shu,
c            in Journal of Computational Physics, v126 (1996), pp.202-228.
c
c The code is set up to run a periodic vortex problem.
c However, it can be easily changed to run other problems
c
c this code runs optimally on Cray.  The second()
c and the ismax functions are special for CRAY (should
c comment and uncomment certain lines around these functions before
c running on CRAY).
c
c the result is in the file fort.8,  in the
c order of   x, y, density, pressure, x-velocity, y-velocity
c and in the file fort.9 for the 1D cut at x=xhalf,  in the order of
c y, density, pressure, x-velocity, y-velocity

      program WENO_LF
      include 'comm.inc'
*****678****************************************************************
* Name:      main.f
* Function:  drive routine
*            System to solve: u_t + f(u)_x + g(u)_y = 0
*   or
*            u_t = RHS = -f(u)_x - g(u)_y
*****678****************************************************************

* readin parameters

      call setup

* initialization

      call init

* 'io' is the main controller which rotates between 0 and "mt-1"
* corresponding to the different stage in Runge-Kutta schemes

      io = 0
      istop = 0
      nt = 0

      if(ntot.eq.0) goto 1001

*****678*************** begin time evolution ***************************

1000  continue

* impose boundary condition

call bc(io)

      if( io.eq.0 .and. istop.eq.1) goto 1001
c       if(io.eq.0 .and. second().gt.550.) goto 1001

* compute -f(u)_x

call fx(io)

* compute -g(u)_y

call gy(io)

* compute time step size "dt"

      if( io.eq.0 ) then

c      tcpu=second()
      if(nt/20*20.eq.nt) then
      write(6,*) nt,'t=',tnum,'accumulated cpu=',tcpu
      write(6,*) 'cpu per 20 time steps',tcpu-t00
      t00=tcpu
      end if
call cflc(aam)
      dt=cfl/aam

c         dt = cfl / em / ( cdx + cdy )
      if( ( tnum + dt ) .ge. tend ) then
      dt = tend - tnum
      istop = 1
      endif
      tnum = tnum + dt
      nt = nt + 1
      if( nt.ge.ntot ) istop = 1
      endif

* Runge-Kutta scheme for advancing in time

call rk(io)

      io = mod( io+1, mt )

      goto 1000

*****678***************  end  time evolution ***************************

1001  continue

      write(*,*) 'nt = ', nt, '  tnum = ', tnum

* output

      call save

      stop
      end


      subroutine setup
      include 'comm.inc'
*****678****************************************************************
* Name:      setup.f
* Function:  set up the # of grid points and read in some parameters
*****678****************************************************************

      write(*,*) ' order in time ? (3,4) '
      read(*,*) mt

c       write(*,*) ' # of equations ? (1,2,3,4) '
c       read(*,*) mn
      mn=4

      write(*,*) ' # of points in x-dir & y-dir ? '
      read(*,*) nx,ny
      nxm = nx + md
      nym = ny + md

      write(*,*) ' cfl ? ( e.g. 0.6 )'
      read(*,*) cfl

      write(*,*) ' maximum # of time steps ? '
      read(*,*) ntot

      write(*,*) ' '
      write(*,*) ' WENO-LF-4 or WENO-LF-5 with RK-3/RK-4 '
      write(*,*) ' '
      write(*,*) 'order in time: ', mt
      write(*,*) '# of equations: ', mn
      write(*,*) '# of points in x-dir: ', nx
      write(*,*) '# of points in y-dir: ', ny
      write(*,*) 'CFL = ', cfl
      write(*,*) 'maximum # of time steps ', ntot

      gamma = 1.4
      gm1 = gamma - 1.0

      epweno = 1.e-6

      return
      end


      subroutine init
      include 'comm.inc'
*****678****************************************************************
* Name:      init_per.f
* Function:  set up the grid and initial condition ( u(x,y) at t= 0 )
*****678****************************************************************

      tend = 0.5
      write(6,*) 'the final time'
      read(5,*) tend

      write(6,*) 'enter problem number'
      write(6,*) '  1 = steady vortex'
      write(6,*) '  2 = horizontally moving vortex'
      write(6,*) '  3 = diagonally moving vortex'
      write(6,*) '  4 = 1D Sod shock tube'
      write(6,*) '  5 = 2D Riemann problem' 
      write(6,*) '  6 = blast wave problem'
      read(5,*) nprob

      xleft = 0.
      xright = 10.0
      yleft = 0.
      yright = 10.0

      dx = ( xright - xleft ) / nx
cdx = 1./dx
      dy = ( yright - yleft ) / ny
cdy = 1./dy

      do i = -md, nxm
      x(i) = xleft + i * dx
      enddo

      do j = -md, nym
      y(j) = yleft + j * dy
      enddo

      pi = atan(1.0) * 4.0

c rkap is the strength of the vortex

      rkap=5.
coefk=rkap/(2.*pi*exp(-0.5))
      gga=1./(gamma-1.)
      ggb=gamma/(gamma-1.)
      ggc=1.0/ggb
      ggd=-1.0/gamma

      if(nprob.eq.1) then

      dri=1.
      uri=0.
      vri=0.
      pri=1.

      else if(nprob.eq.2) then

      dri=1.
      uri=1.
      vri=0.
      pri=1.

      else if(nprob.eq.3) then

      dri=1.
      uri=1.
      vri=1.
      pri=1.


      else if(nprob.eq.4) then
c     1D Sod shock tube (extended to 2D)
      do j=-md,nym
      do i=-md,nxm
      if(x(i) .lt. 5.0) then
c        Left state: high pressure/density
         uc(i,j,1,0) = 1.0
         uc(i,j,2,0) = 0.0
         uc(i,j,3,0) = 0.0
         uc(i,j,4,0) = 1.0/(gamma-1.0)
      else
c        Right state: low pressure/density
         uc(i,j,1,0) = 0.125
         uc(i,j,2,0) = 0.0
         uc(i,j,3,0) = 0.0
         uc(i,j,4,0) = 0.1/(gamma-1.0)
      endif
      enddo
      enddo

      else if(nprob.eq.5) then
c     2D Riemann problem (4-quadrant configuration)
      do j=-md,nym
      do i=-md,nxm
      if(x(i).lt.5.0 .and. y(j).lt.5.0) then
c        Quadrant I: high pressure
         uc(i,j,1,0) = 1.0
         uc(i,j,2,0) = 0.0
         uc(i,j,3,0) = 0.0
         uc(i,j,4,0) = 1.0/(gamma-1.0)
      else if(x(i).ge.5.0 .and. y(j).lt.5.0) then
c        Quadrant II: medium pressure + x-velocity
         uc(i,j,1,0) = 0.5197
         uc(i,j,2,0) = 0.5197 * (-0.7259)
         uc(i,j,3,0) = 0.0
         uc(i,j,4,0) = 0.4/(gamma-1.0) + 0.5*0.5197*(-0.7259)**2
      else if(x(i).lt.5.0 .and. y(j).ge.5.0) then
c        Quadrant III: medium pressure + y-velocity
         uc(i,j,1,0) = 0.5197
         uc(i,j,2,0) = 0.0
         uc(i,j,3,0) = 0.5197 * (-0.7259)
         uc(i,j,4,0) = 0.4/(gamma-1.0) + 0.5*0.5197*(-0.7259)**2
      else
c        Quadrant IV: low pressure
         uc(i,j,1,0) = 0.1379
         uc(i,j,2,0) = 0.0
         uc(i,j,3,0) = 0.0
         uc(i,j,4,0) = 0.029/(gamma-1.0)
      endif
      enddo
      enddo

      else if(nprob.eq.6) then
c     Blast wave problem (high pressure center)
      do j=-md,nym
      do i=-md,nxm
      xr = x(i) - 5.0
      yr = y(j) - 5.0
      rr = sqrt(xr*xr + yr*yr)
      if(rr .lt. 1.0) then
c        High pressure blast region
         uc(i,j,1,0) = 1.0
         uc(i,j,2,0) = 0.0
         uc(i,j,3,0) = 0.0
         uc(i,j,4,0) = 1000.0/(gamma-1.0)
      else
c        Low pressure ambient
         uc(i,j,1,0) = 1.0
         uc(i,j,2,0) = 0.0
         uc(i,j,3,0) = 0.0
         uc(i,j,4,0) = 0.01/(gamma-1.0)
      endif
      enddo
      enddo
      end if

      do 1 j=-md,nym
      do 1 i=-md,nxm
cij=dri
      rij=uri
      sij=vri
      pij=pri

      xr=x(i)-5.
      yr=y(j)-5.
      rr2=xr**2+yr**2
      rr=sqrt(rr2)

      eee=exp(-0.5*rr2)
      pold=pij
      told=pij/cij
      entr=pij/cij**gamma
      rij=rij-coefk*eee*yr
      sij=sij+coefk*eee*xr
      tper=-0.5*coefk**2*eee*eee*ggc
      trin=told+tper

      pij=trin**ggb/entr**gga
      cij=pij/trin

      uc(i,j,1,0)=cij
      uc(i,j,2,0)=cij*rij
      uc(i,j,3,0)=cij*sij
      uc(i,j,4,0)=pij/(gamma-1.)+0.5*cij*(rij**2+sij**2)

1     continue

c     Override with shock tube initial conditions for problems 4,5,6
      if(nprob.eq.4) then
c       1D Sod shock tube - reset all values
        do j=-md,nym
        do i=-md,nxm
        if(x(i) .lt. 5.0) then
c         Left state: high pressure/density
          uc(i,j,1,0) = 1.0
          uc(i,j,2,0) = 0.0
          uc(i,j,3,0) = 0.0
          uc(i,j,4,0) = 1.0/(gamma-1.0)
        else
c         Right state: low pressure/density
          uc(i,j,1,0) = 0.125
          uc(i,j,2,0) = 0.0
          uc(i,j,3,0) = 0.0
          uc(i,j,4,0) = 0.1/(gamma-1.0)
        endif
        enddo
        enddo

      else if(nprob.eq.5) then
c       2D Riemann problem - reset all values
        do j=-md,nym
        do i=-md,nxm
        if(x(i).lt.5.0 .and. y(j).lt.5.0) then
          uc(i,j,1,0) = 1.0
          uc(i,j,2,0) = 0.0
          uc(i,j,3,0) = 0.0
          uc(i,j,4,0) = 1.0/(gamma-1.0)
        else if(x(i).ge.5.0 .and. y(j).lt.5.0) then
          uc(i,j,1,0) = 0.5197
          uc(i,j,2,0) = 0.5197 * (-0.7259)
          uc(i,j,3,0) = 0.0
          uc(i,j,4,0) = 0.4/(gamma-1.0) + 0.5*0.5197*(-0.7259)**2
        else if(x(i).lt.5.0 .and. y(j).ge.5.0) then
          uc(i,j,1,0) = 0.5197
          uc(i,j,2,0) = 0.0
          uc(i,j,3,0) = 0.5197 * (-0.7259)
          uc(i,j,4,0) = 0.4/(gamma-1.0) + 0.5*0.5197*(-0.7259)**2
        else
          uc(i,j,1,0) = 0.1379
          uc(i,j,2,0) = 0.0
          uc(i,j,3,0) = 0.0
          uc(i,j,4,0) = 0.029/(gamma-1.0)
        endif
        enddo
        enddo

      else if(nprob.eq.6) then
c       Blast wave problem - reset all values
        do j=-md,nym
        do i=-md,nxm
        xr = x(i) - 5.0
        yr = y(j) - 5.0
        rr = sqrt(xr*xr + yr*yr)
        if(rr .lt. 1.0) then
          uc(i,j,1,0) = 1.0
          uc(i,j,2,0) = 0.0
          uc(i,j,3,0) = 0.0
          uc(i,j,4,0) = 1000.0/(gamma-1.0)
        else
          uc(i,j,1,0) = 1.0
          uc(i,j,2,0) = 0.0
          uc(i,j,3,0) = 0.0
          uc(i,j,4,0) = 0.01/(gamma-1.0)
        endif
        enddo
        enddo

      endif
c     End of shock initialization override


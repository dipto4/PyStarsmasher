      subroutine GetXYZFromFile(filename, ntotout ,pX,pY,pZ)
      
      include 'starsmasher.h'
      
      character*255 filename
      integer i
      real*8 energy, fract
      real*8 divv(nmax)
      real*8 displacex, displacey, displacez
      integer ndisplace
      real*8 erad
      real*8 pX(nmax), pY(nmax), pZ(nmax)
      integer ntotout

      open(12,file=filename,form="unformatted")
      
      read(12) n,nnopt,hco,hfloor,sep0,tf,dtout,nout,nit,t,
     $  nav,alpha,beta,tjumpahead,ngr,nrelax,trelax,dt,omega2,
     $  ncooling,erad,ndisplace,displacex,displacey,displacez
    
      print*, n

      if (ncooling .eq. 0) then
        do i = 1, n
          read(12) x(i),y(i),z(i),am(i),hp(i),rho(i),
     $      vx(i),vy(i),vz(i),vxdot(i),vydot(i),vzdot(i),
     $      u(i),udot(i), gx(i),gy(i),gz(i),
     $      grpot(i),meanmolecular(i),cc(i),divv(i)
        end do
      else
        do i = 1, n
          read(12) x(i),y(i),z(i),am(i),hp(i),rho(i),
     $      vx(i),vy(i),vz(i),vxdot(i),vydot(i),vzdot(i),
     $      u(i),udot(i), gx(i),gy(i),gz(i),
     $      grpot(i),meanmolecular(i),cc(i),divv(i),
     $      ueq(i),tthermal(i)
        end do
      end if
      
      pX = x
      pY = y
      pZ = z
      ntotout = n  
      end subroutine


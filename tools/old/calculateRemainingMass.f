      SUBROUTINE CalculateBoundMass(filename, boundmass)
      
      include 'starsmasher.h'
      
      character*255 filename
      real*8 boundmass
      integer i
      real*8 energy, fract
      real*8 mass0, mass1
      real*8 px0,px1,py0,py1,pz0,pz1
      real*8 vx0,vx1,vy0,vy1,vz0,vz1
      real*8 divv(nmax)
      real*8 displacex, displacey, displacez
      integer ndisplace
      real*8 erad
c     component 0 is unbound
c     component 1 is bound
      integer icomp(nmax) !keeps track of what component each particle is in
      integer numchanged, count
    

      open(12,file=filename,form="unformatted")
      
      read(12) n,nnopt,hco,hfloor,sep0,tf,dtout,nout,nit,t,
     $  nav,alpha,beta,tjumpahead,ngr,nrelax,trelax,dt,omega2,
     $  ncooling,erad,ndisplace,displacex,displacey,displacez


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


c     first guess at what component particles are in
      do i=1,n
         icomp(i)=1
      enddo
c Work out how to separate what mass is bound to what object 
      count=0
      numchanged=n
      do while (numchanged.ne.0)

         mass0=0d0; mass1=0d0
         px0=0d0; py0=0d0; pz0=0d0
         px1=0d0; py1=0d0; pz1=0d0
         numchanged=0
         
         do i=1,n
          if (icomp(i).EQ.1) then
            mass1=am(i)+mass1
            px1 = px1 + am(i)*vx(i)
            py1 = py1 + am(i)*vy(i)
            pz1 = pz1 + am(i)*vz(i)
          else 
            mass0=am(i)+mass0

          endif 
         enddo
         
         vx1=px1/mass1
         vy1=py1/mass1
         vz1=pz1/mass1
         
         do i=1,n
            energy = u(i)*am(i)+0.5d0*am(i)*((vx(i)-vx1)**2+
     $           (vy(i)-vy1)**2+(vz(i)-vz1)**2)+grpot(i)*am(i)
            
            if(energy .ge. .0) then
c     particle is unbound
               if (icomp(i).ne.0) then
                  icomp(i)=0
                  numchanged=numchanged+1
               endif
            else
c     particle is bound
               if (icomp(i).ne.1) then
                  icomp(i)=1
                  numchanged=numchanged+1
               endif
            endif
            
            
         enddo

         count=count+1



      enddo 
      boundmass = mass1

      RETURN
      END

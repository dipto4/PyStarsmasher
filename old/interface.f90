! Contains subroutine to set all values from Python for now

! Subroutine responsible for setting all of the values that the user puts in python

subroutine PythonSetValues( Pndisplace, Pdisplacex, Pdisplacey, Pdisplacez, Psemimajoraxis&
        &,Pbimpact,Pe0,Pvinf2, Ptf, Pdtout, Pn, Pgflag, Pnnopt, Pnav, Palpha, Pbeta, Pngr,& 
        &Phco, Phfloor,Pnrelax,Ptrelax,Psep0,Pequalmass,Ptreloff,Ptresplintmuoff,Pnitpot,&
        &Ptscanon,Psepfinal,Pnintvar,Pngravprocs,Pqthreads,Pmbh,Prunit,Pmunit,Pcn1,Pcn2,Pcn3,&
        &Pcn4,Pcn5,Pcn6,Pcn7,Pcomputeexclusivemode,Pomega_spin,Pppn,Pneos,Pnselfgravity,Pgam,&
        &Preat, Pstarmass,Pstarradius, Pncooling, Pnkernel, Pteq, Ptjumpahead, Pstartfile1,Pstartfile2,&
        &Peosfile,Popacityfile,Pprofilefile, Pthrowaway,Pstellarevolutioncodetype,Psimulationtype,Pdirname,Pnpoly)

    include 'starsmasher.h'
    common/displace/displacex,displacey,displacez,ndisplace
    common/simtype/simulationtype

    integer ndisplace
    real*8 displacex, displacey, displacez

    ! python variables start here
    integer Pndisplace
    real*8 Pdisplacex,Pdisplacey,Pdisplacez

    real*8 Psemimajoraxis, Pe0, Pbimpact, Ptrelax
    real*8 Pvinf2, Ptf, Pdtout, Pequalmass  
    integer Pn, Pgflag, Pnnopt, Pnav, Pngr, Pnrelax
    real*8 Palpha, Pbeta
    real*8 Pcn1, Pcn2, Pcn3, Pcn4, Pcn5, Pcn6, Pcn7
    real*8 Phco, Phfloor, Ptscanon, Psepfinal, Psep0, Ptreloff 
    real*8 Ptresplintmuoff
    integer Pnitpot, Pnintvar, Pngravprocs, Pqthreads
    real*8 Pmbh
    real*8 Prunit, Pmunit
    integer Pcomputeexclusivemode, Pppn
    real*8 Pomega_spin
    integer Pneos, Pnselfgravity, Pncooling, Pnkernel
    real*8 Pgam, Preat, Pteq, Ptjumpahead
    real*8 Pstarmass, Pstarradius
    character*255 Pstartfile1, Pstartfile2, Peosfile, Popacityfile, Pprofilefile
    logical Pthrowaway
    integer Pstellarevolutioncodetype
    character*3 simulationtype,Psimulationtype
    character*32 dirname,Pdirname
    real*8 npoly, Pnpoly
    common/pydirname/dirname
    common/polytype/npoly

    ndisplace = Pndisplace
    displacex = Pdisplacex
    displacey = Pdisplacey
    displacez = Pdisplacez
    semimajoraxis = Psemimajoraxis
    e0 = Pe0
    bimpact = Pbimpact
    trelax = Ptrelax
    vinf2 = Pvinf2
    tf = Ptf
    dtout = Pdtout
    equalmass = Pequalmass
    n = Pn
    gflag = Pgflag
    nnopt = Pnnopt
    nav = Pnav
    ngr = Pngr
    nrelax = Pnrelax
    alpha = Palpha
    beta = Pbeta
    cn1 = Pcn1
    cn2 = Pcn2
    cn3 = Pcn3
    cn4 = Pcn4
    cn5 = Pcn5
    cn6 = Pcn6
    cn7 = Pcn7
    hco = Phco
    hfloor = Phfloor
    tscanon = Ptscanon
    sepfinal = Psepfinal
    sep0 = Psep0
    treloff = Ptreloff
    tresplintmuoff = Ptresplintmuoff
    nitpot = Pnitpot
    nintvar = Pnintvar
    ngravprocs = Pngravprocs
    qthreads = Pqthreads
    mbh = Pmbh
    runit = Prunit
    munit = Pmunit
    computeexclusivemode=Pcomputeexclusivemode
    ppn = Pppn
    omega_spin = Pomega_spin
    neos = Pneos
    nselfgravity = Pnselfgravity
    ncooling = Pncooling
    nkernel = Pnkernel
    gam = Pgam
    reat = Preat
    teq = Pteq
    tjumpahead = Ptjumpahead
    starmass = Pstarmass
    starradius = Pstarradius
    startfile1 = Pstartfile1
    startfile2 = Pstartfile2
    eosfile = Peosfile
    opacityfile = Popacityfile
    profilefile = Pprofilefile
    throwaway = Pthrowaway
    simulationtype = Psimulationtype
    
    dirname = Pdirname
    npoly = Pnpoly
    ! for testing purposes only
    !nprocs = 16    
end subroutine

subroutine PythonInitializeDouble(Px1, Py1, Pz1, Px2, Py2, Pz2, Pvx1, Pvy1,&
    &Pvz1, Pvx2, Pvy2, Pvz2, Pm1, Pm2)
    
    include 'starsmasher.h'

    common/doubleinitialize/x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2,m1,m2
    
    real*8 Px1, Py1, Pz1, Px2, Py2, Pz2, Pvx1, Pvy1,Pvz1, Pvx2, Pvy2, Pvz2,&
    &Pm1,Pm2

    real*8 x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2,m1,m2
    
    x1 = Px1
    y1 = Py1
    z1 = Pz1

    x2 = Px2
    y2 = Py2
    z2 = Pz2

    vx1 = Pvx1
    vy1 = Pvy1
    vz1 = Pvz1

    vx2 = Pvx2
    vy2 = Pvy2
    vz2 = Pvz2

    m1 = Pm1
    m2 = Pm2


end subroutine

subroutine TestPythonDouble()
    common/doubleinitialize/x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2,m1,m2
    real*8 x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2,m1,m2

    print *, x1,y1,z1,x2,y2,z2
    print *, vx1,vy1,vz1,vx2,vy2,vz2
    print *, m1, m2
end subroutine

subroutine testPrint
    print *," Hello from inside starsmasher library!"

end subroutine

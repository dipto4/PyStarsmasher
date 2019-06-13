!New PythonSetValues subroutine utilies MPI_BCAST and structs 
!instead of setting the values individually

subroutine PythonSetValues(parent)

    use mpi
    include 'starsmasher.h'
    common/displace/displacex,displacey,displacez,ndisplace
    common/simtype/simulationtype
    
    !parent mpi process
    integer parent


    integer ndisplace
    real*8 displacex, displacey, displacez

    ! python variables start here
    character*255 Pstartfile1, Pstartfile2, Peosfile, Popacityfile, Pprofilefile
    character*3 simulationtype,Psimulationtype
    character*32 dirname,Pdirname
    real*8 npoly
    common/pydirname/dirname
    common/polytype/npoly
    
    type , BIND(C) :: inputs
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
        logical Pthrowaway
        integer Pstellarevolutioncodetype
        real*8 Pnpoly

    end type inputs
    ! declaring mpi variables required for MPI_Type_Create_struct
    ! ci = current inputs
    type(inputs) :: ci

    integer :: blocklen(57), typ(57)
    integer :: input_type, ierr
    integer(KIND=MPI_ADDRESS_KIND) :: disp(57), base, lb, extent
    
    integer i

    !geting MPI addreses

    call MPI_GET_ADDRESS(ci%Pndisplace,disp(1),ierr)
    call MPI_GET_ADDRESS(ci%Pdisplacex,disp(2),ierr)
    call MPI_GET_ADDRESS(ci%Pdisplacey,disp(3),ierr)
    call MPI_GET_ADDRESS(ci%Pdisplacez,disp(4),ierr)
    call MPI_GET_ADDRESS(ci%Psemimajoraxis,disp(5),ierr)
    call MPI_GET_ADDRESS(ci%Pe0,disp(6),ierr)
    call MPI_GET_ADDRESS(ci%Pbimpact,disp(7),ierr)
    call MPI_GET_ADDRESS(ci%Ptrelax,disp(8),ierr)
    call MPI_GET_ADDRESS(ci%Pvinf2,disp(9),ierr)
    call MPI_GET_ADDRESS(ci%Ptf,disp(10),ierr)
    call MPI_GET_ADDRESS(ci%Pdtout,disp(11),ierr)
    call MPI_GET_ADDRESS(ci%Pequalmass,disp(12),ierr)
    call MPI_GET_ADDRESS(ci%Pn,disp(13),ierr)
    call MPI_GET_ADDRESS(ci%Pgflag,disp(14),ierr)
    call MPI_GET_ADDRESS(ci%Pnnopt,disp(15),ierr)
    call MPI_GET_ADDRESS(ci%Pnav,disp(16),ierr)
    call MPI_GET_ADDRESS(ci%Pngr,disp(17),ierr)
    call MPI_GET_ADDRESS(ci%Pnrelax,disp(18),ierr)
    call MPI_GET_ADDRESS(ci%Palpha,disp(19),ierr)
    call MPI_GET_ADDRESS(ci%Pbeta,disp(20),ierr)
    call MPI_GET_ADDRESS(ci%Pcn1,disp(21),ierr)
    call MPI_GET_ADDRESS(ci%Pcn2,disp(22),ierr)
    call MPI_GET_ADDRESS(ci%Pcn3,disp(23),ierr)
    call MPI_GET_ADDRESS(ci%Pcn4,disp(24),ierr)
    call MPI_GET_ADDRESS(ci%Pcn5,disp(25),ierr)
    call MPI_GET_ADDRESS(ci%Pch6,disp(26),ierr)
    call MPI_GET_ADDRESS(ci%Pcn7,disp(27),ierr)

    call MPI_GET_ADDRESS(ci%Phco,disp(28),ierr)
    call MPI_GET_ADDRESS(ci%Phfloor,disp(29),ierr)
    call MPI_GET_ADDRESS(ci%Ptscanon,disp(30),ierr)
    call MPI_GET_ADDRESS(ci%Psepfinal,disp(31),ierr)
    call MPI_GET_ADDRESS(ci%Psep0,disp(32),ierr)
    call MPI_GET_ADDRESS(ci%Ptreloff,disp(33),ierr)
    call MPI_GET_ADDRESS(ci%Ptresplintmuoff,disp(34),ierr)
    call MPI_GET_ADDRESS(ci%Pnitpot,disp(35),ierr)
    call MPI_GET_ADDRESS(ci%Pnintvar,disp(36),ierr)
    call MPI_GET_ADDRESS(ci%Pngravprocs,disp(37),ierr)
    call MPI_GET_ADDRESS(ci%Pqthreads,disp(38),ierr)
    call MPI_GET_ADDRESS(ci%Pmbh,disp(39),ierr)
    call MPI_GET_ADDRESS(ci%Prunit,disp(40),ierr)
    call MPI_GET_ADDRESS(ci%Pmunit,disp(41),ierr)
    call MPI_GET_ADDRESS(ci%Pcomputeexclusivemode,disp(42),ierr)
    call MPI_GET_ADDRESS(ci%Ppn,disp(43),ierr)
    call MPI_GET_ADDRESS(ci%Pomega_spin,disp(44),ierr)
    call MPI_GET_ADDRESS(ci%Pneos,disp(45),ierr)
    call MPI_GET_ADDRESS(ci%Pnselfgravity,disp(46),ierr)
    call MPI_GET_ADDRESS(ci%Pncooling,disp(47),ierr)
    call MPI_GET_ADDRESS(ci%Pnkernel,disp(48),ierr)
    call MPI_GET_ADDRESS(ci%Pgam,disp(49),ierr)
    call MPI_GET_ADDRESS(ci%Preat,disp(50),ierr)
    call MPI_GET_ADDRESS(ci%Pteq,disp(51),ierr)
    call MPI_GET_ADDRESS(ci%Ptjumpahead,disp(52),ierr)
    call MPI_GET_ADDRESS(ci%Pstarmass,disp(53),ierr)
    call MPI_GET_ADDRESS(ci%Pstarradius,disp(54),ierr)
    call MPI_GET_ADDRESS(ci%Pthrowaway,disp(55),ierr)
    call MPI_GET_ADDRESS(ci%Pstellarevolutioncodetype,disp(56),ierr)


    call MPI_GET_ADDRESS(ci%Pnpoly,disp(57),ierr)

    base = disp(1)
    
    do i=1, 57
       disp(i) = disp(i) - base 
        blocklen(i) = 1
    end do
    
    typ(1) = MPI_INTEGER
    typ(2) = MPI_DOUBLE
    typ(3) = MPI_DOUBLE
    typ(4) = MPI_DOUBLE
    typ(5) = MPI_DOUBLE
    typ(6) = MPI_DOUBLE
    typ(7) = MPI_DOUBLE
    typ(8) = MPI_DOUBLE
    typ(9) = MPI_DOUBLE
    typ(10) = MPI_DOUBLE
    typ(11) = MPI_DOUBLE
    typ(12) = MPI_DOUBLE
    typ(13) = MPI_INTEGER
    typ(14) = MPI_INTEGER
    typ(15) = MPI_INTEGER
    typ(16) = MPI_INTEGER
    typ(17) = MPI_INTEGER
    typ(18) = MPI_INTEGER
    typ(19) = MPI_DOUBLE
    typ(20) = MPI_DOUBLE
    typ(21) = MPI_DOUBLE
    typ(22) = MPI_DOUBLE
    typ(23) = MPI_DOUBLE
    typ(24) = MPI_DOUBLE
    typ(25) = MPI_DOUBLE
    typ(26) = MPI_DOUBLE
    typ(27) = MPI_DOUBLE
    typ(28) = MPI_DOUBLE
    typ(29) = MPI_DOUBLE
    typ(30) = MPI_DOUBLE
    typ(31) = MPI_DOUBLE
    typ(32) = MPI_DOUBLE
    typ(33) = MPI_DOUBLE
    typ(34) = MPI_DOUBLE
    typ(35) = MPI_INTEGER
    typ(36) = MPI_INTEGER
    typ(37) = MPI_INTEGER
    typ(38) = MPI_INTEGER
    typ(39) = MPI_DOUBLE
    typ(40) = MPI_DOUBLE
    typ(41) = MPI_DOUBLE
    typ(42) = MPI_INTEGER
    typ(43) = MPI_INTEGER
    typ(44) = MPI_DOUBLE
    typ(45) = MPI_INTEGER
    typ(46) = MPI_INTEGER
    typ(47) = MPI_INTEGER
    typ(48) = MPI_INTEGER
    typ(49) = MPI_DOUBLE
    typ(50) = MPI_DOUBLE
    typ(51) = MPI_DOUBLE
    typ(52) = MPI_DOUBLE
    typ(53) = MPI_DOUBLE
    typ(54) = MPI_DOUBLE
    typ(55) = MPI_LOGICAL
    typ(56) = MPI_INTEGER
    typ(57) = MPI_DOUBLE

    

    call MPI_TYPE_CREATE_STRUCT(57,blocklen,disp,typ,input_type,ierr)
    call MPI_TYPE_COMMIT(input_type,ierr)


    call MPI_BCAST(ci,1,input_type,parent,ierr)
    
    call MPI_BCAST(Pstartfile1,255,MPI_CHARACTER,0,parent,ierr)
    call MPI_BCAST(Pstartfile2,255,MPI_CHARACTER,0,parent,ierr)
    call MPI_BCAST(Peosfile,255,MPI_CHARACTER,0,parent,ierr)
    call MPI_BCAST(Popacityfile,255,MPI_CHARACTER,0,parent,ierr)
    call MPI_BCAST(Pprofilefile,255,MPI_CHARACTER,0,parent,ierr)
    
    call MPI_BCAST(Psimulationtype,3,MPI_CHARACTER,0,parent,ierr)
    call MPI_BCAST(Pdirname,32,MPI_CHARACTER,0,parent,ierr)


    
    ndisplace = ci%Pndisplace
    displacex = ci%Pdisplacex
    displacey = ci%Pdisplacey
    displacez = ci%Pdisplacez
    semimajoraxis = ci%Psemimajoraxis
    e0 = ci%Pe0
    bimpact = ci%Pbimpact
    trelax = ci%Ptrelax
    vinf2 = ci%Pvinf2
    tf = ci%Ptf
    dtout = ci%Pdtout
    equalmass = ci%Pequalmass
    n = ci%Pn
    gflag = ci%Pgflag
    nnopt = ci%Pnnopt
    nav = ci%Pnav
    ngr = ci%Pngr
    nrelax = ci%Pnrelax
    alpha = ci%Palpha
    beta = ci%Pbeta
    cn1 = ci%Pcn1
    cn2 = ci%Pcn2
    cn3 = ci%Pcn3
    cn4 = ci%Pcn4
    cn5 = ci%Pcn5
    cn6 = ci%Pcn6
    cn7 = ci%Pcn7
    hco = ci%Phco
    hfloor = ci%Phfloor
    tscanon = ci%Ptscanon
    sepfinal = ci%Psepfinal
    sep0 = ci%Psep0
    treloff = ci%Ptreloff
    tresplintmuoff = ci%Ptresplintmuoff
    nitpot = ci%Pnitpot
    nintvar = ci%Pnintvar
    ngravprocs = ci%Pngravprocs
    qthreads = ci%Pqthreads
    mbh = ci%Pmbh
    runit = ci%Prunit
    munit = ci%Pmunit
    computeexclusivemode=ci%Pcomputeexclusivemode
    ppn = ci%Pppn
    omega_spin = ci%Pomega_spin
    neos = ci%Pneos
    nselfgravity = ci%Pnselfgravity
    ncooling = ci%Pncooling
    nkernel = ci%Pnkernel
    gam = ci%Pgam
    reat = ci%Preat
    teq = ci%Pteq
    tjumpahead = ci%Ptjumpahead
    starmass = ci%Pstarmass
    starradius = Pstarradius
    startfile1 = Pstartfile1
    startfile2 = Pstartfile2
    eosfile = Peosfile
    opacityfile = Popacityfile
    profilefile = Pprofilefile
    throwaway = ci%Pthrowaway
    simulationtype = Psimulationtype
    
    dirname = Pdirname
    npoly = ci%Pnpoly
    
end subroutine

subroutine PythonInitializeDouble(parent)
    use mpi
    include 'starsmasher.h'

    common/doubleinitialize/x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2,m1,m2
    
    type, BIND(C) :: orbital_params
        real*8 x,y,z
        real*8 vx,vy,vz
        real*8 m
    end type orbital_params

    real*8 x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2,m1,m2

    type(orbital_params) :: star1
    type(orbital_params) :: star2
    
    integer :: blocklen(7), typ(7)
    integer :: params_type, ierr
    integer(KIND=MPI_ADDRESS_KIND) :: disp1(7), disp2(7),base1, base2, lb, extent
    

    integer i


    ! star 1

    call MPI_GET_ADDRESS(star1%x,disp1(1),ierr)
    call MPI_GET_ADDRESS(star1%y,disp1(2),ierr)
    call MPI_GET_ADDRESS(star1%z,disp1(3),ierr)
    call MPI_GET_ADDRESS(star1%vx,disp1(4),ierr)
    call MPI_GET_ADDRESS(star1%vy,disp1(5),ierr)
    call MPI_GET_ADDRESS(star1%vz,disp1(6),ierr)
    call MPI_GET_ADDRESS(star1%m,disp1(7),ierr)
    base1 = disp1(1)

    ! star 2

    call MPI_GET_ADDRESS(star2%x,disp2(1),ierr)
    call MPI_GET_ADDRESS(star2%y,disp2(2),ierr)
    call MPI_GET_ADDRESS(star2%z,disp2(3),ierr)
    call MPI_GET_ADDRESS(star2%vx,disp2(4),ierr)
    call MPI_GET_ADDRESS(star2%vy,disp2(5),ierr)
    call MPI_GET_ADDRESS(star2%vz,disp2(6),ierr)
    call MPI_GET_ADDRESS(star2%m,disp2(7),ierr)
    base2 = disp2(1)

    do i = 1,7
        disp1(i) = disp1(i) - base1
        disp2(i) = disp2(i) - base2
        blocklen(i) = 1
        typ(i) = MPI_DOUBLE
    end do

    call MPI_TYPE_CREATE_STRUCT(7,blocklen,disp,typ,params_type,ierr)
    call MPI_TYPE_COMMIT(params_type,ierr)

    call MPI_BCAST(star1,1,params_type,parent,ierr)
    call MPI_BCAST(star2,1,params_type,parent,ierr)




    x1 = star1%x
    y1 = star1%y
    z1 = star1%z

    x2 = sta2%x
    y2 = star2%y
    z2 = star2%z

    vx1 = star1%vx
    vy1 = star1%vy
    vz1 = star1%vz

    vx2 = star2%vx
    vy2 = star2%vy
    vz2 = star2%vz

    m1 = star1%m
    m2 = star2%m

end subroutine

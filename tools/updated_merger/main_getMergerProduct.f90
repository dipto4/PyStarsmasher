subroutine getMergerProduct(filename,ms1,ms2,ms3)
    implicit none
    character*32 filename
    real*8 ms1,ms2,ms3
    call READIT(filename)
    
    call trajectories(ms1,ms2,ms3)
    

end subroutine

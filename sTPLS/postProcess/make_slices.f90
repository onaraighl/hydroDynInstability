!! This code creates cross section of 3d data in botht the x-z plane and
!! the x-y plane.
!! Rory Sullivan

program mainprogram
    implicit none
    
    integer                :: file_interval,file_start,file_end
    integer                :: iteration,l,m,n,j_cent,k_cent
    double precision       :: dx,dy,dz
    character(60)          :: filename
    character(60)          :: filename_out
    
    ! Need to specify how many .dat files will be read in.  Program will
    ! take no inputs so have to ensure that the number of points in the
    ! x,y and z directions are correctly specified as in the code.
    
    parameter (file_start=1000,file_end=400000,file_interval=1000,l=241,m=141,n=121)  
    ! (maxl = 241, maxm = 141, maxn = 121)
    ! Coordinate data sets
    
    double precision :: u3(0:l-2,0:m-2,0:n-2),v3(0:l-2,0:m-2,0:n-2), &
            w3(0:l-2,0:m-2,0:n-2),pres3(0:l-2,0:m-2,0:n-2),visc3(0:l-2,0:m-2,0:n-2)
    
    ! Parameter definiton
    
    dx=8.d0/dble(l-1)
    dy=4.d0/dble(m-1)
    dz=1.d0/dble(n-1)
    
    ! Slice positions.
    
    j_cent=(m-1)/2
    k_cent=(n-1)/2
    
    
    ! Loop over all data files in the directory.
    
    do iteration=file_start,file_end,file_interval
    
        ! Get the filename.
        write(*,*) 'Accessing data file:'
        write(filename,'(A,I9,A)')'sp3dchannel_',iteration,'.dat'
        write(*,*) filename
    
        ! Load data.
        call load_3d_data(u3,v3,w3,pres3,visc3,l,m,n,filename)

        ! Make the slices and write to file.
        write(*,*) 'Generating output file:'
        write(filename_out,'(A,I9,A)')'2dslice_xz_',iteration,'.dat'
        write(*,*) filename_out
        call channel_output_2dslice_xz(u3,v3,w3,pres3,visc3,l,m,n,j_cent,dx,dz,filename_out)
        
        write(filename_out,'(A,I9,A)')'2dslice_xy_',iteration,'.dat'
        write(*,*) filename_out
        call channel_output_2dslice_xy(u3,v3,w3,pres3,visc3,l,m,n,k_cent,dx,dy,filename_out)
        
        write(*,*) 'iteration #',iteration,' done'
        
        
    end do

end program mainprogram

! ****************************************************************************************
! Subroutine to read in 3d data.
! ****************************************************************************************
    
subroutine load_3d_data(u3,v3,w3,pres3,visc3,l,m,n,filename)
    implicit none
    
    integer,intent(in)              :: l,m,n
    character(60),intent(in)        :: filename
    
    double precision ,intent(out)   :: u3(0:l-2,0:m-2,0:n-2)
    double precision ,intent(out)   :: v3(0:l-2,0:m-2,0:n-2)
    double precision ,intent(out)   :: w3(0:l-2,0:m-2,0:n-2)
    double precision ,intent(out)   :: pres3(0:l-2,0:m-2,0:n-2)
    double precision ,intent(out)   :: visc3(0:l-2,0:m-2,0:n-2)
    
    
    integer                         :: i,j,k
    double precision                :: x,y,z
    

    open(unit=9,file=filename,status='old')
    
    do k=0,n-2
      do j=0,m-2
        do i=0,l-2
            read(9,*)x,y,z,u3(i,j,k),v3(i,j,k),w3(i,j,k),pres3(i,j,k),visc3(i,j,k)
            
        end do
      end do
    end do
    
    close(9)
    
    return
    
end subroutine load_3d_data

! ****************************************************************************************
! Subroutine to generate x-z slices and write data to file.
! ****************************************************************************************

subroutine channel_output_2dslice_xz(u3,v3,w3,pres3,visc3,l,m,n,j_cent,dx,dz,filename_out)
    implicit none
    
    integer,intent(in)          :: l,m,n,j_cent
    double precision,intent(in) :: dx,dz
    double precision,intent(in) :: u3(0:l-2,0:m-2,0:n-2)
    double precision,intent(in) :: v3(0:l-2,0:m-2,0:n-2)
    double precision,intent(in) :: w3(0:l-2,0:m-2,0:n-2)
    double precision,intent(in) :: pres3(0:l-2,0:m-2,0:n-2)
    double precision,intent(in) :: visc3(0:l-2,0:m-2,0:n-2)
    character(60),intent(in)    :: filename_out
    
    integer                     :: i,k
    double precision            :: x,z
    double precision            :: u,v,w,pres,visc

    open(unit=9,file=filename_out,status='unknown')
    write(9,*)'VARIABLES="X","Z","U","v","W","Pressure","Viscosity"'
    write(9,*)'ZONE T="Floor", I=',l-1,' K=',n-1,' F=POINT'
    
    
    do k=0,n-2
      z=dz/2.d0+dz*k
      
        do i=0,l-2
          x=dx/2.d0+dx*i
          
          u=u3(i,j_cent,k)
          v=v3(i,j_cent,k)
          w=w3(i,j_cent,k)
          pres=pres3(i,j_cent,k)
          visc=visc3(i,j_cent,k)
          
          write(9,*)x,z,u,v,w,pres,visc
          
        end do
      end do
   
    close(9)
    
    return
    
end subroutine channel_output_2dslice_xz

! ****************************************************************************************
! Subroutine to generate x-y slices and write data to file.
! ****************************************************************************************

subroutine channel_output_2dslice_xy(u3,v3,w3,pres3,visc3,l,m,n,k_cent,dx,dy,filename_out)
    implicit none
    
    integer,intent(in)          :: l,m,n,k_cent
    double precision,intent(in) :: dx,dy
    double precision,intent(in) :: u3(0:l-2,0:m-2,0:n-2)
    double precision,intent(in) :: v3(0:l-2,0:m-2,0:n-2)
    double precision,intent(in) :: w3(0:l-2,0:m-2,0:n-2)
    double precision,intent(in) :: pres3(0:l-2,0:m-2,0:n-2)
    double precision,intent(in) :: visc3(0:l-2,0:m-2,0:n-2)
    character(60),intent(in)    :: filename_out
    
    integer                     :: i,j
    double precision            :: x,y
    double precision            :: u,v,w,pres,visc

    open(unit=9,file=filename_out,status='unknown')
    write(9,*)'VARIABLES="X","Y","U","V","W","Pressure","Viscosity"'
    write(9,*)'ZONE T="Floor", I=',l-1,' J=',m-1,' F=POINT'
    
    
    do j=0,m-2
      y=dy/2.d0+dy*j
      
        do i=0,l-2
          x=dx/2.d0+dx*i
          
          u=u3(i,j,k_cent)
          v=v3(i,j,k_cent)
          w=w3(i,j,k_cent)
          pres=pres3(i,j,k_cent)
          visc=visc3(i,j,k_cent)
          
          write(9,*)x,y,u,v,w,pres,visc
          
        end do
      end do
   
    close(9)
    
    return
    
end subroutine channel_output_2dslice_xy

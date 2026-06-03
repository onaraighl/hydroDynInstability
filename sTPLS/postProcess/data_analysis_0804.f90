program mainprogram
        implicit none

    ! Need to specify how many .dat files will be read in.
    ! Program will take no inputs so have to ensure that the number of points in the
    ! x,y and z directions are correctly specified as in the code.

    integer :: file_start,file_end,file_interval,num_files,iteration,k,l,m,n, &
            file_i,t_equil_i,time_i,ind_2,z_mid
    double precision :: time_val,wnorm_val,Lx,Ly,dt,dx,dy,dz,time_equil,final_time
   ! parameter (file_start=1000,file_end=233000,file_interval=1000,num_files=233,l=81,m=41,n=41) !!! 
    parameter (file_start=1000,file_end=92000,file_interval=1000,num_files=399,l=241,m=141,n=121) !!! 
    character(60) filename

    ! Data sets involved in spatially averaging process.
    double precision :: u_av_xy(0:n-2),v_av_xy(0:n-2),w_av_xy(0:n-2)
    double precision :: uu_av_xy(0:n-2),vv_av_xy(0:n-2),ww_av_xy(0:n-2),uv_av_xy(0:n-2), &
                uw_av_xy(0:n-2),vw_av_xy(0:n-2)
    
    ! Data sets to store spatial averages.
    double precision :: u_av_t(0:(num_files-1),0:n-2),v_av_t(0:(num_files-1),0:n-2),w_av_t(0:(num_files-1),0:n-2), &
        uu_av_t(0:(num_files-1),0:n-2),vv_av_t(0:(num_files-1),0:n-2),ww_av_t(0:(num_files-1),0:n-2), &
        uv_av_t(0:(num_files-1),0:n-2),uw_av_t(0:(num_files-1),0:n-2),vw_av_t(0:(num_files-1),0:n-2)
    
    ! Data sets to store time and spatial averages.
    double precision :: u_stav(0:n-2),v_stav(0:n-2),w_stav(0:n-2),uu_stav(0:n-2),vv_stav(0:n-2), &
            ww_stav(0:n-2),uv_stav(0:n-2),uw_stav(0:n-2),vw_stav(0:n-2)

    double precision :: u_stav_temp,v_stav_temp,w_stav_temp,uu_stav_temp,vv_stav_temp,ww_stav_temp,uv_stav_temp, &
                uw_stav_temp,vw_stav_temp

    ! Data sets to store the RMS values and Reynolds stress.
    double precision :: u_rms(0:n-2),v_rms(0:n-2),w_rms(0:n-2),rey_stress(0:n-2)
    
    ! Data sets for spatially averaged velocities at z=0.5
    double precision :: u_av_zmid_t(0:(num_files-1)),v_av_zmid_t(0:(num_files-1)),w_av_zmid_t(0:(num_files-1))
    double precision :: u_zmid_stav,v_zmid_stav,w_zmid_stav

    ! Coordinate data sets
    double precision :: x_2d(0:l-2,0:n-2),z_2d(0:l-2,0:n-2),u_2d(0:l-2,0:n-2),w_2d(0:l-2,0:n-2)
    ! double precision :: z(0:(n-2)),u3(0:l-2,0:m-2,0:n-2), &
    !         w3(0:l-2,0:m-2,0:n-1),pres(0:l-2,0:m-2,0:n-2),visc(0:l,0:m-2,0:n-2)
	
	double precision :: z(0:(n-2))
	double precision, allocatable, dimension(:,:,:) :: u3,w3,pres,visc

    ! z plus and additional analysis variables.
    double precision :: zplus(0:(n-2)),zplus_xy(0:(n-2)),zplus_t(0:(num_files-1),0:n-2),zplus_temp
    double precision :: time(0:(num_files-1))
	
	! allocate arrays
	
	allocate( 	u3(0:l-2,0:m-2,0:n-2), &
				w3(0:l-2,0:m-2,0:n-1),pres(0:l-2,0:m-2,0:n-2),visc(0:l,0:m-2,0:n-2) )

    ! Parameter definiton
    Lx=8.d0
    Ly=4.d0
    dt=1d-4
        dx=8.d0/dble(l-1)      
        dy=4.d0/dble(m-1)
        dz=1.d0/dble(n-1)
    time_equil=file_start*dt
    ! Needs to be set after running the code with ind_2=0
    t_equil_i=10
    final_time=file_end*dt
    ind_2=0
     z_mid=int((n-1)/2)-1

    ! Open a test file to store the useful information.
    if(ind_2==2)then
        open(unit=8,file='averaged_velocities.dat',status='unknown')
    end if


    open(unit=2,file='u_max.dat',status='unknown')
    open(unit=1,file='filenames.dat',status='unknown')

    ! File counter
    file_i=0
    
    ! Loop over all data files in the directory.
    do iteration=file_start,file_end,file_interval

        ! Get the filename.
        write(1,*) 'Accessing data file:'
        write(1,*) iteration
        write(filename,'(A,I9,A)')'sp3dchannel_',iteration,'.dat'
        write(1,*) filename
    
        ! Read data from file.
        ! call access_data_2dslice(wnorm_val,x_2d,z_2d,u_2d,w_2d,l,m,m,dx,dy,dz,filename)
        time_val=iteration*dt
        time(file_i)=time_val

        ! Access spatially averaged data.
        call average_3ddata_xy(u3,w3,z,zplus_xy,u_av_xy,v_av_xy,w_av_xy,uu_av_xy,vv_av_xy,ww_av_xy,uv_av_xy, &
                    uw_av_xy,vw_av_xy,pres,visc,l,m,n,filename)

        ! Store spatial averages.
        u_av_t(file_i,:)=u_av_xy
        v_av_t(file_i,:)=v_av_xy
        w_av_t(file_i,:)=w_av_xy
        uu_av_t(file_i,:)=uu_av_xy
        vv_av_t(file_i,:)=vv_av_xy
        ww_av_t(file_i,:)=ww_av_xy
        uv_av_t(file_i,:)=uv_av_xy
        uw_av_t(file_i,:)=uw_av_xy
        vw_av_t(file_i,:)=vw_av_xy

        ! Store zplus_xy
        zplus_t(file_i,:)=zplus_xy

        ! Obatin the z=0.5 value of the velocity averages.
        u_av_zmid_t(file_i)=u_av_xy(z_mid)
        v_av_zmid_t(file_i)=v_av_xy(z_mid)
        w_av_zmid_t(file_i)=w_av_xy(z_mid)
        
        ! Iterate counter
        file_i=file_i+1
    
    end do

    ! Time averaging
    if(ind_2.eq.2)then

        ! Now we need to compute the time averages (from equilibrium time onwards).
        ! Also compute zplus
        do k=0,n-2
            u_stav_temp=0.d0
            v_stav_temp=0.d0
            w_stav_temp=0.d0
            uu_stav_temp=0.d0
            vv_stav_temp=0.d0
            ww_stav_temp=0.d0
            uv_stav_temp=0.d0
            uw_stav_temp=0.d0
            vw_stav_temp=0.d0
            zplus_temp=0.d0
            do time_i=t_equil_i,(num_files-1)
                u_stav_temp=u_stav_temp+u_av_t(time_i,k)
                v_stav_temp=v_stav_temp+v_av_t(time_i,k)
                w_stav_temp=w_stav_temp+w_av_t(time_i,k)
                uu_stav_temp=uu_stav_temp+uu_av_t(time_i,k)
                vv_stav_temp=vv_stav_temp+vv_av_t(time_i,k)
                ww_stav_temp=ww_stav_temp+ww_av_t(time_i,k)
                uv_stav_temp=uv_stav_temp+uv_av_t(time_i,k)
                uw_stav_temp=uw_stav_temp+uw_av_t(time_i,k)
                vw_stav_temp=vw_stav_temp+vw_av_t(time_i,k)
                zplus_temp=zplus_temp+zplus_t(time_i,k)
            end do
            u_stav(k)=u_stav_temp
            v_stav(k)=v_stav_temp
            w_stav(k)=w_stav_temp
            uu_stav(k)=uu_stav_temp
            vv_stav(k)=vv_stav_temp
            ww_stav(k)=ww_stav_temp
            uv_stav(k)=uv_stav_temp
            uw_stav(k)=uw_stav_temp
            vw_stav(k)=vw_stav_temp
            zplus(k)=zplus_temp
        end do
        u_stav=u_stav/(num_files-t_equil_i)
        v_stav=v_stav/(num_files-t_equil_i)
        w_stav=w_stav/(num_files-t_equil_i)
        uu_stav=uu_stav/(num_files-t_equil_i)
        vv_stav=vv_stav/(num_files-t_equil_i)
        ww_stav=ww_stav/(num_files-t_equil_i)
        uv_stav=uv_stav/(num_files-t_equil_i)
        uw_stav=uw_stav/(num_files-t_equil_i)
        vw_stav=vw_stav/(num_files-t_equil_i)
        zplus=zplus/(num_files-t_equil_i)

        ! Now we have the space time averages of the velocities.
        ! We compute the root mean square (RMS) of u,v,and w.
        u_rms=sqrt(uu_stav-(u_stav**2))
        v_rms=sqrt(vv_stav-(v_stav**2))
        w_rms=sqrt(ww_stav-(w_stav**2)) 
    
        ! Compute the reynolds stress using the u,w averages.
        rey_stress=-(uw_stav-(u_stav)*(w_stav))

        ! Now we need to compute the time averages of the two point correlation functions.
        ! Get time average of the z=0.5 averaged velocities

        u_zmid_stav=0.d0
        v_zmid_stav=0.d0
        w_zmid_stav=0.d0

        do time_i=t_equil_i,(num_files-1)
          u_zmid_stav=u_zmid_stav+u_av_zmid_t(time_i)
          v_zmid_stav=v_zmid_stav+v_av_zmid_t(time_i)
          w_zmid_stav=w_zmid_stav+w_av_zmid_t(time_i)       
        end do

        u_zmid_stav=u_zmid_stav/(num_files-t_equil_i)   
        v_zmid_stav=v_zmid_stav/(num_files-t_equil_i)
        w_zmid_stav=w_zmid_stav/(num_files-t_equil_i)

        !Now write data to file.
        write(8,*)'VARIABLES="Z ","Z_PLUS","U_STAV ","U_RMS ","V_RMS ","W_RMS ","T_REY"'
            write(8,*)'ZONE T="Floor", I=',l-1,' J=',m-1,' K=',n-1,' F=POINT'
        do k=0,(n-2)
            write(8,*) z(k),zplus(k),u_stav(k),u_rms(k),v_rms(k),w_rms(k),rey_stress(k)
        end do  

        close(8)

    end if

    write(2,*)'VARIABaug="Time ","U_AVG_ZMID "'
    write(2,*)'ZONE T="Floor", I=',l-1,' J=',m-1,' K=',n-1,' F=POINT'
    do k=0,(num_files-1)
        write(2,*) time(k),u_av_zmid_t(k)
    end do
    close(2)

    close(1)
    

end program mainprogram

! ****************************************************************************************

    subroutine access_data_2dslice(wnorm_val,x_2d,z_2d,u_2d,w_2d,l,m,n,dx,dy,dz,filename)
    implicit none

    integer :: l,m,n,i,k
    double precision :: dx,dy,dz
    double precision :: wnorm_val,temp_val
    double precision :: x_2d(0:l-2,0:n-2),z_2d(0:l-2,0:n-2),u_2d(0:l-2,0:n-2),w_2d(0:l-2,0:n-2), &
            w_2d_temp(0:l-2)

    character(60)    filename

    open(unit=9,file=filename,status='old')

    do k=0,n-2
            do i=0,l-2
            read(9,*)x_2d(i,k),z_2d(i,k),u_2d(i,k),w_2d(i,k)
            end do
        enddo

    close(9)

    w_2d_temp=sum(w_2d**2)
    temp_val=sum(w_2d_temp)*dx*dz
    wnorm_val=sqrt(temp_val)

        return
        end subroutine access_data_2dslice

! *****************************************************************************************************

    subroutine average_3ddata_xy(u3,w3,z,zplus_xy,u_av_xy,v_av_xy,w_av_xy,uu_av_xy,vv_av_xy,ww_av_xy,uv_av_xy, &
                    uw_av_xy,vw_av_xy,pres,visc,l,m,n,filename)
    implicit none

    integer :: l,m,n,i,j,k
    double precision :: u3(0:l-2,0:m-2,0:n-2),v3(0:l-2,0:m-2,0:n-2),w3(0:l-2,0:m-2,0:n-2), &
            pres(0:l-2,0:m-2,0:n-2),visc(0:l-2,0:m-2,0:n-2)
    double precision :: x,y,z(0:n-2)

    ! Averaging processes.
    double precision :: u2d_slice(0:l-2,0:m-2),v2d_slice(0:l-2,0:m-2),w2d_slice(0:l-2,0:m-2)
    double precision :: u_av_xy(0:n-2),v_av_xy(0:n-2),w_av_xy(0:n-2)
    double precision :: uu3(0:l-2,0:m-2,0:n-2),vv3(0:l-2,0:m-2,0:n-2),ww3(0:l-2,0:m-2,0:n-2), &
                uv3(0:l-2,0:m-2,0:n-2),uw3(0:l-2,0:m-2,0:n-2),vw3(0:l-2,0:m-2,0:n-2)
    double precision :: uu3_slice(0:l-2,0:m-2),vv3_slice(0:l-2,0:m-2),ww3_slice(0:l-2,0:m-2), &
                uv3_slice(0:l-2,0:m-2),uw3_slice(0:l-2,0:m-2),vw3_slice(0:l-2,0:m-2)
    double precision :: uu_av_xy(0:n-2),vv_av_xy(0:n-2),ww_av_xy(0:n-2),uv_av_xy(0:n-2), &
                uw_av_xy(0:n-2),vw_av_xy(0:n-2)

    ! zplus variables
    double precision :: visc_slice(0:l-2,0:m-2), zplus_xy(0:n-2)
    

    character(60) filename

    open(unit=9,file=filename,status='old')
    
    do k=0,n-2
      do j=0,m-2
        do i=0,l-2
        read(9,*)x,y,z(k),u3(i,j,k),v3(i,j,k),w3(i,j,k),pres(i,j,k),visc(i,j,k)
        end do
      end do
    end do

    close(9)

    ! Take products of the velocities (needed for correlation functions)
    uu3=u3*u3
    vv3=v3*v3
    ww3=w3*w3
    uv3=u3*v3
    uw3=u3*w3
    vw3=v3*w3

    ! Compute space average (in the x and y direction) of the velocities.
    ! Pre allocate zplus_xy
    zplus_xy=0.d0

    do k=0,n-2
        
        ! Average in x and y of the velocities
        u2d_slice=u3(:,:,k)
        v2d_slice=v3(:,:,k)
        w2d_slice=w3(:,:,k)
        u_av_xy(k)=sum(u2d_slice)/((l-1)*(m-1))
        v_av_xy(k)=sum(v2d_slice)/((l-1)*(m-1))
        w_av_xy(k)=sum(w2d_slice)/((l-1)*(m-1))
        
        ! Average in x and y of the squared velocities
        uu3_slice=uu3(:,:,k)
        vv3_slice=vv3(:,:,k)
        ww3_slice=ww3(:,:,k)
        uv3_slice=uv3(:,:,k)
        uw3_slice=uw3(:,:,k)
        vw3_slice=vw3(:,:,k)
        uu_av_xy(k)=sum(uu3_slice)/((l-1)*(m-1))
        vv_av_xy(k)=sum(vv3_slice)/((l-1)*(m-1))
        ww_av_xy(k)=sum(ww3_slice)/((l-1)*(m-1))
        uv_av_xy(k)=sum(uv3_slice)/((l-1)*(m-1))
        uw_av_xy(k)=sum(uw3_slice)/((l-1)*(m-1))
        vw_av_xy(k)=sum(vw3_slice)/((l-1)*(m-1))

        ! Average the viscosity in the x and y direction.
        visc_slice=visc(:,:,k)
        do j=0,m-2
          do i=0,l-2
            zplus_xy(k)=zplus_xy(k) + (z(k)/(visc_slice(i,j)))
          end do
        end do
        zplus_xy(k)=zplus_xy(k)/((l-1)*(m-1))

    end do


    end subroutine average_3ddata_xy

! *****************************************************************************************************

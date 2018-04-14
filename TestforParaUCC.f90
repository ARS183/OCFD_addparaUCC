program TestforParaUCC
implicit none
include 'openNS3d.h'
real(kind=OCFD_REAL_KIND),allocatable :: x(:),u(:),d2f(:),df(:),dfexat(:)
real(kind=OCFD_REAL_KIND),allocatable :: xx(:),yy(:),temp_x(:),temp_y(:)
real(kind=OCFD_REAL_KIND) :: time_start,time_end,ctime
real(kind=OCFD_REAL_KIND):: aa, pi,erri,errmax1,errmax2,rate_globle
!real(kind=OCFD_REAL_KIND) :: slx,hx
integer :: ka,k,i,kid,i_global,ii,iii,npx1,npx2,my_mod1
character*100 filename2
!integer :: npx0, nx_global, npx1, npx2
!real(kind=OCFD_REAL_KIND),allocatable:: i_nn, i_offset

pi=datan(1.d0)*4.d0
aa=0.d0

call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,np_size,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)

if (my_id==0) then
	open(unit=11,file='rate.dat')

	write(*,"(A7,3X,2(A10,3X))")"No.","Err.w.b","Rate"
	write(11,"(A7,3X,2(A10,3X))")"No.","Err.w.b","Rate"
  ! Err.w.b和Err.wo.b分别表示包括、不包括边界点值时的误差。由于边界处的格式与内点格式不同，故均于此列出
  ! 其后的阶数Rate分别与其之前的Err对应，不再赘写后缀
end if

!iii=1
do iii=1,8
	! call read_parameter()
nx_global=8*2**iii+1
slx=1.d0
npx0=np_size
hx=slx/dble(nx_global-1)

!call part3d()
npx=mod(my_id,npx0)
nx=nx_global/npx0
if(npx .lt. mod(nx_global,npx0)) nx=nx+1

do k=0,npx0-1
    ka=min(k,mod(nx_global,npx0))
    i_offset(k)=int(nx_global/npx0)*k+ka+1
    i_nn(k)=nx_global/npx0
    if(k .lt. mod(nx_global,npx0)) i_nn(k)=i_nn(k)+1
enddo
npx1=my_mod1(npx-1,npx0)
npx2=my_mod1(npx+1,npx0)
ID_XM1=npx1    ! -1 proc in x-direction
ID_XP1=npx2 
if(npx .eq. 0) ID_XM1=MPI_PROC_NULL     ! if not periodic, 0 node donot send mesg to npx0-1 node
if(npx .eq. npx0-1) ID_XP1=MPI_PROC_NULL

allocate(x(1-LAP:nx+LAP),u(1-LAP:nx+LAP))
allocate(d2f(1-LAP:nx+LAP),df(1-LAP:nx+LAP),dfexat(1-LAP:nx+LAP))


!call define_grid(x)
if (my_id .eq. 0) then

    allocate(xx(1:nx_global))
    do i=1,nx_global
        xx(i)=dble(i-1)*hx
    enddo
    
    do 100 i=0,npx0-1		
		ka=i
			
		if (ka .eq. 0) then
			do ii=1,nx
					i_global=i_offset(i)-1+ii
					x(ii)=xx(i_global)
			enddo

		else
			allocate(temp_x(1:i_nn(i)))
			do ii=1,i_nn(i)
				i_global=i_offset(i)+ii-1
				temp_x(ii)=xx(i_global)
			enddo
					
	   	    call MPI_SEND(temp_x,i_nn(i),OCFD_DATA_TYPE,  &
			ka,1,MPI_COMM_WORLD,ierr) 

		    deallocate(temp_x)
		endif
100	continue

	deallocate(xx)
	
	else
	    allocate(temp_x(1:i_nn(npx)))

		call MPI_RECV(temp_x,i_nn(npx),OCFD_DATA_TYPE,  &
        0,1,MPI_COMM_WORLD,status,ierr) 
        
		do ii=1,i_nn(npx)
			x(ii)=temp_x(ii)
		enddo	

		deallocate(temp_x)

endif

call check_x1d(x)

call testfunc(x,u,dfexat)

call check_x1d(u)


time_start=MPI_WTIME()
if (npx==0) then
    call OCFD_D2F_BOUND_PADE4(u,d2f,nx,hx,1)
    call OCFD_D2F_SB_PADE4(u,d2f,nx,hx,1)
    call D2F_PADE4(u,d2f,nx,hx)
    call OCFD_DF_BOUND_UCC45(u,df,nx,hx,1)
    call DF_UCC45_P(u,d2f,df,nx,hx,1)
elseif (npx==npx0-1) then
    call OCFD_D2F_BOUND_PADE4(u,d2f,nx,hx,2)
    call OCFD_D2F_SB_PADE4(u,d2f,nx,hx,2)
    call D2F_PADE4(u,d2f,nx,hx)
    call OCFD_DF_SB_UCC45(u,df,nx,hx)
    call OCFD_DF_BOUND_UCC45(u,df,nx,hx,2)
    call DF_UCC45_P(u,d2f,df,nx,hx,2)
else
    call OCFD_D2F_SB_PADE4(u,d2f,nx,hx,0)
    call D2F_PADE4(u,d2f,nx,hx)
    call OCFD_DF_SB_UCC45(u,df,nx,hx)
    call DF_UCC45_P(u,d2f,df,nx,hx,0)
endif

erri = maxval(dabs(df(1:nx)-dfexat(1:nx)))

call MPI_Reduce(erri,errmax1,1,OCFD_DATA_TYPE,MPI_MAX,0,MPI_COMM_WORLD,status,ierr)

       


            !write (filename2,"('file-' I4.4 '.dat')") my_id
            !open(unit=my_id,file=filename2)
            !call display(df,nx,1)
            !close(my_id)

if (my_id==0) then
	if (iii>1) then
		rate_globle=dlog(errmax2/errmax1)/dlog(2.d0)
	end if
	errmax2=errmax1

	if (iii==1) then
		write(*,"(I7,3X,E12.4,3X,A10,3X)")nx_global,errmax1,"/"
		write(11,"(I7,3X,E12.4,3X,A10,3X)")nx_global,errmax1,"/"
	else
		write(*,"(I7,3X,2(E12.4,3X))")nx_global,errmax1,rate_globle
		write(11,"(I7,3X,2(E12.4,3X))")nx_global,errmax1,rate_globle
	end if
end if
deallocate(x,u,d2f,df,dfexat)
enddo

time_end=MPI_WTIME()

call MPI_BARRIER(MPI_COMM_WORLD,ierr)

ctime=time_end-time_start
if (my_id==0) then
        
	write(*,*)"Walltime=",ctime
	write(11,*)"Walltime=",ctime
	close(11)
end if



!deallocate()
call MPI_FINALIZE(ierr)

end program TestforParaUCC



function my_mod1(i,n)
          
    integer my_mod1,i,n
    if(i.lt.0) then
        my_mod1=i+n
    else if (i.gt.n-1) then
        my_mod1=i-n
    else
        my_mod1=i
    endif
end







subroutine testfunc(x,u,dfexat)
    include 'openNS3d.h'
    integer i,k1,npx1,npx2,my_mod1
    real(kind=OCFD_REAL_KIND):: x(1-LAP:nx+LAP),u(1-LAP:nx+LAP),dfexat(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND):: pi
    
    pi=4.d0*datan(1.d0)	
    
    do i=1,nx
        u(i)=dsin(2.d0*pi*x(i))+dcos(2.d0*pi*x(i))
        dfexat(i)=2.d0*pi*( dcos(2.d0*pi*x(i))-dsin(2.d0*pi*x(i)) )
    enddo

end subroutine


!send and recv mesg, to check array in x direction.
subroutine check_x1d(f)
	include 'openNS3d.h'
	integer i,j,k1,npx1,npx2,my_mod1
	real(kind=OCFD_REAL_KIND):: f(1-LAP:nx+LAP),  &
		    tmp_send1(1:LAP),tmp_send2(1:LAP), tmp_recv1(1:LAP),tmp_recv2(1:LAP)
	  
		
		do i=1,LAP
		k1=i
		tmp_send1(k1)=f(i)
		tmp_send2(k1)=f(nx-LAP+i)
		enddo
    
!!

		call MPI_Sendrecv(tmp_send1,LAP,OCFD_DATA_TYPE, ID_XM1, 9000, &
		    tmp_recv2, LAP,  OCFD_DATA_TYPE,ID_XP1, 9000,MPI_COMM_WORLD,Status,ierr)
		call MPI_Sendrecv(tmp_send2,LAP, OCFD_DATA_TYPE,ID_XP1,8000,   &
		    tmp_recv1,LAP, OCFD_DATA_TYPE,ID_XM1,  8000,MPI_COMM_WORLD,Status,ierr)

		if(ID_XM1 .ne. MPI_PROC_NULL) then
		 do i=1,LAP
		   k1=i
		   f(i-LAP)=tmp_recv1(k1)
		 enddo
		endif

		if(ID_XP1 .ne. MPI_PROC_NULL) then
		 do i=1,LAP
		   k1=i
		   f(nx+i)=tmp_recv2(k1)
		 enddo
        endif
        
end subroutine


subroutine display(tar,n,n0)
	include 'openNS3d.h'
    real(kind=OCFD_REAL_KIND)::tar(1-LAP:nx+LAP)
    integer :: j,n,n0
    do j=n0,n
        write(my_id,*) tar(j)
    end do
    !write(*,*)

end subroutine
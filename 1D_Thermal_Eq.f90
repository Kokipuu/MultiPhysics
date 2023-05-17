module constants
   implicit none
   integer::i,nt
   integer, parameter :: imax = 100,ntmax = 2500
   real(8):: T(0:imax),xmax,x,dx,dt,kp
   integer,parameter :: SP = kind(1.0)
   integer,parameter :: DP = selected_real_kind(2*precision(1.0_SP))
end module constants
   
  module print
    use constants
    implicit none
   
   contains
    subroutine print__profile(f)
      real(DP), dimension(0:imax), intent(in) :: f
      integer :: counter = 0
      character(len=*), parameter :: base = "./data/temp.j=middle."
      character(len=4) :: serial_num
      character(len=40) filename
    
      write(serial_num,'(i4.4)') counter
   
      write(filename,"(a,i5.5,a)") base//serial_num
      open(10,file=filename)
      do i = 0 , imax
         write(10,*) i*dx, f(i)
      end do
      close(10)
      counter = counter + 1
      
    end subroutine print__profile
  end module print
   
  program thermal_explicit
    use constants
    use print
   
    implicit none
   
    dt   = 0.004
    kp   = 1
    xmax = 10
    dx   = xmax/imax
    
    open(10,file="initial.txt")
    open(11,file="final.txt")
   
     !initial state
     do i = 0, imax
        x = dx*i
        if (i < imax/2) then
           T(i) = x
         else
           T(i) = xmax - x
        endif
        write(10,*) i,x,T(i)
     enddo
     call print__profile(T(:))
   
     
     !Time step
     do nt = 1, ntmax
        do i = 1, imax-1
           T(i) = T(i) + kp*dt/(dx*dx)*(T(i+1)- 2*T(i) +T(i-1))
        end do
        if (mod(nt,100)==0) call print__profile(T(:))
        write(*,*)nt
     end do
     
     !Final state
     do i = 0, imax
        x = dx*i
        write(11,*) i,x,T(i)
     end do
  end program thermal_explicit
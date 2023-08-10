program control_file
        implicit none

        ! input related variables
        character(len=100) :: buffer, label
        character(len=100) :: arg
        character(len=100) :: filename
        integer :: pos
        integer, parameter :: fh=15
        integer :: ios = 0
        integer :: line = 0
        integer :: i
        ! control file variables
        real :: pi = 0
        integer, dimension(5) :: vector


        i = 0
        do

                call get_command_argument(i,arg)
                if (len_trim(arg) .eq. 0) exit
                print*,i,trim(arg)
                i = i + 1

        enddo

        call get_command_argument(1,filename)
        print*,trim(filename)
        open(fh, file=trim(filename))

        ! ios is negative if an end of record condition is encountered
        ! or if an endfile condition was detected. It is positive if 
        ! an error was detected, ios is zero otherwise

        do while (ios .eq. 0)
          read(fh,'(A)',iostat=ios) buffer
          if (ios .eq. 0) then
                  line = line + 1

                  ! find the first instance of whitespace, split label and data
                  pos = scan(buffer,' ')
                  label = buffer(1:pos)
                  buffer = buffer(pos+1:)

                  select case(label)
                  case ('pi')
                          read(buffer, *, iostat=ios) pi
                          print *, 'Read pi:', pi
                  case ('vector')
                          read(buffer, *, iostat=ios) vector
                          print *, 'Read vector:',vector
                  case default
                          print *, "Skipping invalid label at line", line
                  end select
          end if
        enddo
end program control_file



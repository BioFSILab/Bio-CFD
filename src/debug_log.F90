module biocfd_debug_log
       implicit none
       private

       integer :: log_file_unit

#if CREATE_LOGFILE == 1
       logical, parameter :: write_to_logfile=.true.
#else
       logical, parameter :: write_to_logfile=.false.
#endif
       
       public :: log_file_open, log_file_close, log_file_write
       
       contains

         subroutine log_file_open()
           if(write_to_logfile) open(newunit=log_file_unit, file="log_file.txt")
         end subroutine log_file_open

         subroutine log_file_write(subroutine_name, log_message)
           character(len=*), intent(in) :: subroutine_name, log_message
           if(write_to_logfile) write(log_file_unit,*) subroutine_name," ",log_message           
         end subroutine log_file_write

         subroutine log_file_close()
           if(write_to_logfile) close(log_file_unit)
         end subroutine log_file_close

end module biocfd_debug_log


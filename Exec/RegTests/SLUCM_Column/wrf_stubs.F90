!===============================================================================
! Minimal stand-ins for the two WRF modules that phys/module_sf_urban.F uses.
!
! These exist so the PRISTINE upstream file can be compiled without a WRF tree.
! The upstream source is never edited -- that is the whole point of the
! comparison -- so its `use module_wrf_error` and
! `use module_model_constants, only : piconst` must be satisfiable.
!
! piconst is copied verbatim from WRF share/module_model_constants.F:69,
! including the bare (unsuffixed) literal, so that the promotion flag decides
! its kind exactly as it does upstream. ERF's NoahmpUrbanShimMod writes the same
! digits with an explicit kind; both yield the same double.
!
! wrf_message / wrf_error_fatal mirror NoahmpUrbanShimMod's behaviour. They are
! reached only on a table-parsing failure, never on the trajectory being
! compared, so they cannot influence the numbers.
!===============================================================================

module module_model_constants

   implicit none
   public

   REAL, PARAMETER :: piconst = 3.1415926535897932384626433    ! constant of PI

end module module_model_constants


module module_wrf_error

   implicit none
   public

contains

   subroutine wrf_message(msg)
      character(len=*), intent(in) :: msg
      write(*,'(a)') trim(msg)
   end subroutine wrf_message

   subroutine wrf_error_fatal(msg)
      character(len=*), intent(in) :: msg
      write(*,'(a)') 'FATAL (upstream urban()): '//trim(msg)
      error stop 'wrf_stubs: wrf_error_fatal'
   end subroutine wrf_error_fatal

end module module_wrf_error

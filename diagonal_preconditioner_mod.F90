!-------------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-------------------------------------------------------------------------------


module diagonal_preconditioner_mod
  use constants_mod, only : r_def, i_def
  use preconditioner_mod, only : abstract_preconditioner_type
  use vector_mod,         only : abstract_vector_type
  use linear_operator_mod, only : abstract_linear_operator_type
  use random_operator_mod, only : random_operator_type
  use line_vector_mod,    only : line_vector_type
  use log_mod,            only : log_event, LOG_LEVEL_ERROR, &
                                 log_scratch_space

  implicit none
  private
  type, public, extends(abstract_preconditioner_type) :: diagonal_preconditioner_type
     private
     integer(kind=i_def) :: ndata
     real(kind=r_def), allocatable, dimension(:) :: op_data
   contains
     procedure, public  :: apply => apply_diagonal
     procedure, private :: apply_diagonal
     final              :: destroy_diag_precon
  end type diagonal_preconditioner_type

  interface diagonal_preconditioner_type
     module procedure diagonal_preconditioner_constructor
  end interface diagonal_preconditioner_type

contains

  function diagonal_preconditioner_constructor(nsize,op) result(self)
    implicit none
    integer(kind=i_def), intent(in) :: nsize
    class(abstract_linear_operator_type), intent(in) :: op
    type(diagonal_preconditioner_type) :: self
    self%ndata = nsize
    allocate(self%op_data(self%ndata))
    select type(op)
    type is(random_operator_type)
       call op%extract_diagonal(self%op_data)
    class default
       self%op_data(:) = 1.0_r_def
    end select
  end function diagonal_preconditioner_constructor

  subroutine apply_diagonal(self, x, y)
    implicit none    
    class(diagonal_preconditioner_type), intent(in)    :: self
    class(abstract_vector_type),         intent(inout) :: x
    class(abstract_vector_type),         intent(inout) :: y

    integer(kind = i_def) :: nlp

    select type(x)
    type is(line_vector_type)
       select type(y)
       type is(line_vector_type)
          do nlp = 1, self%ndata
             y%vdata(nlp) = self%op_data(nlp) * x%vdata(nlp)
          end do
       class default
          write(log_scratch_space,'(A)') "diagonal_preconditioner: type of y is not line_vector_type"
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
       end select
    class default
       write(log_scratch_space,'(A)') "diagonal_preconditioner: type of x is not line_vector_type"
       call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end select
  end subroutine apply_diagonal
  
  subroutine destroy_diag_precon(self)
    implicit none
    type(diagonal_preconditioner_type), intent(inout) :: self
    if(allocated(self%op_data)) then
       deallocate(self%op_data)
    end if
  end subroutine destroy_diag_precon
end module diagonal_preconditioner_mod

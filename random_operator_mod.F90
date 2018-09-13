!-------------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-------------------------------------------------------------------------------
module random_operator_mod
  use constants_mod,       only : r_def, i_def
  use linear_operator_mod, only : abstract_linear_operator_type
  use vector_mod,          only : abstract_vector_type
  use line_vector_mod,     only : line_vector_type
  use log_mod,             only : log_event, LOG_LEVEL_ERROR, log_scratch_space
  implicit none
  private

  type, public, extends(abstract_linear_operator_type) :: random_operator_type
     private
     integer(kind=i_def) :: ndata
     real(kind=r_def), allocatable, dimension(:,:) :: op_data
   contains
     procedure, public  :: apply => apply_random_op
     procedure, private :: apply_random_op
     procedure, public  :: extract_diagonal
     final              :: destroy_random_op
  end type random_operator_type

  interface random_operator_type
     module procedure random_operator_constructor
  end interface random_operator_type

contains
  function random_operator_constructor(nsize) result(self)
    implicit none
    integer(kind=i_def), intent(in) :: nsize
    type(random_operator_type) :: self
    integer(kind=i_def) :: row, col

    self%ndata = nsize
    allocate( self%op_data(self%ndata, self%ndata) )

    self%op_data(:,:) = 0.0_r_def
    do row = 1, self%ndata
       do col = row, self%ndata
          self%op_data(row,col) = real(rand(), kind=r_def)
          if(col>row) then
             self%op_data(col,row) = self%op_data(row,col)
          end if
       end do
    end do
    
  end function random_operator_constructor

  subroutine destroy_random_op(self)
    implicit none
    type(random_operator_type), intent(inout) :: self
    if(allocated(self%op_data)) then
       deallocate(self%op_data)
    end if
  end subroutine destroy_random_op

  subroutine apply_random_op(self, x, y)

    implicit none    
    class(random_operator_type), intent(in)    :: self
    class(abstract_vector_type),          intent(in)    :: x
    class(abstract_vector_type),          intent(inout) :: y

    integer(kind = i_def) :: row, col

    select type(x)
    type is(line_vector_type)
       select type(y)
       type is(line_vector_type)
          call y%set_scalar(0.0_r_def)
          do row = 1, self%ndata
             do col = 1, self%ndata
                y%vdata(row) = y%vdata(row) + self%op_data(row, col) * x%vdata(col)
             end do
          end do
       class default
          write(log_scratch_space,'(A)') "dense_operator type of y is not line_vector_type"
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
       end select
    class default
       write(log_scratch_space,'(A)') "dense_operator: type of x is not line_vector_type"
       call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end select
    
  end subroutine apply_random_op

  subroutine extract_diagonal(self, x)

    implicit none    
    class(random_operator_type), intent(in)    :: self
    !    class(abstract_vector_type), intent(inout) :: x
    real(kind=r_def), intent(inout):: x(self%ndata)

    integer(kind = i_def) :: row

       do row = 1, self%ndata
          x(row) = self%op_data(row, row)
       end do

  end subroutine extract_diagonal
end module random_operator_mod

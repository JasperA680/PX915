module vehicle_mod
    !--------------------------------------------------------------------
    ! Vehicle / cell state codes for the networked TASEP.
    !
    ! A single integer per lattice cell encodes both occupancy and
    ! turning intention.  V_EMPTY means the cell is empty; any other
    ! code means the cell is occupied by a vehicle with that intention.
    !--------------------------------------------------------------------
    implicit none
    private

    integer, parameter, public :: V_EMPTY    = 0
    integer, parameter, public :: V_STRAIGHT = 1
    integer, parameter, public :: V_LEFT     = 2
    integer, parameter, public :: V_RIGHT    = 3

    public :: is_occupied, sample_indicator

contains

    pure function is_occupied(cell) result(occ)
        integer, intent(in) :: cell
        logical :: occ
        occ = (cell /= V_EMPTY)
    end function is_occupied

    function sample_indicator(p_left, p_right) result(code)
        ! Stochastic indicator assignment used at lane entry.
        ! P(straight) = 1 - p_left - p_right.
        real, intent(in) :: p_left, p_right
        integer :: code
        real :: r

        call random_number(r)
        if (r < p_left) then
            code = V_LEFT
        else if (r < p_left + p_right) then
            code = V_RIGHT
        else
            code = V_STRAIGHT
        end if
    end function sample_indicator

end module vehicle_mod

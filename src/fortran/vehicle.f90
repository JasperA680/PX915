module vehicle_mod
    !--------------------------------------------------------------------
    ! Vehicle turning-intent codes for the networked model.
    !
    ! Occupancy is stored separately in cell%has_car; cell%turning_intent
    ! uses these codes when a vehicle is present.
    !--------------------------------------------------------------------
    implicit none
    private

    integer, parameter, public :: V_EMPTY    = 0
    integer, parameter, public :: V_OCCUPIED = 1

    public :: is_occupied

contains

    pure function is_occupied(cell) result(occ)
        integer, intent(in) :: cell
        logical :: occ
        occ = (cell /= V_EMPTY)
    end function is_occupied

end module vehicle_mod

module vehicle_mod
    !--------------------------------------------------------------------
    ! Vehicle occupancy codes for the networked model.
    ! Occupancy is tracked by cell%has_car; V_EMPTY/V_OCCUPIED are used
    ! by the junction holding-cell detection logic.
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

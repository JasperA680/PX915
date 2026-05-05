module vehicle_mod
    !--------------------------------------------------------------------
    ! Vehicle / cell state codes for the networked TASEP.
    !
    ! V_EMPTY = 0 (cell vacant); V_OCCUPIED = 1 (cell occupied).
    ! Routing decisions live on the junction (in_routes), not the cell.
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

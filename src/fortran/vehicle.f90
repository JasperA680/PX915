module vehicle_mod
    ! Vehicle occupancy constants for the networked cellular automaton model.
    !
    ! This module defines simple integer codes used when junction routines need
    ! to reason about whether a holding cell is empty or occupied.
    !
    ! Most of the network model stores occupancy using the logical field
    ! ``cell%has_car``. These integer constants are mainly used when junction
    ! logic builds compact arrays of holding-cell states, for example:
    !
    ! .. math::
    !
    !    H_k =
    !    \begin{cases}
    !    0, & \text{if inbound holding cell } k \text{ is empty}, \\
    !    1, & \text{if inbound holding cell } k \text{ is occupied}.
    !    \end{cases}
    !
    ! The helper function ``is_occupied`` converts one of these integer occupancy
    ! codes into a logical value.

    implicit none
    private

    integer, parameter, public :: V_EMPTY = 0
    ! Integer code representing an empty holding cell.
    integer, parameter, public :: V_OCCUPIED = 1
    ! Integer code representing an occupied holding cell.

    public :: is_occupied

contains

    pure function is_occupied(cell) result(occ)
        ! Return whether an integer occupancy code represents an occupied cell.
        !
        ! The function returns ``.true.`` for any value that is not ``V_EMPTY``.
        ! This makes it robust to future extensions where more than one non-empty
        ! vehicle code might be introduced.
        integer, intent(in) :: cell
        ! Integer occupancy code to test.
        logical :: occ
        ! True if ``cell`` is not equal to ``V_EMPTY``.

        occ = (cell /= V_EMPTY)
    end function is_occupied

end module vehicle_mod
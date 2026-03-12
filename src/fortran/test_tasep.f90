program test_tasep 
    use tasep_model 
    implicit none 

    integer, parameter :: L = 10 
    integer :: state(L)
    integer :: step, exit_count 

    call initialise_lattice(state, L)

    do step = 1, 20 
        call tasep_step(state, L, 0.5, 0.5, exit_count)
        print *, "Step:", step, " State:", state, "Density:", compute_density(state, L), "Exited:", exit_count
    end do

end program test_tasep
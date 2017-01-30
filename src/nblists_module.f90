module nblists

type, bind(C) :: tOpts
  logical(lb) :: thirdLaw
  logical(lb) :: jointXYZ
  logical(lb) :: zeroBase
end type tOpts

type, bind(C) :: nbList
  type(c_ptr) :: start                     ! Location of the first neighbor of each atom
  type(c_ptr) :: end                       ! Location of the last neighbor of each atom
  type(c_ptr) :: item                      ! Indices of neighbor atoms
  integer(ib) :: nitems                    ! Number of atoms in neighbor list
  integer(ib) :: nmax                      ! Maximum number of atoms in neighbor list
  integer(ib) :: builds                    ! Number of neighbor list builds
  real(rb)    :: time                      ! Time taken in neighbor list handling
  type(tOpts) :: Options                   ! List of options
  type(c_ptr) :: Data                      ! Pointer to system data
end type nbList

interface

  function neighbor_list( threads, rc, skin, N, group ) bind(C,name="neighbor_list")
    import
    integer(ib), value :: threads, N
    real(rb),    value :: rc, skin
    type(c_ptr), value :: group
    type(nbList)       :: neighbor_list
  end function neighbor_list

  subroutine neighbor_allocate_array( list, array, jointXYZ ) bind(C,name="neighbor_allocate_array")
    import
    type(nbList), intent(inout) :: list
    type(c_ptr),  intent(inout) :: array
    logical(lb),  value         :: jointXYZ
  end subroutine neighbor_allocate_array

  function neighbor_list_outdated( list, coords, N, atoms ) bind(C,name="neighbor_list_outdated")
    import
    type(nbList), intent(inout) :: list
    type(c_ptr),  value         :: coords
    integer(ib),  value         :: N
    integer(ib),  intent(in)    :: atoms(*)
    logical(lb)                 :: neighbor_list_outdated
  end function neighbor_list_outdated

  subroutine neighbor_list_build( list, Lbox, positions ) bind(C,name="neighbor_list_build")
    import
    type(nbList), intent(inout) :: list
    real(rb),     value         :: L
    type(c_ptr),  intent(in)    :: coords
  end subroutine neighbor_list_build

end interface

end module nblists

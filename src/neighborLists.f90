module neighborLists

use global
use lists
use omp_lib

implicit none

integer, parameter, private :: extra = 500

integer, parameter, private :: ndiv = 2
integer, parameter, private :: nbcells = 62
integer, parameter, private :: nb(3,nbcells) = reshape( [ &
   1, 0, 0,    2, 0, 0,   -2, 1, 0,   -1, 1, 0,    0, 1, 0,    1, 1, 0,    2, 1, 0,   -2, 2, 0,  &
  -1, 2, 0,    0, 2, 0,    1, 2, 0,    2, 2, 0,   -2,-2, 1,   -1,-2, 1,    0,-2, 1,    1,-2, 1,  &
   2,-2, 1,   -2,-1, 1,   -1,-1, 1,    0,-1, 1,    1,-1, 1,    2,-1, 1,   -2, 0, 1,   -1, 0, 1,  &
   0, 0, 1,    1, 0, 1,    2, 0, 1,   -2, 1, 1,   -1, 1, 1,    0, 1, 1,    1, 1, 1,    2, 1, 1,  &
  -2, 2, 1,   -1, 2, 1,    0, 2, 1,    1, 2, 1,    2, 2, 1,   -2,-2, 2,   -1,-2, 2,    0,-2, 2,  &
   1,-2, 2,    2,-2, 2,   -2,-1, 2,   -1,-1, 2,    0,-1, 2,    1,-1, 2,    2,-1, 2,   -2, 0, 2,  &
  -1, 0, 2,    0, 0, 2,    1, 0, 2,    2, 0, 2,   -2, 1, 2,   -1, 1, 2,    0, 1, 2,    1, 1, 2,  &
   2, 1, 2,   -2, 2, 2,   -1, 2, 2,    0, 2, 2,    1, 2, 2,    2, 2, 2 ], shape(nb) )

type, private :: tCell
  integer :: neighbor(nbcells)
end type tCell

type, private :: tData

  logical :: MonteCarlo                   ! True if Monte Carlo instead of Molecular Dynamics
  logical :: zeroBasedIndexing            ! True/False if atom indexing starts with 0/1
  integer :: natoms                       ! Number of atoms in the system
  integer :: mcells = 0                   ! Number of cells at each dimension
  integer :: ncells = 0                   ! Total number of cells
  integer :: maxcells = 0                 ! Maximum number of cells
  integer :: maxatoms = 0                 ! Maximum number of atoms in a cell
  integer :: maxpairs = 0                 ! Maximum number of pairs formed by all atoms of a cell
  integer :: nthreads                     ! Number of parallel openmp threads
  integer :: threadAtoms                  ! Number of atoms per parallel thread

  real(rb) :: Rc                          ! Cut-off distance
  real(rb) :: RcSq                        ! Cut-off distance squared
  real(rb) :: xRc                         ! Extended cutoff distance (including skin)
  real(rb) :: xRcSq                       ! Extended cutoff distance squared
  real(rb) :: skinSq                      ! Square of the neighbor list skin width

  type(tList) :: cellAtom                 ! List of atoms belonging to each cell
  type(tList) :: threadCell               ! List of cells to be dealt with in each parallel thread

  integer,     allocatable :: atomCell(:) ! Array containing the current cell of each atom
  integer,     allocatable :: group(:)    ! Array containing the group index of each atom
  real(rb),    allocatable :: R0(:,:)     ! Position of each atom at latest neighbor list building
  type(tCell), allocatable :: cell(:)     ! Array containing all neighbor cells of each cell
  type(tList), allocatable :: neighbor(:) ! Pointer to neighbor lists

end type tData

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

contains

!===================================================================================================
!                                L I B R A R Y   P R O C E D U R E S
!===================================================================================================

  function neighbor_list( threads, rc, skin, N, group ) bind(C,name="neighbor_list")
    integer(ib), value :: threads, N
    real(rb),    value :: rc, skin
    type(c_ptr), value :: group
    type(nbList)       :: neighbor_list

    integer :: i

    type(tData), pointer :: me
    integer,     pointer :: pgroup(:), start(:), end(:), item(:)

    ! Allocate data structure:
    allocate( me )
    allocate( start(N), end(N), item(N), source = 0 )
    allocate( me%R0(3,N), source = zero )

    ! Set up fixed entities:
    me%nthreads = threads
    me%Rc = rc
    me%RcSq = rc*rc
    me%xRc = rc + skin
    me%xRcSq = me%xRc**2
    me%skinSq = skin*skin
    me%natoms = N
    me%threadAtoms = (N + threads - 1)/threads

    ! Set up group:
    if (c_associated(group)) then
      call c_f_pointer( group, pgroup, [N] )
      allocate( me%group(N), source = pgroup )
    else
      allocate( me%group(N), source = [(i,i=1,N)] )
    end if

    ! Initialize counters and other mutable entities:
    allocate( me%cell(0), me%atomCell(N) )

    ! Allocate memory for list of atoms per cell:
    call me % cellAtom % allocate( N, 0 )

    ! Allocate memory for neighbor lists:
    allocate( me%neighbor(threads) )
    call me % neighbor % allocate( extra, N )

    ! Deploy data structure:
    neighbor_list % data = c_loc(me)
    neighbor_list % start = c_loc(start(1))
    neighbor_list % end = c_loc(end(1))
    neighbor_list % item = c_loc(item(1))
    neighbor_list % nitems = N
    neighbor_list % nmax = N
    neighbor_list % builds = 0
    neighbor_list % time = zero
    neighbor_list % options % thirdLaw = .true.
    neighbor_list % options % jointXYZ = .true.
    neighbor_list % options % zeroBase = .true.

  end function neighbor_list

!===================================================================================================

  subroutine neighbor_allocate_2d_array( list, array ) bind(C,name="neighbor_allocate_2d_array")
    type(nbList), value         :: list
    type(c_ptr),  intent(inout) :: array

    integer :: i
    type(tData), pointer :: me
    real(rb),    pointer :: R(:,:)
    type(c_ptr), pointer :: ptr(:)

    call c_f_pointer( list%data, me )
    if (list % options % jointXYZ) then
      allocate( R(3,me%natoms), ptr(me%natoms) )
      forall (i=1:me%natoms) ptr(i) = c_loc(R(1,i))
    else
      allocate( R(me%natoms,3), ptr(3) )
      forall (i=1:me%natoms) ptr(i) = c_loc(R(i,1))
    end if
    array = c_loc(ptr(1))

  end subroutine neighbor_allocate_2d_array

!===================================================================================================

  function neighbor_list_outdated( list, coords ) bind(C,name="neighbor_list_outdated")
    type(nbList), value :: list
    type(c_ptr),  value :: coords
    logical(lb)         :: neighbor_list_outdated

    real(rb),    pointer :: R(:,:)
    type(tData), pointer :: me

    call c_f_pointer( list%data, me )
    if (list % options % jointXYZ) then
      call c_f_pointer( coords, R, [3,me%natoms] )
      neighbor_list_outdated = maximum_approach_sq( me%natoms, R - me%R0 ) > me%skinSq
    else
      call c_f_pointer( coords, R, [me%natoms,3] )
      neighbor_list_outdated = maximum_approach_sq( me%natoms, transpose(R) - me%R0 ) > me%skinSq
    end if

  end function neighbor_list_outdated

!===================================================================================================

  subroutine neighbor_list_build( list, L, coords ) bind(C,name="neighbor_list_build")
    type(nbList), intent(inout) :: list
    real(rb),     value         :: L
    type(c_ptr),  value         :: coords

    integer  :: M, i, n
    real(rb) :: invL, time

    integer,  allocatable :: neighbors(:,:)
    real(rb), allocatable :: Rs(:,:)

    integer,     pointer :: start(:), end(:), item(:)
    real(rb),    pointer :: R(:,:)
    type(tData), pointer :: me

    list%time = list%time - omp_get_wtime()

    call c_f_pointer( list%data, me )
    allocate( Rs(3,me%natoms), neighbors(me%natoms,me%nthreads) )

    if (list % options % jointXYZ) then
      call c_f_pointer( coords, R, [3,me%natoms] )
      me%R0 = R
    else
      call c_f_pointer( coords, R, [me%natoms,3] )
      me%R0 = transpose(R)
    end if

    invL = one/L
    Rs = invL*me%R0

    me%MonteCarlo = .not.list%options%thirdLaw
    me%zeroBasedIndexing = list%options%zeroBase

    list%builds = list%builds + 1

    M = floor(ndiv*L/me%xRc)
    call distribute_atoms( me, max(M,2*ndiv+1), Rs )

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread
      thread = omp_get_thread_num() + 1
      call find_pairs( me, thread, InvL, Rs )
      call count_neighbors( me, thread, neighbors(:,thread) )
    end block
    !$omp end parallel

    call c_f_pointer( list%start, start, [me%natoms] )
    call c_f_pointer( list%end, end, [me%natoms] )
    n = 0
    do i = 1, me%natoms
      start(i) = n + 1
      n = n + sum(neighbors(i,:))
      end(i) = n
    end do
    if (list%nmax < n) call reallocate_item_array( n )
    list%nitems = n
    call c_f_pointer( list%item, item, [n] )

    call deploy_neighbor_list( me, start, end, item )

    time = omp_get_wtime()
    list%time = list%time + time

    contains

      subroutine reallocate_item_array( n )
        integer, intent(in) :: n
        integer, pointer :: aux(:)
        call c_f_pointer( list%item, aux, [list%nmax] )
        deallocate( aux )
        allocate( aux(n) )
        list%item = c_loc(aux(1))
        list%nmax = n
      end subroutine reallocate_item_array

  end subroutine neighbor_list_build

!===================================================================================================
!                              A U X I L I A R Y   P R O C E D U R E S
!===================================================================================================

  pure subroutine count_neighbors( me, thread, ncount )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: thread
    integer,     intent(out)   :: ncount(me%natoms)

    integer :: firstAtom, lastAtom, m, i, j, k

    ncount = 0
    firstAtom = me%cellAtom%first(me%threadCell%first(thread))
    lastAtom = me%cellAtom%last(me%threadCell%last(thread))
    associate (neighbor => me%neighbor(thread) )
      do m = firstAtom, lastAtom
        i = me%cellAtom%item(m)
        ncount(i) = ncount(i) + neighbor%last(i) - neighbor%first(i) + 1
        if (me%MonteCarlo) then
          do k = neighbor%first(i), neighbor%last(i)
            j = neighbor%item(k)
            ncount(j) = ncount(j) + 1
          end do
        end if
      end do
    end associate

  end subroutine count_neighbors

!===================================================================================================

  subroutine deploy_neighbor_list( me, start, end, item )
    type(tData), intent(inout) :: me
    integer,     intent(inout) :: start(me%natoms), end(me%natoms)
    integer,     intent(out)   :: item(:)

    integer :: position(me%natoms)

    position = start - 1

    !$omp parallel num_threads(me%nthreads)
    call build_partial_list
    !$omp end parallel

    contains

      subroutine build_partial_list
        integer :: thread, first, last, m, i, j, k, n, ipos
        thread = omp_get_thread_num() + 1

        first = me%cellAtom%first(me%threadCell%first(thread))
        last = me%cellAtom%last(me%threadCell%last(thread))
        associate (neighbor => me%neighbor(thread) )
          do m = first, last
            i = me%cellAtom%item(m)
            n = neighbor%last(i) - neighbor%first(i) + 1
            !$omp atomic capture
            ipos = position(i)
            position(i) = position(i) + n
            !$omp end atomic
            do k = neighbor%first(i), neighbor%last(i)
              j = neighbor%item(k)
              ipos = ipos + 1
              item(ipos) = j
              if (me%MonteCarlo) then
                !$omp atomic update
                position(j) = position(j) + 1
                !$omp end atomic
                item(position(j)) = i
              end if
            end do
          end do
        end associate
        !$omp barrier

        if (me%zeroBasedIndexing) then
          first = (thread - 1)*me%threadAtoms + 1
          last = min(thread*me%threadAtoms, me%natoms)
          forall (i=start(first):end(last)) item(i) = item(i) - 1
          start(first:last) = start(first:last) - 1
          end(first:last) = end(first:last) - 1
        end if

      end subroutine build_partial_list

  end subroutine deploy_neighbor_list

!===================================================================================================

  real(rb) function maximum_approach_sq( N, delta )
    integer,  intent(in) :: N
    real(rb), intent(in) :: delta(3,N)

    integer  :: i
    real(rb) :: maximum, next, deltaSq

    maximum = sum(delta(:,1)**2)
    next = maximum
    do i = 2, N
      deltaSq = sum(delta(:,i)**2)
      if (deltaSq > maximum) then
        next = maximum
        maximum = deltaSq
      end if
    end do
    maximum_approach_sq = maximum + 2*sqrt(maximum*next) + next

  end function maximum_approach_sq

!===================================================================================================

  pure subroutine insertion_sort( val, ind )
    real(rb), intent(inout) :: val(:)
    integer,  intent(inout) :: ind(:)

    integer :: i, j, itemp
    real(rb) :: rtemp
 
    do i = 2, size(val)
      j = i - 1
      rtemp = val(i)
      itemp = ind(i)
      do while ( (j >= 1).and.(val(j) > rtemp) )
        val(j+1) = val(j)
        ind(j+1) = ind(j)
        j = j - 1
      end do
      val(j+1) = rtemp
      ind(j+1) = itemp
    end do
  end subroutine insertion_sort

!===================================================================================================

  subroutine distribute_atoms( me, M, Rs )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: M
    real(rb),    intent(in)    :: Rs(3,me%natoms)

    integer :: MM, cells_per_thread, maxNatoms, threadNatoms(me%nthreads), next(me%natoms)
    logical :: make_cells
    integer, allocatable :: natoms(:)

    MM = M*M
    make_cells = M /= me%mcells
    if (make_cells) then
      me%mcells = M
      me%ncells = M*MM
      if (me%ncells > me%maxcells) then
        deallocate( me%cell, me%cellAtom%first, me%cellAtom%last )
        allocate( me%cell(me%ncells), me%cellAtom%first(me%ncells), me%cellAtom%last(me%ncells) )
        call me % threadCell % allocate( 0, me%nthreads )
        me%maxcells = me%ncells
      end if
      cells_per_thread = (me%ncells + me%nthreads - 1)/me%nthreads
    end if

    allocate( natoms(me%ncells) )

    !$omp parallel num_threads(me%nthreads) reduction(max:maxNatoms)
    call distribute( omp_get_thread_num() + 1, maxNatoms )
    !$omp end parallel
    me%maxatoms = maxNatoms
    me%maxpairs = (maxNatoms*((2*nbcells + 1)*maxNatoms - 1))/2

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      subroutine distribute( thread, maxNatoms )
        integer, intent(in)  :: thread
        integer, intent(out) :: maxNatoms

        integer :: i, k, icell, ix, iy, iz, first, last
        integer :: icoord(3)
        integer, allocatable :: head(:)

        if (make_cells) then
          first = (thread - 1)*cells_per_thread + 1
          last = min( thread*cells_per_thread, me%ncells )
          do icell = first, last
            k = icell - 1
            iz = k/MM
            iy = (k - iz*MM)/M
            ix = k - (iy*M + iz*MM)
            me%cell(icell)%neighbor = 1 + pbc(ix+nb(1,:)) + pbc(iy+nb(2,:))*M + pbc(iz+nb(3,:))*MM
          end do
          me%threadCell%first(thread) = first
          me%threadCell%last(thread) = last
        else
          first = me%threadCell%first(thread)
          last = me%threadCell%last(thread)
        end if

        do i = (thread - 1)*me%threadAtoms + 1, min( thread*me%threadAtoms, me%natoms )
          icoord = int(M*(Rs(:,i) - floor(Rs(:,i))),ib)
          me%atomCell(i) = 1 + icoord(1) + M*icoord(2) + MM*icoord(3)
        end do
        !$omp barrier

        allocate( head(first:last) )
        head = 0
        natoms(first:last) = 0
        do i = 1, me%natoms
          icell = me%atomCell(i)
          if ((icell >= first).and.(icell <= last)) then
            next(i) = head(icell)
            head(icell) = i
            natoms(icell) = natoms(icell) + 1
          end if
        end do
        threadNatoms(thread) = sum(natoms(first:last))
        !$omp barrier

        maxNatoms = 0
        k = sum(threadNatoms(1:thread-1))
        do icell = first, last
          me%cellAtom%first(icell) = k + 1
          i = head(icell)
          do while (i /= 0)
            k = k + 1
            me%cellAtom%item(k) = i
            i = next(i)
          end do
          me%cellAtom%last(icell) = k
          if (natoms(icell) > maxNatoms) maxNatoms = natoms(icell)
        end do
      end subroutine distribute
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      elemental integer function pbc( x )
        integer, intent(in) :: x
        if (x < 0) then
          pbc = x + M
        else if (x >= M) then
          pbc = x - M
        else
          pbc = x
        end if
      end function pbc
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine distribute_atoms

!===================================================================================================

  subroutine find_pairs( me, thread, invL, Rs )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: thread
    real(rb),    intent(in)    :: invL, Rs(3,me%natoms)

    integer  :: i, j, k, m, n, icell, jcell, npairs, igroup
    integer  :: nlocal, ntotal, first, last
    real(rb) :: xRc2, Rc2, r2
    integer  :: atom(me%maxpairs)
    real(rb) :: Ri(3), Rij(3), invL2

    invL2 = invL*invL
    xRc2 = me%xRcSq*invL2
    Rc2 = me%RcSq*invL2

    npairs = 0
    associate (neighbor => me%neighbor(thread))
      do icell = me%threadCell%first(thread), me%threadCell%last(thread)

        if (neighbor%nitems < npairs + me%maxpairs) then
          call neighbor % resize( npairs + me%maxpairs + extra )
        end if

        first = me%cellAtom%first(icell)
        last = me%cellAtom%last(icell)
        nlocal = last - first + 1
        atom(1:nlocal) = me%cellAtom%item(first:last)

        ntotal = nlocal
        do m = 1, nbcells
          jcell = me%cell(icell)%neighbor(m)
          first = me%cellAtom%first(jcell)
          last = me%cellAtom%last(jcell)
          n = ntotal + 1
          ntotal = n + last - first
          atom(n:ntotal) = me%cellAtom%item(first:last)
        end do

        do k = 1, nlocal
          i = atom(k)
          igroup = me%group(i)
          neighbor%first(i) = npairs + 1
          Ri = Rs(:,i)
          do m = k + 1, ntotal
            j = atom(m)
            if (me%group(j) /= igroup) then
              Rij = Ri - Rs(:,j)
              Rij = Rij - anint(Rij)
              r2 = sum(Rij*Rij)
              if (r2 < xRc2) then
                npairs = npairs + 1
                neighbor%item(npairs) = j
              end if
            end if
          end do
          neighbor%last(i) = npairs
        end do

      end do
      neighbor%count = npairs
    end associate

  end subroutine find_pairs

!===================================================================================================

end module neighborLists

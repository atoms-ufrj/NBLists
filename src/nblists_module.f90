module nblists

type, bind(C) :: tOpts
  logical(lb) :: thirdLaw
  logical(lb) :: jointXYZ
end type tOpts

type, bind(C) :: nbList
  type(c_ptr) :: first                     ! Location of the first neighbor of each atom
  type(c_ptr) :: last                      ! Location of the last neighbor of each atom
  type(c_ptr) :: item                      ! Indices of neighbor atoms
  type(c_ptr) :: R0                        ! Atom positions at latest neighbor list build
  integer(ib) :: nitems                    ! Number of atoms in neighbor list
  integer(ib) :: nmax                      ! Maximum number of atoms in neighbor list
  integer(ib) :: builds                    ! Number of neighbor list builds
  real(rb)    :: time                      ! Time taken in neighbor list handling
  type(c_ptr) :: Data                      ! Pointer to system data
  type(tOpts) :: Options                   ! List of options
end type nbList

interface

  function neighbor_list( threads, rc, skin, N, body ) bind(C,name="neighbor_list")
    import
    integer(ib), value :: threads, N
    real(rb),    value :: rc, skin
    type(c_ptr), value :: body
    type(nbList)       :: neighbor_list
  end function neighbor_list

!===================================================================================================

  function neighbor_list_outdated( list, positions ) result( outdated ) &
    bind(C,name="neighbor_list_outdated")
    import
    type(nbList), intent(inout) :: list
    type(c_ptr),  value         :: positions
    logical(lb)                 :: outdated

    type(tData), pointer :: me
    real(rb),    pointer :: R(:,:)

    call c_f_pointer( list%data, me )
    call c_f_pointer( positions, R, [3,me%natoms] )
    if (list % options % jointXYZ) then
      outdated = maximum_approach_sq( me%natoms, R - me%R0 ) > me%skinSq
    else
      outdated = maximum_approach_sq( me%natoms, transpose(R) - me%R0 ) > me%skinSq
    end if

  end function neighbor_list_outdated

!===================================================================================================

  subroutine neighbor_list_build( list, Lbox, positions ) &
    bind(C,name="neighbor_list_build")
    type(nbList), intent(inout) :: list
    real(rb),     value         :: Lbox
    type(c_ptr),  value         :: positions

    integer  :: M, i, n
    real(rb) :: invL, time
    integer, allocatable  :: neighbors(:,:)
    integer,  pointer     :: first(:), last(:), item(:)
    real(rb), pointer     :: R(:,:)
    real(rb), allocatable :: Rs(:,:)
    type(tData),  pointer :: me

    list%time = list%time - omp_get_wtime()

    call c_f_pointer( list%data, me )
    call c_f_pointer( positions, R, [3,me%natoms] )
    call c_f_pointer( list%first, first, [me%natoms] )
    call c_f_pointer( list%last, last, [me%natoms] )
    call c_f_pointer( list%item, item, [list%nmax] )
    allocate( Rs(3,me%natoms), neighbors(me%natoms,me%nthreads) )
    me%MonteCarlo = .not.list%options%thirdLaw

    list%builds = list%builds + 1
    if (list%options%jointXYZ) then
      me%R0 = R
    else
      me%R0 = transpose(R)
    end if

    M = floor(ndiv*Lbox/me%xRc)
    invL = one/Lbox
    Rs = invL*R
    call distribute_atoms( me, max(M,2*ndiv+1), Rs )

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread
      thread = omp_get_thread_num() + 1
      call find_pairs( me, thread, InvL, Rs )
      call count_neighbors( me, thread, neighbors(:,thread) )
    end block
    !$omp end parallel

    n = 0
    do i = 1, me%natoms
      first(i) = n + 1
      n = n + sum(neighbors(i,:))
      last(i) = n
    end do

    if (list%nmax < n) then
      deallocate( item )
      allocate( item(n) )
      list%nmax = n
      list%item = c_loc(item)
    end if

    call deploy_neighbor_list( me, first, item )

    time = omp_get_wtime()
    list%time = list%time + time

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

  subroutine deploy_neighbor_list( me, first, item )
    type(tData), intent(inout) :: me
    integer,     intent(in)    :: first(me%natoms)
    integer,     intent(out)   :: item(:)

    integer :: position(me%natoms)

    position = first - 1

    !$omp parallel num_threads(me%nthreads)
    block
      integer :: thread, firstAtom, lastAtom, m, i, j, k
      thread = omp_get_thread_num() + 1

      firstAtom = me%cellAtom%first(me%threadCell%first(thread))
      lastAtom = me%cellAtom%last(me%threadCell%last(thread))
      associate (neighbor => me%neighbor(thread) )
        do m = firstAtom, lastAtom
          i = me%cellAtom%item(m)
          do k = neighbor%first(i), neighbor%last(i)
            j = neighbor%item(k)
            !$omp atomic update
            position(i) = position(i) + 1
            !$omp end atomic
            item(position(i)) = j
            if (me%MonteCarlo) then
              !$omp atomic update
              position(j) = position(j) + 1
              !$omp end atomic
              item(position(j)) = i
            end if
          end do
        end do
      end associate
    end block
    !$omp end parallel

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

    integer  :: i, j, k, m, n, icell, jcell, npairs, ibody
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
          ibody = me%body(i)
          neighbor%first(i) = npairs + 1
          Ri = Rs(:,i)
          do m = k + 1, ntotal
            j = atom(m)
            if (me%body(j) /= ibody) then
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

end module nblists


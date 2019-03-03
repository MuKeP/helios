    module harvester

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MODULES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    use glob,      only: r8kind,i4kind,void,uch,date_time
    use fcontrol,  only: fcNewID,fcNullID
    use txtParser, only: tpAdjustc
    use printmod,  only: prJoin

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ CONSTANTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    character (len=*), parameter :: hrVersion='1.000'
    character (len=*), parameter :: hrDate   ='20-July-2018'
    character (len=*), parameter :: hrAuthor ='Anton B. Zakharov'

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TYPES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    type hpolariz
        character(len=25)    :: date
        character(len=10)    :: version

        character(len=15)    :: method
        integer(kind=i4kind) :: state

        integer(kind=i4kind) :: nAtoms,nElectrons

        integer(kind=i4kind) :: points
        real(kind=r8kind)    :: field
        character(len=3)     :: scales

        real(kind=r8kind)    :: m(4)  !mx,my,mz,mmod
        real(kind=r8kind)    :: a(4)  !axx,xyy,azz,ava
        real(kind=r8kind)    :: b(10) !bxxx,byyy,bzzz,bxyy,byxx,bxzz,bzxx,byzz,bzyy,bmod
        real(kind=r8kind)    :: g(10) !gxxxx,gyyyy,gzzzz,gxxyy,gyyxx,gxxzz,gzzxx,gyyzz,gzzyy,avg

        real(kind=r8kind)    :: sigmaContributionAlpha,sigmaContributionGamma

        real(kind=r8kind), allocatable :: Energies(:,:,:)
    end type hpolariz

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ARRAYS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIABLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=1) :: probe

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ INTERFACES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    interface init_harvester
        module procedure init_harvester_polariz
    end interface init_harvester

    interface store_harvester
        module procedure store_harvester_polariz
    end interface store_harvester

    ! interface read_harvester
    !     module procedure read_harvester_polariz
    ! end interface read_harvester

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ACCESS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    private

    public :: hpolariz
    public :: init_harvester,store_harvester,collect_harvester

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    contains

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine init_harvester_polariz(harv,uver,umeth,ustate,un,une,upnt,ufield,uscales)
    implicit none

    type(hpolariz)                   :: harv
    character (len=*)   , intent(in) :: uver,umeth,uscales
    integer(kind=i4kind), intent(in) :: ustate,un,une,upnt
    real(kind=r8kind)   , intent(in) :: ufield

    integer(kind=i4kind) :: sta,sto


    harv%date=date_time()

    harv%version=repeat(' ',len(harv%version)); harv%version=trim(uver)
    harv%method=repeat(' ',len(harv%method)); harv%method=umeth
    harv%scales=repeat(' ',len(harv%scales)); harv%scales=uscales

    harv%state=ustate
    harv%nAtoms=un
    harv%nElectrons=une
    harv%points=upnt
    harv%field=ufield

    sta=-(upnt-1)/2; sto=(upnt-1)/2

    if (allocated(harv%Energies)) deallocate(harv%Energies)

    allocate (harv%Energies(sta:sto,sta:sto,sta:sto))
    harv%Energies = 0

    harv%m=0; harv%a=0; harv%b=0; harv%g=0
    harv%sigmaContributionAlpha=0
    harv%sigmaContributionGamma=0

    return
    end subroutine init_harvester_polariz

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine store_harvester_polariz(harv,file)
    implicit none

    type(hpolariz)      , intent(in) :: harv
    character(len=*)    , intent(in) :: file
    integer(kind=i4kind)             :: id,err,ln,kind


    if (.NOT.allocated(harv%Energies)) return

    id=fcNewID()
    open(id,file=file,form='binary')

    read(id,iostat=err) probe; rewind(id)
    if (err.NE.0) then
        write (id) 'hrv'
    endif

    close(id)

    open(id,file=file,form='binary',access='append')

    ln=len(harv%date)+len(harv%version)+len(harv%method)+4*i4kind+3+3*r8kind+&
       sizeof(harv%Energies)+sizeof(harv%m)+sizeof(harv%a)+sizeof(harv%b)+sizeof(harv%g)

    kind=int(100, i4kind)

    write (id) kind
    write (id) ln

    write (id) harv%date,harv%version
    write (id) harv%method,harv%state,harv%nAtoms,harv%nElectrons,harv%scales,harv%points,harv%field

    write (id) harv%Energies

    write (id) harv%m
    write (id) harv%a
    write (id) harv%b
    write (id) harv%g

    write (id) harv%sigmaContributionAlpha,harv%sigmaContributionGamma

    close (id)
    void=fcNullID(id)

    return
    end subroutine store_harvester_polariz

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    integer(kind=i4kind) function count_records_harvester(file) result(ret)
    implicit none

    character(len=*)    , intent(in) :: file
    integer(kind=i4kind)             :: id,k,kind,rlen
    character(len=3)                 :: header


    id=fcNewID()
    open(id,file=file,form='binary')
    read (id) header

    if (header.NE.'hrv') then
        ret=0
        close(id)
        void=fcNullID(id)
        return
    endif

    ! hpolariz: kind=100

    ret=0
    do
        if (eof(id)) exit

        read (id) kind
        read (id) rlen
        read (id) (probe, k = 1,rlen)
        ret=ret+1
    enddo

    close(id)
    void=fcNullID(id)

    return
    end function count_records_harvester

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine collect_harvester(file,kind,ou,accuracy)

    implicit none

    character (len=*)   , intent(in) :: file,kind
    integer(kind=i4kind), intent(in) :: ou,accuracy


    select case (kind)
        case('polarizability')
            call collect_harvester_polariz(file,ou,accuracy)
        case('density')
            continue
    end select

    return
    end subroutine collect_harvester

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine collect_harvester_polariz(file,ou,accuracy)
    implicit none

    character(len=*)    , intent(in) :: file
    integer(kind=i4kind), intent(in) :: ou,accuracy

    character(len=3)                 :: head
    integer(kind=i4kind)             :: i,ln,mln,nrecs,kind,id
    type(hpolariz), allocatable      :: hall(:)
    character(len=*), parameter      :: header(29)=['E0','X','Y','Z','|D|','XX','YY','ZZ','<A>',&
                                                    'XXX','YYY','ZZZ','XYY','YXX','XZZ','ZXX','YZZ','ZYY','|B|',&
                                                    'XXXX','YYYY','ZZZZ','XXYY','YYXX','XXZZ','ZZXX','YYZZ','ZZYY','<G>']


    nrecs=count_records_harvester(file)

    allocate (hall(nrecs))

    id=fcNewID()
    open(id,file=file,form='binary')
    read(id) head

    do i = 1,nrecs
        call read_harvester_polariz(hall(i),id)
    enddo

    close(id)
    void=fcNullID(id)

    ln=8+accuracy
    mln=len(hall(1)%method)

    write (ou,'(A)') '^'//repeat(' ',mln+1)//tpAdjustc(header(1), 22)//prJoin([(tpAdjustc(header(i), ln), i=2,29)], ' ')
    write (ou,100) '#',hall(1)%nAtoms,hall(1)%nElectrons,hall(1)%scales,hall(1)%points,hall(1)%field
    do i = 1,nrecs
        write (ou,101) '*'//tpAdjustc(hall(i)%method,mln),hall(i)%Energies(0,0,0),&
                       hall(i)%m,hall(i)%a,hall(i)%b,hall(i)%g
    enddo
    write (ou,'(/)')

    do i = 1,nrecs
        deallocate (hall(i)%Energies)
    enddo
    deallocate(hall)


100 format(A,1X,'Atoms',1X,i4,1X,'Electrons',1X,i4,1X,'Scale',1X,A,1X,'Points',1X,i2,1X,'FieldStep',1X,F6.4)
101 format(A,1X,ES20.12,28(1X,ES<8+accuracy>.<accuracy>))

    return
    end subroutine collect_harvester_polariz

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    subroutine read_harvester_polariz(harv,id)
    implicit none

    type(hpolariz)                   :: harv
    integer(kind=i4kind), intent(in) :: id
    integer(kind=i4kind)             :: ln,kind

    integer(kind=i4kind)             :: sta,sto


    read (id) kind
    read (id) ln
    read (id) harv%date,harv%version
    read (id) harv%method,harv%state,harv%nAtoms,harv%nElectrons,harv%scales,harv%points,harv%field

    sta=-(harv%points-1)/2; sto= (harv%points-1)/2

    allocate (harv%Energies(sta:sto,sta:sto,sta:sto))

    read (id) harv%Energies

    read (id) harv%m
    read (id) harv%a
    read (id) harv%b
    read (id) harv%g

    read (id) harv%sigmaContributionAlpha,harv%sigmaContributionGamma

    return
    end subroutine read_harvester_polariz

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   !

    end module harvester
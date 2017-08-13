	subroutine parseInput

	use hdb
	use datablock , only: bdParseFile
	use argsParser, only: parseArgs,apArgumentFound

	implicit none

	integer(kind=iglu)     :: err,k,rwidth=150,pwidth=70


	void=parseArgs(true)

	if (apArgumentFound('-i')) then
		!call ShowInputHelp
		stop
	endif

	if ( tpIndex( uchGet(generalbd%fname),'.inp').GT.0 ) then
		generalbd%outfile=uchSet( uchGet(generalbd%fname,1,tpIndex( uchGet(generalbd%fname),'.inp',rev=true)-1)//'.out' )
	elseif (tpIndex( uchGet(generalbd%fname),'.',rev=true).GT.0) then
		generalbd%outfile=uchSet( uchGet(generalbd%fname,1,tpIndex( uchGet(generalbd%fname),'.'   ,rev=true)-1)//'.out' )
	else
		generalbd%outfile=uchSet( uchGet(generalbd%fname)//'.out' )
	endif

	!call prEchoFile(init,ou,'input %numeration >>')
	!void=bdPrintBlockData(bdNames(k),ou,width=pwidth,ignoreFreeBlock=true)

	do k = 1,UBound(bdNames,1)
		!write (*,'(4X,i2.2,A,i2.2,1X,A)') k,'/',UBound(bdNames,1),'Parsing block data "'//trim(bdNames(k))//'"'
		err=bdParseFile(bdNames(k),uchGet(generalbd%fname),rwidth)
		if (err.NE.0) then
			write (eu,*) bdNames(k),': Error occured while parsing input file.'; stop
		endif
	enddo

	generalbd%bondsAlternated=abs(generalbd%alternation).GT.0

	! todo for all accuracies.
	geometrybd%symmetryTolerance=real(10,rglu)**-geometrybd%symmetryTolerance
	scfbd%accuracy=real(10,rglu)**-scfbd%accuracy
	fcibd%accuracy=real(10,rglu)**-fcibd%accuracy
	fcibd%zeroThreshold=real(10,rglu)**-fcibd%zeroThreshold

	

	open (ou,file=uchGet(generalbd%outfile))
	open (init,file=uchGet(generalbd%fname))

	return
	end subroutine parseInput
    subroutine controlVersions

    use approx,       only: axVersion,axDate,axAuthor
    use datablock,    only: bdVersion,bdDate,bdAuthor
    use derivat,      only: deVersion,deDate,deAuthor
    use fci,          only: ciVersion,ciDate,ciAuthor
    use fcontrol,     only: fcVersion,fcDate,fcAuthor
    use glob,         only: glVersion,glDate,glAuthor,iglu
    use hdb,          only: heVersion,heDate,heAuthor,ou
    use lrccsdModule, only: lrVersion,lrDate,lrAuthor
    use math,         only: maVersion,maDate,maAuthor
    use orientation,  only: orVersion,orDate,orAuthor
    use painter,      only: tdVersion,tdDate,tdAuthor
    use argsParser,   only: apVersion,apDate,apAuthor
    use printmod,     only: prVersion,prDate,prAuthor
    use txtParser,    only: tpVersion,tpDate,tpAuthor

    implicit none

    integer(kind=iglu) :: ln


    ln=max(&
            len(axAuthor),&
            len(bdAuthor),&
            len(deAuthor),&
            len(ciAuthor),&
            len(fcAuthor),&
            len(glAuthor),&
            len(heAuthor),&
            len(lrAuthor),&
            len(maAuthor),&
            len(orAuthor),&
            len(tdAuthor),&
            len(apAuthor),&
            len(prAuthor),&
            len(tpAuthor) &
    )

    write (ou,100) heVersion,heDate,adjustl(heAuthor),&
    &               glVersion,glDate,adjustl(glAuthor),&
    &               apVersion,apDate,adjustl(apAuthor),&
    &               maVersion,maDate,adjustl(maAuthor),&
    &               fcVersion,fcDate,adjustl(fcAuthor),&
    &               bdVersion,bdDate,adjustl(bdAuthor),&
    &               tdVersion,tdDate,adjustl(tdAuthor),&
    &               orVersion,orDate,adjustl(orAuthor),&
    &               prVersion,prDate,adjustl(prAuthor),&
    &              ciVersion,ciDate,adjustl(ciAuthor),&
    &              lrVersion,lrDate,adjustl(lrAuthor),&
    &              deVersion,deDate,adjustl(deAuthor),&
    &              axVersion,axDate,adjustl(axAuthor),&
    &              tpVersion,tpDate,adjustl(tpAuthor)

100 format (4X,<67+ln>('-')/&
&            4X,'|',1X,"Libraries' info",<49+ln>X,'|'/&
&            4X,<67+ln>('-')/&
&            4X,'|',1X,'Core                  :',1X,'v.',A5,1X,'Last modified at',1X,A10,1X,'by',1X,A<ln>,1X,'|'/&
&            4X,'|',1X,'Global                :',1X,'v.',A5,1X,'Last modified at',1X,A10,1X,'by',1X,A<ln>,1X,'|'/&
&            4X,'|',1X,'Argument parser       :',1X,'v.',A5,1X,'Last modified at',1X,A10,1X,'by',1X,A<ln>,1X,'|'/&
&            4X,'|',1X,'Internal math module  :',1X,'v.',A5,1X,'Last modified at',1X,A10,1X,'by',1X,A<ln>,1X,'|'/&
&            4X,'|',1X,'File control          :',1X,'v.',A5,1X,'Last modified at',1X,A10,1X,'by',1X,A<ln>,1X,'|'/&
&            4X,'|',1X,'Block data            :',1X,'v.',A5,1X,'Last modified at',1X,A10,1X,'by',1X,A<ln>,1X,'|'/&
&            4X,'|',1X,'Geometric module      :',1X,'v.',A5,1X,'Last modified at',1X,A10,1X,'by',1X,A<ln>,1X,'|'/&
&            4X,'|',1X,'Orientation control   :',1X,'v.',A5,1X,'Last modified at',1X,A10,1X,'by',1X,A<ln>,1X,'|'/&
&            4X,'|',1X,'Printer               :',1X,'v.',A5,1X,'Last modified at',1X,A10,1X,'by',1X,A<ln>,1X,'|'/&
&            4X,'|',1X,'FCI module            :',1X,'v.',A5,1X,'Last modified at',1X,A10,1X,'by',1X,A<ln>,1X,'|'/&
&            4X,'|',1X,'Linear response module:',1X,'v.',A5,1X,'Last modified at',1X,A10,1X,'by',1X,A<ln>,1X,'|'/&
&            4X,'|',1X,'Numerical derivatives :',1X,'v.',A5,1X,'Last modified at',1X,A10,1X,'by',1X,A<ln>,1X,'|'/&
&            4X,'|',1X,'Approximation         :',1X,'v.',A5,1X,'Last modified at',1X,A10,1X,'by',1X,A<ln>,1X,'|'/&
&            4X,'|',1X,'Text parser           :',1X,'v.',A5,1X,'Last modified at',1X,A10,1X,'by',1X,A<ln>,1X,'|'/&
&            4X,<67+ln>('-'))

    return
    end subroutine controlVersions
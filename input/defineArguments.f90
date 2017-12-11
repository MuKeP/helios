    subroutine defineArguments

    use glob,       only: void,true,false
    use hdb,        only: generalbd,showInputHelp
    use argsParser, only: addArg,apShare,addGroup,appendRule
    use txtParser,  only: tpTranslateEscapes,tpFill

    implicit none


    void=addArg('-f|--file','input file name',&
                var=apShare(generalbd%infile),expect=true,default=false,required=true ,group=addGroup('General'))

    void=addArg('-i|--help-input',tpTranslateEscapes('show help for blocks in input file.\n'//&
                tpFill(26)//'Type "all" to see the list.'),&
                var=apShare(showInputHelp),expect=true,default=false,required=true,group=addGroup('General'))

    void=appendRule(['--help-input','--file'],'xor')

    return
    end subroutine defineArguments
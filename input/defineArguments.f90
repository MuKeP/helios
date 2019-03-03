    subroutine defineArguments

    use glob,       only: void,true,false,iglu
    use hdb,        only: generalbd,showInputHelp,cparamstring
    use argsParser, only: addArg,apShare,addGroup,appendRule
    use txtParser,  only: tpTranslateEscapes,tpFill

    implicit none

    integer(kind=iglu) :: null

    void=addArg('-f|--file','input file name',&
                var=apShare(generalbd%infile),expect=true,default=false,required=true ,group=addGroup('General'))

    void=addArg('-i|--help-input',tpTranslateEscapes('show help for blocks in input file.\n'//&
                tpFill(26)//'Type "all" to see the list.'),&
                var=apShare(showInputHelp),expect=true,default=false,required=true,group=addGroup('General'))

    void=addArg('-t|--show-template','generate template for input file.',&
                var=apShare(null),expect=true,default=false,required=true,group=addGroup('General'))

    void=addArg('-p|--change-parameter','change parameter of input file. use bd1:param1=x1|bd2:param2=x2...',&
                var=apShare(cparamstring),expect=true,default=false,required=false,group=addGroup('General'))

    void=appendRule(['--help-input','--file', '--show-template'],'xor')

    return
    end subroutine defineArguments
    subroutine defineArguments

    use glob      , only: void,true,false
    use hdb       , only: generalbd,showInputHelp
    use argsParser, only: addArg,apShare,addGroup

    implicit none


    void=addArg('-f|--file'      ,'input file name'         ,var=apShare(generalbd%infile),expect=true,&
                default=false,required=true ,group=addGroup('General'))

    void=addArg('-i|--help-input','show help for input file',var=apShare(showInputHelp)     ,expect=false,&
                default=true ,required=false,group=addGroup('General'))

    return
    end subroutine defineArguments
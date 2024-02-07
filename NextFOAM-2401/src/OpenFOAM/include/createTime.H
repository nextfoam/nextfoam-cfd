Foam::Info<< "Create time\n" << Foam::endl;

Foam::Time runTime(Foam::Time::controlDictName, args);

#ifdef fvCFD_H
if (Pstream::master())
{
    args.createProcInfoFile();
}
#endif

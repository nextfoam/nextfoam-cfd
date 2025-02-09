/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION 

    PtrList<fv::options> solidFvOptions(solidRegions.size());
    forAll(solidRegions, i)
    {
        solidFvOptions.set(i, new fv::options(solidRegions[i]));

        if (!solidFvOptions[i].optionList::size())
        {
            Info << "No finite volume options present" << endl;
        }
    }

    PtrList<fv::options> fluidFvOptions(fluidRegions.size());
    forAll(fluidRegions, i)
    {
        Info<< "    Adding fvOptions\n" << endl;
        fluidFvOptions.set(i, new fv::options(fluidRegions[i]));

        if (!fluidFvOptions[i].optionList::size())
        {
            Info << "No finite volume options present" << endl;
        }
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

    fv::options& fvOptions(fv::options::New(mesh));

    if (!fvOptions.optionList::size())
    {
        Info << "No finite volume options present" << endl;
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

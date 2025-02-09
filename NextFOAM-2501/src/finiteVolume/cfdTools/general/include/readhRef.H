/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION 

    PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        Info<< "    Adding to hRefFluid\n" << endl;
        hRefFluid.set
        (
            i,
            new uniformDimensionedScalarField
            (
                IOobject
                (
                    "hRef",
                    runTime.constant(),
                    fluidRegions[i].thisDb(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                dimensionedScalar(word::null, dimLength, Zero)
            )
        );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

    Info<< "\nReading hRef" << endl;
    uniformDimensionedScalarField hRef
    (
        IOobject
        (
            "hRef",
            runTime.constant(),
            mesh.thisDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dimensionedScalar(word::null, dimLength, Zero)
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

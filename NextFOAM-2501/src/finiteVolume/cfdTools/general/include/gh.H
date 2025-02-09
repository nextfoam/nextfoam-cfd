/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION 

    PtrList<volScalarField> ghFluid(fluidRegions.size());
    PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
    PtrList<dimensionedScalar> ghRefFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        ghRefFluid.set
        (
            i,
            new dimensionedScalar
            (
                "ghRef",
                mag(g.value()) > SMALL
              ? g & (cmptMag(g.value())/mag(g.value()))*hRefFluid[i]
              : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
            )
        );

        const int oldLocal = volVectorField::Boundary::localConsistency;
        volVectorField::Boundary::localConsistency = 0;

        Info<< "    Adding to ghFluid\n" << endl;
        ghFluid.set
        (
            i,
            new volScalarField
            (
                "gh",
                (g & fluidRegions[i].C()) - ghRefFluid[i]
            )
        );

        Info<< "    Adding to ghfFluid\n" << endl;
        ghfFluid.set
        (
            i,
            new surfaceScalarField
            (
                "ghf",
                (g & fluidRegions[i].Cf()) - ghRefFluid[i]
            )
        );

        volVectorField::Boundary::localConsistency = oldLocal;
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

    Info<< "Calculating field g.h\n" << endl;
    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRef
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );
    const int oldLocal = volVectorField::Boundary::localConsistency;
    volVectorField::Boundary::localConsistency = 0;

    volScalarField gh("gh", (g & mesh.C()) - ghRef);
    surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);

    volVectorField::Boundary::localConsistency = oldLocal;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

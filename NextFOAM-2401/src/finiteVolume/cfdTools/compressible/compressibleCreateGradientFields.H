/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION 

    PtrList<volVectorField> gradp_rghFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        gradp_rghFluid.set
        (
            i,
            new volVectorField
            (
                IOobject("grad(p_rgh)", runTime.timeName(), fluidRegions[i]),
                fluidRegions[i],
                dimensionedVector("grad(p_rgh)", dimPressure/dimLength, Zero)
            )
        );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

    volVectorField gradp_rgh
    (
        IOobject("grad(p_rgh)", runTime.timeName(), mesh),
        mesh,
        dimensionedVector("grad(p_rgh)", dimPressure/dimLength, Zero)
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

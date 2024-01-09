/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION 

    PtrList<volScalarField> AUFluid(fluidRegions.size());
    PtrList<volScalarField> rAUFluid(fluidRegions.size());
    PtrList<tmp<volScalarField>> drAUFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        AUFluid.set
        (
            i,
            new volScalarField
            (
                IOobject("UAp", runTime.timeName(), fluidRegions[i]),
                fluidRegions[i],
                dimensionedScalar("UAp", dimMass/dimTime, One)
            )
        );

        rAUFluid.set
        (
            i,
            new volScalarField
            (
                IOobject("rUAp", runTime.timeName(), fluidRegions[i]),
                fluidRegions[i],
                dimensionedScalar("rUAp", dimTime/dimDensity, One)
            )
        );

        drAUFluid.set(i, new tmp<volScalarField>(rAUFluid[i]));
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

    volScalarField AU
    (
        IOobject("UAp", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("UAp", dimMass/dimTime, One)
    );

    volScalarField rAU
    (
        IOobject("rUAp", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("rUAp", dimTime/dimDensity, One)
    );

    tmp<volScalarField> drAU(rAU);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION 

    PtrList<const dimensionedScalar> initialMassFluid(fluidRegions.size());
    List<scalar> pCorrLimitFluid(fluidRegions.size());
    PtrList<dimensionedScalar> pCorrMaxFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        uniformDimensionedScalarField initialMassIO
        (
            IOobject
            (
                "initialMass",
                runTime.timeName(),
                "uniform",
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::domainIntegrate(rhoFluid[i])
        );

        initialMassFluid.set
        (
            i,
            new const dimensionedScalar
            (
                "initialMass",
                dimMass,
                initialMassIO.value()
            )
        );

        pCorrLimitFluid[i] =
            pvCoupling.dict(fluidRegions.get(i)).getOrDefault<scalar>
            (
                "pCorrLimit", 
                GREAT
            );

        pCorrMaxFluid.set
        (
            i,
            new dimensionedScalar
            (
                pCorrLimitFluid[i]*max(thermoFluid[i].p().internalField())
            )
        );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

    uniformDimensionedScalarField initialMassIO
    (
        IOobject
        (
            "initialMass",
            runTime.timeName(),
            "uniform",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::domainIntegrate(rho)
    );

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassIO.value()
    );

    scalar pCorrLimit
    (
        pvCoupling.dict().getOrDefault<scalar>("pCorrLimit", GREAT)
    );

    dimensionedScalar pCorrMax(pCorrLimit*max(p.internalField()));

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

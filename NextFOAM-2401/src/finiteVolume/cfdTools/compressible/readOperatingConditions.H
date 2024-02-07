/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION 

    PtrList<IOdictionary> operatingConditionsFluid(fluidRegions.size());
    PtrList<dimensionedScalar> operatingPressureFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        operatingConditionsFluid.set
        (
            i,  
            new IOdictionary
            (
                IOobject
                (
                    "operatingConditions",
                    runTime.constant(),
                    fluidRegions[i],
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE
                )
            )
        );

        operatingPressureFluid.set
        (
            i,
            new dimensionedScalar
            (
                dimensionedScalar::getOrDefault
                (
                    "operatingPressure",
                    operatingConditionsFluid[i],
                    dimPressure,
                    101325
                )
            )
        );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

    IOdictionary operatingConditions
    (
        IOobject
        (
            "operatingConditions",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar operatingPressure
    (
        dimensionedScalar::getOrDefault
        (
            "operatingPressure",
            operatingConditions,
            dimPressure,
            101325
        )
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

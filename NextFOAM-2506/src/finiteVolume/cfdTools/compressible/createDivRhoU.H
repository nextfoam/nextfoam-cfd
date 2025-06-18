/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION 

    PtrList<autoPtr<volScalarField>> divrhoUFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        divrhoUFluid.set(i, new autoPtr<volScalarField>);

        if (correctPhiFluid[i])
        {
            divrhoUFluid[i].reset
            (
                new volScalarField
                (
                    "divrhoU",
                    fvc::div
                    (   
                        fvc::absolute(phiFluid[i], rhoFluid[i], UFluid[i])
                    )
                )
            );
        }
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

    autoPtr<volScalarField> divrhoU;
    if (correctPhi)
    {
        divrhoU.reset
        (
            new volScalarField
            (
                "divrhoU",
                fvc::div(fvc::absolute(phi, rho, U))
            )
        );
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

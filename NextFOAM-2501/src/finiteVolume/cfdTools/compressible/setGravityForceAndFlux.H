bool gravityEnabled(mag(g.value()) > SMALL ? true : false);

/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION 

    PtrList<volScalarField> hydroStaticPressureFluid(fluidRegions.size());
    #if defined(PIMPLE_CONTROL)
    PtrList<volScalarField> hydroStaticDensityFluid(fluidRegions.size());
    #endif
    PtrList<surfaceScalarField> gravityFluxPotentialFluid(fluidRegions.size());
    PtrList<volVectorField> gravityForceFluid(fluidRegions.size());
    PtrList<surfaceScalarField> gravityFluxFluid(fluidRegions.size());

    forAll(fluidRegions, i)
    {
        hydroStaticPressureFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "hydroStaticPressure", 
                    runTime.timeName(), 
                    fluidRegions[i]
                ),
                fluidRegions[i],
                dimensionedScalar("hydroStaticPressure", dimPressure, Zero)
            )
        );

        #if defined(PIMPLE_CONTROL)
        hydroStaticDensityFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "hydroStaticDensity", 
                    runTime.timeName(), 
                    fluidRegions[i]
                ),
                fluidRegions[i],
                dimensionedScalar("hydroStaticDensity", dimDensity, Zero)
            )
        );
        #endif

        gravityFluxPotentialFluid.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "gravityFluxPotential",
                    runTime.timeName(),
                    fluidRegions[i]
                ),
                fluidRegions[i],
                dimensionedScalar
                (
                    "gravityFluxPotential",
                    dimForce/dimLength,
                    Zero
                )
            )
        );

        gravityForceFluid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "gravityForce",
                    runTime.timeName(),
                    fluidRegions[i]
                ),
                fluidRegions[i],
                dimensionedVector
                (
                    "gravityForce",
                    dimDensity*dimAcceleration,
                    Zero
                )
            )
        );

        gravityFluxFluid.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "gravityFlux",
                    runTime.timeName(),
                    fluidRegions[i]
                ),
                fluidRegions[i],
                dimensionedScalar
                (
                    "gravityFlux",
                    dimDensity*dimVelocity*dimArea,
                    Zero
                )
            )
        );

        if (gravityEnabled)
        {
            hydroStaticPressureFluid[i] = rhoFluid[i]*ghFluid[i];
            #if defined(PIMPLE_CONTROL)
            hydroStaticDensityFluid[i] = 
                thermoFluid[i].psi()
                *(hydroStaticPressureFluid[i] + operatingPressureFluid[i]);
            #endif
            gravityForceFluid[i] = ghFluid[i]*fvc::grad(rhoFluid[i]);
            gravityFluxPotentialFluid[i] = 
                ghfFluid[i]*fvc::snGrad(rhoFluid[i])*fluidRegions[i].magSf();
            //gravityForceFluid[i] = fvc::reconstruct(gravityFluxPotentialFluid[i]);
        }
        else
        {
            hydroStaticPressureFluid[i].clear();
            #if defined(PIMPLE_CONTROL)
            hydroStaticDensityFluid[i].clear();
            #endif
            gravityFluxPotentialFluid[i].clear();
            gravityForceFluid[i].clear();
            gravityFluxFluid[i].clear();
        }
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

    tmp<volScalarField> tHydroStaticPressure
    (
        new volScalarField
        (
            IOobject("hydroStaticPressure", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar("hydroStaticPressure", dimPressure, Zero)
        )
    );
    volScalarField& hydroStaticPressure = tHydroStaticPressure.ref();

    #if defined(PIMPLE_CONTROL)
    tmp<volScalarField> tHydroStaticDensity
    (
        new volScalarField
        (
            IOobject("hydroStaticDensity", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar("hydroStaticDensity", dimDensity, Zero)
        )
    );
    volScalarField& hydroStaticDensity = tHydroStaticDensity.ref();
    #endif

    tmp<surfaceScalarField> tGravityFluxPotential
    (
        new surfaceScalarField
        (
            IOobject("gravityFluxPotential", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar("gravityFluxPotential", dimForce/dimLength, Zero)
        )
    );
    surfaceScalarField& gravityFluxPotential = tGravityFluxPotential.ref();

    tmp<volVectorField> tGravityForce
    (
        new volVectorField
        (
            IOobject("gravityForce", runTime.timeName(), mesh),
            mesh,
            dimensionedVector("gravityForce", dimDensity*dimAcceleration, Zero)
        )
    );
    volVectorField& gravityForce = tGravityForce.ref();

    tmp<surfaceScalarField> tGravityFlux
    (
        new surfaceScalarField
        (
            IOobject("gravityFlux", runTime.timeName(), mesh),
            mesh,
            dimensionedScalar
            (
                "gravityFlux", 
                dimDensity*dimVelocity*dimArea, 
                Zero
            )
        )
    );
    #ifdef NO_POST_PROCESS
    surfaceScalarField& gravityFlux = tGravityFlux.ref();
    #endif

    if (gravityEnabled)
    {
        hydroStaticPressure = rho*gh;
        #if defined(PIMPLE_CONTROL)
        hydroStaticDensity = psi*(hydroStaticPressure + operatingPressure);
        #endif
        gravityForce = gh*fvc::grad(rho);
        gravityFluxPotential = ghf*fvc::snGrad(rho)*mesh.magSf();
        //gravityForce = fvc::reconstruct(gravityFluxPotential)
    }
    else
    {
        tHydroStaticPressure.clear();
        #if defined(PIMPLE_CONTROL)
        tHydroStaticDensity.clear();
        #endif
        tGravityFluxPotential.clear();
        tGravityForce.clear();
        tGravityFlux.clear();
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

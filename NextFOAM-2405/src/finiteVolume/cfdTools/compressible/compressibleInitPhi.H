/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION 

    Info<< "    Adding to phiFluid\n" << endl;

    tmp<surfaceScalarField> initPhi
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fluidRegions[i],
            dimensionedScalar(dimMass/dimTime, Zero)
        )
    );

    initPhi.ref().setOriented();

    //initPhi.ref() = 
    //    linearInterpolate(rhoFluid[i]*UFluid[i]) & fluidRegions[i].Sf();
    //initPhi.ref() = 
    //    midPointInterpolate(rhoFluid[i]*UFluid[i]) & fluidRegions[i].Sf();

    //- determine if the case is inviscid
    if (gAverage(thermoFluid[i].nu()().primitiveField()) < SMALL)
    {
        initPhi.ref() = 
            linearInterpolate(rhoFluid[i]*UFluid[i]) & fluidRegions[i].Sf();
    }

    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            initPhi
        )
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

    Info<< "Reading/calculating face flux field phi\n" << endl;

    tmp<surfaceScalarField> initPhi
    (
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(dimMass/dimTime, Zero)
        )
    );

    initPhi.ref().setOriented();

    //initPhi.ref() = linearInterpolate(rho*U) & mesh.Sf();
    //initPhi.ref() = midPointInterpolate(rho*U) & mesh.Sf();

    //- determine if the case is inviscid
    if (gAverage(thermo.nu()().primitiveField()) < SMALL)
    {
        initPhi.ref() = linearInterpolate(rho*U) & mesh.Sf();
    }

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        initPhi
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

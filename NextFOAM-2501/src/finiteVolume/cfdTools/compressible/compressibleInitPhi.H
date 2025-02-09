/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION 

    Info<< "    Adding to phiFluid\n" << endl;

    tmp<volVectorField> tInitU(new volVectorField(UFluid[i]));

    volVectorField& initU = tInitU.ref();

    forAll(initU.boundaryField(), patchi)
    {
        if (!isA<wallFvPatch>(fluidRegions[i].boundary()[patchi]))
        {
            initU.boundaryFieldRef()[patchi] ==
                initU.boundaryField()[patchi].patchInternalField();
        }
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
            linearInterpolate(rhoFluid[i]*initU) & fluidRegions[i].Sf()
        )
    );

    tInitU.clear();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

    Info<< "Reading/calculating face flux field phi\n" << endl;

    tmp<volVectorField> tInitU(new volVectorField(U));

    volVectorField& initU = tInitU.ref();

    forAll(initU.boundaryField(), patchi)
    {
        if (!isA<wallFvPatch>(mesh.boundary()[patchi]))
        {
            initU.boundaryFieldRef()[patchi] ==
                initU.boundaryField()[patchi].patchInternalField();
        }
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
        linearInterpolate(rho*initU) & mesh.Sf()
    );

    tInitU.clear();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

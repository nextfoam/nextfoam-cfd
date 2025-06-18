/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION 

bool LTS(true);

forAll(fluidRegions, i)
{
    if (!fv::localEulerDdt::enabled(fluidRegions[i]))
    {
        LTS = false;
    }
}

PtrList<tmp<volScalarField>> trDeltaTFluid(fluidRegions.size());

if (LTS)
{
    Info<< "Using LTS" << endl;

    forAll(fluidRegions, i)
    {
        trDeltaTFluid.set
        (
            i,
            new tmp<volScalarField>
            (
                tmp<volScalarField>::New
                (
                    IOobject
                    (
                        fv::localEulerDdt::rDeltaTName,
                        runTime.timeName(),
                        fluidRegions[i],
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    fluidRegions[i],
                    dimensionedScalar("one", dimless/dimTime, 1),
                    extrapolatedCalculatedFvPatchScalarField::typeName
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

bool LTS = fv::localEulerDdt::enabled(mesh);

tmp<volScalarField> trDeltaT;

if (LTS)
{
    Info<< "Using LTS" << endl;

    trDeltaT = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                fv::localEulerDdt::rDeltaTName,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("one", dimless/dimTime, 1),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

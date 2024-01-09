/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION

    PtrList<autoPtr<CombustionModel<rhoReactionThermo>>> reactionFluid
    (
        fluidRegions.size()
    );

    PtrList<multivariateSurfaceInterpolationScheme<scalar>::fieldTable>
        fieldsFluid(fluidRegions.size());

    PtrList<volScalarField> QdotFluid(fluidRegions.size());

    List<bool> reactionEnabledFluid(fluidRegions.size(), false);

    forAll(fluidRegions, i)
    {
        Info<< "    Adding to reactionFluid\n" << endl;
        reactionFluid.set
        (
            i,
            new autoPtr<CombustionModel<rhoReactionThermo>>
            (
                CombustionModel<rhoReactionThermo>::New
                (
                    thermoFluid[i],
                    turbulenceFluid[i]
                )
            )
        );

        Info<< "    Adding to fieldsFluid\n" << endl;
        fieldsFluid.set
        (
            i,
            new multivariateSurfaceInterpolationScheme<scalar>::fieldTable
        );
        forAll(thermoFluid[i].composition().Y(), j)
        {
            fieldsFluid[i].add(thermoFluid[i].composition().Y()[j]);
        }
        fieldsFluid[i].add(thermoFluid[i].he());

        Info<< "    Adding to QdotFluid\n" << endl;
        QdotFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Qdot",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
            )
        );

        #ifdef NO_POST_PROCESS
        if 
        (
            thermoFluid[i].composition().Y().size() 
         && reactionFluid[i]->active()
        )
        {
            reactionEnabledFluid[i] = true;
        }
        #endif
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

    basicSpecieMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();

    label inertIndex = -1;
    if (Y.size())
    {
        const word inertSpecie(thermo.get<word>("inertSpecie"));
        inertIndex = composition.species().find(inertSpecie);
        if (inertIndex < 0)
        {
            FatalIOErrorIn(args.executable().c_str(), thermo)
                << "Inert specie " << inertSpecie
                << " not found in available species "
                << flatOutput(composition.species())
                << exit(FatalIOError);
        }
    }

    Info<< "Creating reaction model\n" << endl;
    autoPtr<CombustionModel<rhoReactionThermo>> reaction
    (
        CombustionModel<rhoReactionThermo>::New(thermo, turbulence())
    );

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    forAll(Y, i)
    {
        fields.add(Y[i]);
    }
    fields.add(thermo.he());

    volScalarField Qdot
    (
        IOobject
        (
            "Qdot",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
    );

    #ifdef NO_POST_PROCESS
    bool reactionEnabled(Y.size() && reaction->active());
    #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

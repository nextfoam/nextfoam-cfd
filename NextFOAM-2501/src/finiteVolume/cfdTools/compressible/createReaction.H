/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION

    PtrList<autoPtr<CombustionModel<rhoReactionThermo>>> reactionFluid
    (
        fluidRegions.size()
    );

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
    autoPtr<volScalarField> Dt(nullptr);
    PtrList<dimensionedScalar> Dm(Y.size());

    autoPtr<volScalarField> speciesDiffusionEnthalpy(nullptr);

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

    autoPtr<volScalarField> Qdot(nullptr);

    #ifdef NO_POST_PROCESS
    bool reactionEnabled(Y.size() && reaction->active());
    #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

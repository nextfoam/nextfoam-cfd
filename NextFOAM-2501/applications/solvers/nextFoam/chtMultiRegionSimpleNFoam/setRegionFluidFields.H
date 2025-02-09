    label i(pvCoupling.fluidIndex());

    #if defined(PIMPLE_CONTROL)
    dynamicFvMesh& mesh = fluidRegions[i];
    #else
    fvMesh& mesh = fluidRegions[i];
    #endif

    volScalarField& gh = ghFluid[i];
    surfaceScalarField& ghf = ghfFluid[i];

    dimensionedScalar& operatingPressure = operatingPressureFluid[i];

    volVectorField& gradp_rgh = gradp_rghFluid[i];

    volScalarField& rAU = rAUFluid[i];
    tmp<volScalarField>& drAU = drAUFluid[i];

    rhoReactionThermo& thermo = thermoFluid[i];

    volScalarField& rho = rhoFluid[i];
    volScalarField& p = thermo.p();

    #if defined(PIMPLE_CONTROL)
    autoPtr<volScalarField>& divrhoU = divrhoUFluid[i];
    autoPtr<volVectorField>& rhoU = rhoUFluid[i];
    dimensionedScalar& ghRef = ghRefFluid[i];
    const volScalarField& psi = thermo.psi();
    volScalarField& K = KFluid[i];
    //pressureControl& pressureControl = pressureControlFluid[i];
    autoPtr<surfaceVectorField>& rhoUf = rhoUfFluid[i];
    bool& correctPhi = correctPhiFluid[i];
    bool& moveMeshOuterCorrectors = moveMeshOuterCorrectorsFluid[i];
    bool& checkMeshCourantNo = checkMeshCourantNoFluid[i];
    #endif

    volVectorField& U = UFluid[i];
    volScalarField& Umag = UmagFluid[i];
    surfaceScalarField& phi = phiFluid[i];
    surfaceScalarField& rhof = rhofFluid[i];

    autoPtr<compressible::turbulenceModel>& turbulence = turbulenceFluid[i];

    volScalarField& p_rgh = p_rghFluid[i];

    label& pRefCell = pRefCellFluid[i];
    scalar& pRefValue = pRefValueFluid[i];

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

    autoPtr<CombustionModel<rhoReactionThermo>>& reaction = reactionFluid[i];

    autoPtr<volScalarField> Qdot(nullptr);

    bool& reactionEnabled = reactionEnabledFluid[i];

    volScalarField& hydroStaticPressure
    (
        gravityEnabled
      ? hydroStaticPressureFluid[i]
      : p // dummy placeholder
    );

    #if defined(PIMPLE_CONTROL)
    volScalarField& hydroStaticDensity
    (
        gravityEnabled
      ? hydroStaticDensityFluid[i]
      : rho // dummy placeholder
    );
    #endif

    surfaceScalarField& gravityFluxPotential
    (
        gravityEnabled
      ? gravityFluxPotentialFluid[i]
      : phi // dummy placeholder
    );

    volVectorField& gravityForce
    (
        gravityEnabled
      ? gravityForceFluid[i]
      : U // dummy placeholder
    );

    surfaceScalarField& gravityFlux
    (
        gravityEnabled
      ? gravityFluxFluid[i]
      : phi // dummy placeholder
    );

    IOMRFZoneList& MRF = MRFFluid[i];
    bool& MRFEnabled = MRFEnabledFluid[i];

    autoPtr<radiation::radiationModel>& radiation = radiationFluid[i];
    bool& radiationEnabled = radiationEnabledFluid[i];

    fv::options& fvOptions = fluidFvOptions[i];

    const dimensionedScalar& initialMass = initialMassFluid[i];

    scalar& pCorrLimit = pCorrLimitFluid[i];

    dimensionedScalar& pCorrMax = pCorrMaxFluid[i];

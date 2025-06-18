{
    const dimensionedScalar rhoMax
    (
        "rhoMax",
        dimDensity,
        GREAT,
        pvCoupling.dict()
    );

    const dimensionedScalar rhoMin
    (
        "rhoMin",
        dimDensity,
        SMALL,
        pvCoupling.dict()
    );

    thermo.limitRho(rhoMin, rhoMax);

    rho = thermo.rho();

    rho.relax();

    // A temporary, unjustified way to stabilize mass flow rate. 
    // Needs to be implemented as a method with a good rationale.
    {
        rhof.storePrevIter();

        scalar UUrf = mesh.equationRelaxationFactor(U.name());
        scalar pUrf = mesh.fieldRelaxationFactor(p_rgh.name());
        scalar hUrf = mesh.equationRelaxationFactor(thermo.he().name());
        scalar rhofUrf = sqr(UUrf)*pUrf*hUrf;

        rhof = fvc::interpolate(rho);
        rhof.relax(rhofUrf);
    }

    if (gravityEnabled)
    {
        hydroStaticPressure = rho*gh;
        gravityForce = gh*fvc::grad(rho);
        //gravityForce = fvc::reconstruct(gravityFluxPotential);
        gravityFluxPotential = ghf*fvc::snGrad(rho)*mesh.magSf();
    }

    //----------------------------------------------------------------------//
    // Treatment to maintain initial mass when simulating closed domain
    // problems. ---> What if the mass source exists??
    bool maintainInitialMass
    (
        pvCoupling.dict().lookupOrDefault<bool>("maintainInitialMass", false)
    );

    if (p_rgh.needReference() && maintainInitialMass)
    {
        // Changing the density doesn't really affect the mass conservation
        //, but explicitly changing the pressure solution never causes 
        // the pressure equation to converge.

        const volScalarField& psi = thermo.psi();
            
        dimensionedScalar compressibility = fvc::domainIntegrate(psi);

        if (compressibility.value() > SMALL)
        {
            dimensionedScalar totalMass(fvc::domainIntegrate(rho));
            //rho += psi*(initialMass - totalMass)/compressibility;
            p += (initialMass - totalMass)/compressibility;
            p_rgh = p - operatingPressure - hydroStaticPressure;
            p_rgh.correctBoundaryConditions();

            Info << "Mass Imbalance = " 
                << (initialMass - totalMass).value()
                << nl << endl;
        }
    }
    //----------------------------------------------------------------------//
}

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

    #include "compressibleContinuityErrs.H"

    rho = thermo.rho();

    rho.relax();

    rhof = fvc::interpolate(rho);

    if (gravityEnabled)
    {
        hydroStaticPressure = rho*gh;
        hydroStaticDensity = psi*(hydroStaticPressure + operatingPressure);
        gravityForce = gh*fvc::grad(rho);
        gravityFluxPotential = ghf*fvc::snGrad(rho)*mesh.magSf();
        //gravityForce = fvc::reconstruct(gravityFluxPotential);
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

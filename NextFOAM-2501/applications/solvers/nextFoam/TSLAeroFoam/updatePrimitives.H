    // bound density
    //boundMinMax(rho,rhoMin,rhoMax);
    bound(rho, rhoMin);

    // bound rhoE
    //boundMinMax(rhoE,rhoEMin,rhoEMax);
    bound(rhoE, rhoEMin);

    // Compute internal field of U
    U.primitiveFieldRef() = rhoU.primitiveField()/rho.primitiveField();

    // Compute internal field of static enthalpy
    he.primitiveFieldRef() =
    (
        gamma.primitiveField()*(rhoE.primitiveField()/rho.primitiveField()
      - 0.5*magSqr(U.primitiveField()))
      //- turbulence->k()().primitiveField() //cwo
    );

    bound(he, hMin);

    // Compute internal field of p
    p.primitiveFieldRef() = 
    (
        (1.0 - 1.0/gamma.primitiveField())*rho.primitiveField()
        *he.primitiveField()
    );

    // Correct internal field of T, psi, mu and alpha
    thermo.correctInternal();

    // Update boundary conditions of U
    U.correctBoundaryConditions();

    // Bound the velocity
    volScalarField magU(mag(U));

    if (max(magU) > UMax)
    {
        Info<< "bounding " << U.name()
            << " max: " << max(magU).value()
            << endl;

        volScalarField Ulimiter
        (
            pos(magU - UMax)*UMax/(magU + smallU) + neg(magU - UMax)
        );

        Ulimiter.max(scalar(0));
        Ulimiter.min(scalar(1));

        U *= Ulimiter;
        U.correctBoundaryConditions();
    }

    // Correct boundary conditions of static enthalpy
    he.correctBoundaryConditions();
    //boundMinMax(h,hMin,hMax);
    
    // Correct boundary conditions of p
    p.correctBoundaryConditions();
    //boundMinMax(p,pMin,pMax);

    // Correct boundary field of psi, mu and alpha
    thermo.correctBoundary();

    // Update boundary field of rho
    rho.boundaryFieldRef() == p.boundaryField()*thermo.psi().boundaryField();

    // Update boundary field of rhoU
    rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

    // Update boundary field of rhoE
    rhoE.boundaryFieldRef() ==
    (
        (
            0.5*magSqr(U.boundaryField())
          //+ turbulence->k()().boundaryField() //cwo
        )*rho.boundaryField() 
      + p.boundaryField()/(gamma.boundaryField()-1.0)
    );  

    // Needed for turbulence and CoEuler ddt scheme and maybe elsewhere;
    phi = Godunov->rhoFlux();

    // Update mach number for post-processing
    Mach = mag(U)/sqrt(thermo.gamma()*p/rho);


    #include "createSnGradCorrLimiter.H" // for snGradCorr limiting

    Info<< "Reading thermophysical properties\n" << endl;
    
    // In etc/controlDict, Tstd is defined as 298.15
    setTstd(0.0);
    
    autoPtr<fluidThermo> pThermo
    (
        fluidThermo::New(mesh)
    );
    fluidThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h", "e");

    volScalarField& he = thermo.he();
    volScalarField& p = thermo.p();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );
    rho.oldTime(); //by Gill : for allocation of field0Ptr pointer

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "compressibleCreatePhi.H"

    Info<< "Creating turbulence model" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        ));

    Info<< nl << "Creating Riemann solver" << endl;
    autoPtr<compressible::godunovFlux> Godunov
    (
        compressible::godunovFlux::New
        (
            rho,
            U,
            thermo, 
            turbulence()
        )
    );

    Info<< nl << "Creating Spectral radius" << nl << endl;
    compressible::spectralRadius spectralRadius
    (
        rho, 
        U,
        thermo, 
        turbulence()
    );

    volVectorField rhoU
    (
        IOobject
        (
            "rhoU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            //IOobject::AUTO_WRITE
            IOobject::NO_WRITE
        ),
        rho*U
    );
    rhoU.oldTime(); //by Gill : for allocation of field0Ptr pointer

    volScalarField rhoE
    (
        IOobject
        (
            "rhoE",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            //IOobject::AUTO_WRITE
            IOobject::NO_WRITE
        ),
        //rho*(h + 0.5*magSqr(U) + turbulence->k()) - p //cwo
        rho*(he + 0.5*magSqr(U)) - p
    );
    rhoE.oldTime(); //by Gill : for allocation of field0Ptr pointer

    // Mach number field for post process
    volScalarField Mach
    (
        IOobject
        (
            "Mach",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mag(U)/sqrt(thermo.gamma()*p/rho)
    );

    compressible::localTimeStep localTimeStep(U, turbulence());

    IOdictionary RMSResDict
    (
        IOobject
        (
            "RMSResidual",
            runTime.timeName(),
            "uniform",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE,
            true
        )
    );

    scalar RMSResRhoMax5 =    
        RMSResDict.lookupOrAddDefault("RMSResRhoMax5", 1e-16);

    vector RMSResRhoUMax5 =
        RMSResDict.lookupOrAddDefault
        (
            "RMSResRhoUMax5", 
            vector(1e-16, 1e-16, 1e-16)
        );

    scalar RMSResRhoEMax5 =
        RMSResDict.lookupOrAddDefault("RMSResRhoEMax5", 1e-16);

    scalar RMSResRhoOld =    
        RMSResDict.lookupOrAddDefault("RMSResRhoOld", 1.0);

    vector RMSResRhoUOld =
        RMSResDict.lookupOrAddDefault
        (
            "RMSResRhoUOld", 
            vector(1.0, 1.0, 1.0)
        );

    scalar RMSResRhoEOld =
        RMSResDict.lookupOrAddDefault("RMSResRhoEOld", 1.0);


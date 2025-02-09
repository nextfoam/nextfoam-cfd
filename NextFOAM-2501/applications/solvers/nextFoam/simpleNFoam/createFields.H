// ------------------------------------------------------------------------- //

#include "createSnGradCorrLimiter.H" // for snGradCorr limiting

volVectorField gradp
(
    IOobject("grad(p)", runTime.timeName(), mesh),
    mesh,
    dimensionedVector("grad(p)", dimAcceleration, Zero)
);

volScalarField AU
(
    IOobject("UAp", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("UAp", dimVolume/dimTime, One)
);

volScalarField rAU
(
    IOobject("rUAp", runTime.timeName(), mesh),
    mesh,
    dimensionedScalar("rUAp", dimTime, One)
);

tmp<volScalarField> drAU(rAU);

// ------------------------------------------------------------------------- //

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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

volScalarField Umag
(
    IOobject
    (
        "Umag",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(U)
);

#include "initPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pvCoupling.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"

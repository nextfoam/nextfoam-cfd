#ifdef PIMPLE_CONTROL
if (!thermo.isochoric())
{
    label debug(Foam::solverPerformance::debug);
    Foam::solverPerformance::debug = 0;
    #include "rhoEqn.H"
    Foam::solverPerformance::debug = debug;
}
#endif

// Update transport properties for porous media including porosity effect
volSymmTensorField alphaEff
(
    "alphaEff",
    turbulence->alphaEff()*symmTensor::I 
);

volSymmTensorField kappaEff
(
    "kappaEff",
    turbulence->kappaEff()*symmTensor::I 
);

fvOptions.adjustTransport(kappaEff, alphaEff, thermo.Cp());

volScalarField& he = thermo.he();

tmp<fvScalarMatrix> diffusionEEqn
(
    Y.size()
  ? 
    correction(fvm::laplacian(alphaEff, he)) 
  + fvc::laplacian(kappaEff, thermo.T())
  + speciesDiffusionEnthalpy()
  : 
    fvm::laplacian(alphaEff, he)
);

fvScalarMatrix EEqn(fvm::div(phi, he) - diffusionEEqn == fvOptions(rho, he));

#ifdef PIMPLE_CONTROL
EEqn += fvm::ddt(rho, he);
#endif

if (thermo.includeViscousDissipation())
{
    EEqn += fvc::div(turbulence->devRhoReff() & U);
}

if (thermo.includeKineticEnergy() || !thermo.incompressible())
{
    #ifdef PIMPLE_CONTROL
    K = 0.5*magSqr(U);
    EEqn += fvc::ddt(rho, K);
    EEqn += fvc::div(phi, K, "div(phi," + he.name() + ")");
    #else
    EEqn += fvc::div(phi, 0.5*magSqr(U), "div(phi," + he.name() + ")");
    #endif
}

#include "pressureWork.H"

if (reactionEnabled)
{
    EEqn -= Qdot();
}

if (gravityEnabled)
{
    EEqn -= rho*(U&g);
}

if (radiationEnabled)
{
    EEqn -= radiation->Sh(thermo, he);
}

EEqn.relax();

fvOptions.constrain(EEqn);

EEqn.solve();

fvOptions.correct(he);

thermo.correct();

if (radiationEnabled)
{
    radiation->correct();
}

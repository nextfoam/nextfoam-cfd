K = 0.5*magSqr(U);

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);

    if (mesh.moving())
    {
        dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
    }
}

if (!thermo.isochoric())
{
    label debug(Foam::solverPerformance::debug);
    Foam::solverPerformance::debug = 0;
    #include "rhoEqn.H"
    Foam::solverPerformance::debug = debug;
}

volScalarField& he = thermo.he();

fvScalarMatrix EEqn
(
    fvm::ddt(rho, he) + fvm::div(phi, he)
  + fvc::ddt(rho, K) + fvc::div(phi, K, "div(phi," + he.name() + ")")
  + (
        he.name() == "e"
      ? fvc::div(fvc::absolute(phi, U), p/rho, "div(phi,e)")
      : -dpdt
    )
  - fvm::laplacian(turbulence->alphaEff(), he)
 ==
    fvOptions(rho, he)
);

if (reactionEnabled)
{
    EEqn -= Qdot;
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

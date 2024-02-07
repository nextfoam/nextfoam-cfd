volScalarField& he = thermo.he();

fvScalarMatrix EEqn
(
    fvm::div(phi, he)
  + (
        he.name() == "e"
      ? fvc::div(phi, 0.5*magSqr(U) + p/rho, "div(phi,e)")
      : fvc::div(phi, 0.5*magSqr(U), "div(phi,h)")
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

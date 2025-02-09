// Solve the Momentum equation
if (MRFEnabled)
{
    MRF.correctBoundaryVelocity(U);
}

tmp<fvVectorMatrix> tUEqn
(
    fvm::div(phi, U)
  + turbulence->divDevRhoReff(U)
 ==
    fvOptions(rho, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.updateCoeffs();

#ifdef PIMPLE_CONTROL
UEqn += fvm::ddt(rho, U);
#endif

UEqn.relax();

fvOptions.constrain(UEqn);

// Update pressure gradient
if (pvCoupling.momentumPredictor() || pvCoupling.consistent())
{
    gradp_rgh = fvc::gaussGrad(p_rgh);
}

if (pvCoupling.momentumPredictor())
{
    volVectorField UEqnSource(gradp_rgh);

    if (gravityEnabled)
    {
        UEqnSource += gravityForce;
    }

    if (MRFEnabled)
    {
        UEqnSource += MRF.coriolisForce(rho, U);
    }

    solve(UEqn == -UEqnSource);

    fvOptions.correct(U);
}

rAU = UEqn.VbyA(pvCoupling.consistent() || fvOptions.hasActivePorousZone());

tmp<volScalarField> rAtU(rAU);

if (pvCoupling.consistent())
{
    const volScalarField& AbyV(UEqn.AbyV());
    rAtU = 1.0/max(AbyV - UEqn.H1(), 0.1*AbyV);
    drAU = rAtU() - rAU;
}

surfaceScalarField rhorAUf("rhorAUf", rhof*(momentumInterpolate(rAtU())));

if (gravityEnabled)
{
    gravityFlux = -rhorAUf*gravityFluxPotential;
}

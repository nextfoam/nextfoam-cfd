if (MRFEnabled)
{
    MRF.updateForce(phi, U);
}

volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

tUEqn.clear();

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (momentumInterpolate(HbyA) & mesh.Sf())
  + fvc::relaxCorr(phi, U, p, U.name())
);

#include "setConsistentFlux.H"

// Add Coriolis flux and make it relative to MRF
if (MRFEnabled)
{
    MRF.addCoriolisFlux(rAUf, phiHbyA);
}

adjustPhi(phiHbyA, U, p);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAtU(), MRF);

// Non-orthogonal pressure corrector loop
while (pvCoupling.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve();

    if (pvCoupling.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

// Momentum corrector
#include "correctVelocity.H"
fvOptions.correct(U);

// Explicitly relax pressure for next momentum predictor
p.relax();
p.correctBoundaryConditions();

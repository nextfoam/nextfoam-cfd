if (MRFEnabled)
{
    MRF.updateForce(rhof, phi, U);
}

volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

tUEqn.clear();

surfaceScalarField phiHbyA
(
    "phiHbyA",
    rhof*(momentumInterpolate(HbyA) & mesh.Sf())
  + fvc::relaxCorr(rhof, phi, U, p, U.name())
);

#include "compressibleSetConsistentFlux.H"

// Add Coriolis flux and make it relative to MRF
if (MRFEnabled)
{
    MRF.addCoriolisFlux(rhof, rhorAUf, phiHbyA);
}

adjustPhi(phiHbyA, U, p_rgh);

if (gravityEnabled)
{
    phiHbyA += constrainPhig(gravityFlux, U);
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

while (pvCoupling.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
        fvm::laplacian(rhorAUf, p_rgh) == fvc::div(phiHbyA)
    );

    p_rghEqn.setReference(pRefCell, pRefValue);

    p_rghEqn.solve();

    if (pvCoupling.finalNonOrthogonalIter())
    {   
        // Correct mass flow rate
        phi = phiHbyA - p_rghEqn.flux();
    }
}

#include "compressibleCorrectVelocity.H"
fvOptions.correct(U);

// Explicitly relax pressure for momentum corrector
p_rgh.relax();
p_rgh.correctBoundaryConditions();

#include "updateAbsolutePressure.H"

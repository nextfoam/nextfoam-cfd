if (MRFEnabled)
{
    MRF.updateForce(rhof, phi, U);
}

volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

if (pvCoupling.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    rhof*(momentumInterpolate(HbyA) & mesh.Sf())
  + fvc::relaxCorr(rhof, phi, U, p, U.select(pvCoupling.finalIter()))
);

// ddtCorr(rho, U, phi) should be handled in each ddt schemes
// currnet implimentation is not valid for Crank-Nicolson 
if (pvCoupling.ddtCorr())
{
    //phiHbyA += MRF.zeroFilter(fvc::ddtCorr(rho, U, phi));
    phiHbyA += fvc::ddtCorr(rho, U, phi);
}

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

fvc::makeRelative(phiHbyA, rho, U);

while (pvCoupling.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
        fvm::ddt(psi, p_rgh)
      - fvm::laplacian(rhorAUf, p_rgh) 
     == 
        fvOptions(psi, p_rgh, rho.name())
      - fvc::div(phiHbyA)
    );

    if (gravityEnabled)
    {
        p_rghEqn += fvc::ddt(hydroStaticDensity);
    }

    if (thermo.incompressible())
    {
        p_rghEqn.setReference(pRefCell, pRefValue);
    }

    p_rghEqn.solve(mesh.solver(p_rgh.select(pvCoupling.finalInnerIter())));

    if (pvCoupling.finalNonOrthogonalIter())
    {
        phi = phiHbyA + p_rghEqn.flux();
    }
}

#include "compressibleCorrectVelocity.H"
fvOptions.correct(U);

// Explicitly relax pressure for momentum corrector
p_rgh.relax();
p_rgh.correctBoundaryConditions();

#include "updateAbsolutePressure.H"

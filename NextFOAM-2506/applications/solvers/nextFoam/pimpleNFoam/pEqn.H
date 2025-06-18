if (MRFEnabled)
{
    MRF.updateForce(phi, U);
}

volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

if (pvCoupling.corrPISO() == pvCoupling.nCorrPISO())
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (momentumInterpolate(HbyA) & mesh.Sf())
  + fvc::relaxCorr(phi, U, p, U.select(pvCoupling.finalIter()))
);

if (pvCoupling.ddtCorr())
{
//   phiHbyA += MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi, Uf));
    //phiHbyA += fvc::ddtCorr(U, phi, Uf);
    //phiHbyA += MRF.zeroFilter(fvc::ddtCorr(U, phi, Uf));
    //phiHbyA += fvc::interpolate(rAU)*fvc::ddtCorr(U, phi, Uf);
    //phiHbyA += fvc::interpolate(rAU)*fvc::ddtCorr(U, Uf());
    //phiHbyA += momentumInterpolate(rAU)*fvc::ddtCorr(U, Uf());
    //phiHbyA += fvc::ddtCorr(U, Uf());
    phiHbyA += fvc::ddtCorr(U, phi);
}

if (p.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p);
    fvc::makeAbsolute(phiHbyA, U);
}


#include "setConsistentFlux.H"

// Add Coriolis flux and make it relative to MRF
if (MRFEnabled)
{
    MRF.addCoriolisFlux(rAUf, phiHbyA);
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAUf, MRF);

// Non-orthogonal pressure corrector loop
while (pvCoupling.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve(mesh.solver(p.select(pvCoupling.finalInnerIter())));

    if (pvCoupling.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

// Momentum corrector
#include "correctVelocity.H"
fvOptions.correct(U);

// Explicitly relax pressure for momentum corrector
p.relax();
p.correctBoundaryConditions();

if (pvCoupling.corrPISO() == pvCoupling.nCorrPISO())
{
    #include "continuityErrs.H"

    fvc::correctUf(Uf, U, phi);

    fvc::makeRelative(phi, U);
}

if (pvCoupling.consistent())
{
    volVectorField velocityCorrection(gradp);
    surfaceScalarField fluxCorrection(fvc::snGrad(p)*mesh.magSf());

    if (MRFEnabled)
    {
        velocityCorrection += MRF.coriolisForce(U);
        fluxCorrection += MRF.volumeAcceleration();
    }

    velocityCorrection *= drAU();
    fluxCorrection *= momentumInterpolate(drAU());

    HbyA = constrainHbyA(HbyA + velocityCorrection, U, p);
    phiHbyA += constrainPhi(fluxCorrection, U, p);
}

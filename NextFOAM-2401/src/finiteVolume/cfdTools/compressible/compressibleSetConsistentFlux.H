if (pvCoupling.consistent())
{
    volVectorField velocityCorrection(gradp_rgh);
    surfaceScalarField fluxCorrection(fvc::snGrad(p_rgh)*mesh.magSf());

    if (gravityEnabled)
    {
        velocityCorrection += gravityForce;
        fluxCorrection += gravityFluxPotential;
    }

    if (MRFEnabled)
    {
        velocityCorrection +=  MRF.coriolisForce(rho, U);
        fluxCorrection += rhof*MRF.volumeAcceleration();
    }

    velocityCorrection *= drAU();
    fluxCorrection *= rhof*momentumInterpolate(drAU());

    HbyA = constrainHbyA(HbyA + velocityCorrection, U, p_rgh);
    phiHbyA += constrainPhi(fluxCorrection, U, p_rgh);
}

gradp_rgh = fvc::gaussGrad(p_rgh);

volVectorField velocityCorrection(gradp_rgh);

if (gravityEnabled)
{
    velocityCorrection += gravityForce;
}

if (MRFEnabled)
{
    velocityCorrection +=  MRF.coriolisForce(rho, U);
}

velocityCorrection *= rAtU();

U = HbyA - velocityCorrection;
U.correctBoundaryConditions();

Umag = postMag(U);


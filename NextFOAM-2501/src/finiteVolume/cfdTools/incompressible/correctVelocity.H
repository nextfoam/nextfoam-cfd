gradp = fvc::gaussGrad(p);

volVectorField velocityCorrection(gradp);

if (MRFEnabled)
{
    velocityCorrection += MRF.coriolisForce(U);
}

velocityCorrection *= rAtU();

U = HbyA - velocityCorrection;
U.correctBoundaryConditions();

Umag = postMag(U);

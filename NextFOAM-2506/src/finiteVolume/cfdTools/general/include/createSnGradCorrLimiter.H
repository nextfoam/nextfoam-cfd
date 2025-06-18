// This file is added by Gill
scalar maxAngle = 80;

scalar kmax = Foam::sqrt
(
    1.0/sqr(Foam::cos(maxAngle*constant::mathematical::pi/180)) - 1.0
);

#if defined(MULTIREGION)
PtrList<surfaceScalarField> snGradCorrLimiterFluid(fluidRegions.size());
PtrList<surfaceScalarField> snGradCorrLimiterSolid(solidRegions.size());

forAll(fluidRegions, i)
{
    snGradCorrLimiterFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "snGradCorrLimiter",
                runTime.timeName(),
                fluidRegions[i]
            ),
            kmax/max(kmax, mag(fluidRegions[i].nonOrthCorrectionVectors()))
        )
    );
}

forAll(solidRegions, i)
{
    snGradCorrLimiterSolid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "snGradCorrLimiter",
                runTime.timeName(),
                solidRegions[i]
            ),
            kmax/max(kmax, mag(solidRegions[i].nonOrthCorrectionVectors()))
        )
    );
}
#else
surfaceScalarField snGradCorrLimiter
(
    IOobject
    (
        "snGradCorrLimiter",
        runTime.timeName(),
        mesh
    ),
    kmax/max(kmax, mag(mesh.nonOrthCorrectionVectors()))
    //mesh,
    //dimensionedScalar("one", dimless, One)
);
#endif

/*
    label limitedFace = 0;
    
    forAll (snGradCorrLimiter.internalField(), i)
    {   
        if (snGradCorrLimiter.internalField()[i] < 1)
        {   
            ++limitedFace;
        }
    }
    
    forAll (snGradCorrLimiter.boundaryField(), patchi)
    {   
        scalarField limiter = snGradCorrLimiter.boundaryField()[patchi];
        
        forAll (limiter, j)
        {   
            if (limiter[j] < 1)
            {
                ++limitedFace;
            }
        }
    }
    
    reduce(limitedFace, sumOp<label>());
*/
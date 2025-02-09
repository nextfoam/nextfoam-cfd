// Limit pressure correction
if (pvCoupling.dict().found("pCorrLimit"))
{
    pCorrLimit = pvCoupling.dict().get<scalar>("pCorrLimit");
    volScalarField pCorr = (p_rgh - p_rgh.prevIter())();
    pCorr.ref() = max(min(pCorr.ref(), pCorrMax), -pCorrMax);
    p_rgh = p_rgh.prevIter() + pCorr;
    p_rgh.correctBoundaryConditions();

    pCorrMax = pCorrLimit*(max(p_rgh.internalField()) + operatingPressure);
}

p = p_rgh + operatingPressure;

if (!thermo.incompressible())
{
    // bounding p
    dimensionedScalar absPressureMin
    (
        "absPressureMin", 
        dimPressure, 
        One, 
        pvCoupling.dict()
    );

    const scalar minP = min(p).value();

    if (minP < absPressureMin.value())
    {
        Info<< "bounding " << p.name()
            << ", min: " << minP
            << " max: " << max(p).value()
            << " average: " << gAverage(p.primitiveField())
            << endl;

        p.primitiveFieldRef() = max
        (
            max
            (
                p.primitiveField(),
                fvc::average(max(p, absPressureMin))().primitiveField()
              * pos0(-p.primitiveField())
            ),
            absPressureMin.value()
        );

        p.boundaryFieldRef() = max(p.boundaryField(), absPressureMin.value());
    }
}


if (gravityEnabled)
{
    p += hydroStaticPressure;
}

if (pvCoupling.dict().lookupOrDefault<Switch>("reportMinMax", false))
{
    scalar Umax = max(mag(U)).value();
    scalar rhoMin = min(rho).value();

    Info << "Umax = " << Umax << endl;
    Info << "rhoMin = " << rhoMin << endl;
}

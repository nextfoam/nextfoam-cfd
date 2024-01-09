// Limit pressure correction
if (pvCoupling.dict().found("pCorrLimit"))
{
    pCorrLimit = pvCoupling.dict().get<scalar>("pCorrLimit");
    volScalarField pCorr = (p_rgh - p_rgh.prevIter())();

    const scalar maxPcorr(max(mag(pCorr)).value());
    if (maxPcorr > pCorrMax.value())
    {
        Info << "Pressure correction has been limited...." << nl << endl;
        pCorr.ref() = max(min(pCorr.ref(), pCorrMax), -pCorrMax);
        p_rgh = p_rgh.prevIter() + pCorr;
        //p_rgh = p_rgh.oldTime() + pCorr;
        p_rgh.correctBoundaryConditions();
    }

    pCorrMax = pCorrLimit*(max(p_rgh.internalField()) + operatingPressure);
}

if (pvCoupling.corrPISO() == pvCoupling.nCorrPISO())
{
    p = p_rgh + operatingPressure;

    if (!thermo.incompressible())
    {
        // bounding p
        const dimensionedScalar absPressureMin
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

            p.boundaryFieldRef() = 
                max(p.boundaryField(), absPressureMin.value());
        }
    }
    
    if (gravityEnabled)
    {   
        p += hydroStaticPressure;
    }
    
    /*if (pressureControl.limit(p))
    {
        p_rgh = p - operatingPressure;

        if (gravityEnabled)
        {   
            p_rgh -= hydroStaticPressure;
        }
        
        p_rgh.correctBoundaryConditions();
    }*/
}

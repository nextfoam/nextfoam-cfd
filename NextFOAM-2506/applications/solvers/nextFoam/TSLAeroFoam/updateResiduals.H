    //- Residual definition from FLUENT
    scalar RMSResRho
    (
        Foam::sqrt
        (
            gSum(sqr((rho - rho.prevIter())().primitiveField()*rDeltaTau))
        )
    );

    vector MSResRhoU
    (
        gSum(cmptMagSqr((rhoU - rhoU.prevIter())().primitiveField()*rDeltaTau))
    );

    vector RMSResRhoU
    (
        Foam::sqrt(MSResRhoU.x()),
        Foam::sqrt(MSResRhoU.y()),
        Foam::sqrt(MSResRhoU.z())
    );

    scalar RMSResRhoE
    (
        Foam::sqrt
        (
            gSum(sqr((rhoE - rhoE.prevIter())().primitiveField()*rDeltaTau))
        )
    );

#ifdef DBNFOAM
    scalar nIter = runTime.value();
#else
    scalar nIter = lusgs.corr();
#endif

    if (nIter <= 5)
    {
        RMSResRhoMax5 = Foam::max(RMSResRho, RMSResRhoMax5);
        RMSResRhoUMax5 = Foam::max(RMSResRhoU, RMSResRhoUMax5);
        RMSResRhoEMax5 = Foam::max(RMSResRhoE, RMSResRhoEMax5);

#ifdef DBNFOAM
        RMSResDict.set("RMSResRhoMax5", RMSResRhoMax5);
        RMSResDict.set("RMSResRhoUMax5", RMSResRhoUMax5);
        RMSResDict.set("RMSResRhoEMax5", RMSResRhoEMax5);
#endif
    }

    RMSResRho = RMSResRho/RMSResRhoMax5;

    Info << "LU-SGS Solver:  Solving for " << "rho"
         << ", Initial residual = " << RMSResRhoOld
         << ", Final residual = " << RMSResRho
         << ", No Iterations " << scalar(1)
         << endl;

    if (dim.x())
    {
        RMSResRhoU.x() = RMSResRhoU.x()/RMSResRhoUMax5.x();

        Info << "LU-SGS Solver:  Solving for " << "rhoUx"
             << ", Initial residual = " << RMSResRhoUOld.x()
             << ", Final residual = " << RMSResRhoU.x()
             << ", No Iterations " << scalar(1)
             << endl;
    }
    else
    {
        RMSResRhoU.x() = 0.0;
    }

    if (dim.y())
    {
        RMSResRhoU.y() = RMSResRhoU.y()/RMSResRhoUMax5.y();

        Info << "LU-SGS Solver:  Solving for " << "rhoUy"
             << ", Initial residual = " << RMSResRhoUOld.y()
             << ", Final residual = " << RMSResRhoU.y()
             << ", No Iterations " << scalar(1)
             << endl;
    }
    else
    {
        RMSResRhoU.y() = 0.0;
    }

    if (dim.z())
    {
        RMSResRhoU.z() = RMSResRhoU.z()/RMSResRhoUMax5.z();

        Info << "LU-SGS Solver:  Solving for " << "rhoUz"
             << ", Initial residual = " << RMSResRhoUOld.z()
             << ", Final residual = " << RMSResRhoU.z()
             << ", No Iterations " << scalar(1)
             << endl;
    }
    else
    {
        RMSResRhoU.z() = 0.0;
    }

    RMSResRhoE = RMSResRhoE/RMSResRhoEMax5;

    Info << "LU-SGS Solver:  Solving for " << "rhoE"
         << ", Initial residual = " << RMSResRhoEOld
         << ", Final residual = " << RMSResRhoE
         << ", No Iterations " << scalar(1)
         << endl;

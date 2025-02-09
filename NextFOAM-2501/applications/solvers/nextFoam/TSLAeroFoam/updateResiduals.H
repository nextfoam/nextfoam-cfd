//Residual definition from FLUENT
scalar RMSResRho
(
    Foam::sqrt
    (
        gSum(sqr((rho - rho.oldTime())().primitiveField()/deltaT))
    )
);

vector MSResRhoU
(
    gSum(cmptMagSqr((rhoU - rhoU.oldTime())().primitiveField()/deltaT))
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
        gSum(sqr((rhoE - rhoE.oldTime())().primitiveField()/deltaT))
    )
);

if (runTime.value() <= 6)
{
    RMSResRhoMax5 = Foam::max(RMSResRho, RMSResRhoMax5);
    RMSResRhoUMax5 = Foam::max(RMSResRhoU, RMSResRhoUMax5);
    RMSResRhoEMax5 = Foam::max(RMSResRhoE, RMSResRhoEMax5);

    RMSResDict.set("RMSResRhoMax5", RMSResRhoMax5);
    RMSResDict.set("RMSResRhoUMax5", RMSResRhoUMax5);
    RMSResDict.set("RMSResRhoEMax5", RMSResRhoEMax5);
}

RMSResRho = RMSResRho/RMSResRhoMax5;
RMSResRhoU.x() = RMSResRhoU.x()/RMSResRhoUMax5.x();
RMSResRhoU.y() = RMSResRhoU.y()/RMSResRhoUMax5.y();
if (mesh.nGeometricD() == 3)
{
    RMSResRhoU.z() = RMSResRhoU.z()/RMSResRhoUMax5.z();
}
RMSResRhoE = RMSResRhoE/RMSResRhoEMax5;


Info << "LU-SGS Solver:  Solving for " << "rho"
        << ", Initial residual = " << RMSResRhoOld 
        << ", Final residual = " << RMSResRho
        << ", No Iterations " << scalar(1)
        << endl;

Info << "LU-SGS Solver:  Solving for " << "rhoUx"
        << ", Initial residual = " << RMSResRhoUOld.x()
        << ", Final residual = " << RMSResRhoU.x()
        << ", No Iterations " << scalar(1)
        << endl;

Info << "LU-SGS Solver:  Solving for " << "rhoUy"
        << ", Initial residual = " << RMSResRhoUOld.y()
        << ", Final residual = " << RMSResRhoU.y()
        << ", No Iterations " << scalar(1)
        << endl;

if (mesh.nGeometricD() == 3)
{
    Info << "LU-SGS Solver:  Solving for " << "rhoUz"
        << ", Initial residual = " << RMSResRhoUOld.z()
        << ", Final residual = " << RMSResRhoU.z()
        << ", No Iterations " << scalar(1)
        << endl;
}

Info << "LU-SGS Solver:  Solving for " << "rhoE"
        << ", Initial residual = " << RMSResRhoEOld
        << ", Final residual = " << RMSResRhoE
        << ", No Iterations " << scalar(1)
        << endl;

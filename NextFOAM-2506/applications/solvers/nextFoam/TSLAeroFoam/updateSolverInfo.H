    //- rho
    List<solverPerformance> rhoSPs;
    solverPerformance spRho
    (
        "LU-SGS",
        "rho",
        RMSResRhoOld,
        RMSResRho,
        1,
        true
    );
    rhoSPs.setSize(rhoSPs.size() + 1, spRho);
    const_cast<dictionary&>
    (
        mesh.data().solverPerformanceDict()
    ).set("rho", rhoSPs);

    //- rhoU
    List<SolverPerformance<vector>> rhoUSPs;
    SolverPerformance<vector> spRhoU
    (
        "LU-SGS",
        "rhoU",
        RMSResRhoUOld,
        RMSResRhoU ,
        vector(1, 1, 1),
        true
    );
    rhoUSPs.setSize(rhoUSPs.size() + 1, spRhoU);
    const_cast<dictionary&>
    (
        mesh.data().solverPerformanceDict()
    ).set("rhoU", rhoUSPs);

    //- rhoE
    List<solverPerformance> rhoESPs;
    solverPerformance spRhoE
    (
        "LU-SGS",
        "rhoE",
        RMSResRhoEOld,
        RMSResRhoE,
        1,
        true
    );
    rhoESPs.setSize(rhoESPs.size() + 1, spRhoE);
    const_cast<dictionary&>
    (
        mesh.data().solverPerformanceDict()
    ).set("rhoE", rhoESPs);

    RMSResRhoOld = RMSResRho;
    RMSResRhoUOld = RMSResRhoU;
    RMSResRhoEOld = RMSResRhoE;

#ifdef DBNFOAM
    RMSResDict.set("RMSResRhoOld", RMSResRho);
    RMSResDict.set("RMSResRhoUOld", RMSResRhoU);
    RMSResDict.set("RMSResRhoEOld", RMSResRhoE);
#endif


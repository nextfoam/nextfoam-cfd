/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION

    PtrList<autoPtr<radiation::radiationModel>> radiationFluid
    (
        fluidRegions.size()
    );
    List<bool> radiationEnabledFluid(fluidRegions.size(), false);

    forAll(fluidRegions, i)
    {
        radiationFluid.set
        (
            i,
            new autoPtr<radiation::radiationModel>
            (
                radiation::radiationModel::New(thermoFluid[i].T())
            )
        );

        #ifdef NO_POST_PROCESS
        radiationEnabledFluid[i] = radiationFluid[i]->radiation();
        #endif
    }

    PtrList<autoPtr<radiation::radiationModel>> radiationSolid
    (
        solidRegions.size()
    );
    List<bool> radiationEnabledSolid(solidRegions.size(), false);

    forAll(solidRegions, i)
    {
        radiationSolid.set
        (
            i,
            new autoPtr<radiation::radiationModel>
            (
                radiation::radiationModel::New(thermoSolid[i].T())
            )
        );

        #ifdef NO_POST_PROCESS
        radiationEnabledSolid[i] = radiationSolid[i]->radiation();
        #endif
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

    autoPtr<radiation::radiationModel> radiation
    (
        radiation::radiationModel::New(thermo.T())
    );

    #ifdef NO_POST_PROCESS
    bool radiationEnabled(radiation->radiation());
    #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

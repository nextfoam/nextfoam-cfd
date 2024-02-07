/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION 

    PtrList<IOMRFZoneList> MRFFluid(fluidRegions.size());
    List<bool> MRFEnabledFluid(fluidRegions.size(), false);

    forAll(fluidRegions, i)
    {
        Info<< "    Adding MRF\n" << endl;
        MRFFluid.set(i, new IOMRFZoneList(fluidRegions[i]));

        #ifdef NO_POST_PROCESS
        MRFEnabledFluid[i] = MRFFluid[i].active();
        #endif
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

    IOMRFZoneList MRF(mesh);

    #ifdef NO_POST_PROCESS
    bool MRFEnabled(MRF.active());
    #endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif 

// ************************************************************************* //

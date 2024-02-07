#include "readTimeControls.H"

/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION

forAll(fluidRegions, i)
{
    correctPhiFluid[i] =
        pvCoupling.dict
        (
            fluidRegions.get(i)
        ).getOrDefault("correctPhi", correctPhiFluid[i]);

    checkMeshCourantNoFluid[i] =
        pvCoupling.dict
        (
            fluidRegions.get(i)
        ).getOrDefault("checkMeshCourantNo", checkMeshCourantNoFluid[i]);

    moveMeshOuterCorrectorsFluid[i] =
        pvCoupling.dict
        (
            fluidRegions.get(i)
        ).getOrDefault
        (
            "moveMeshOuterCorrectors", 
            moveMeshOuterCorrectorsFluid[i]
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

correctPhi = pimple.dict().getOrDefault
(
    "correctPhi",
    correctPhi
);

checkMeshCourantNo = pimple.dict().getOrDefault
(
    "checkMeshCourantNo",
    checkMeshCourantNo
);

moveMeshOuterCorrectors = pimple.dict().getOrDefault
(
    "moveMeshOuterCorrectors",
    moveMeshOuterCorrectors
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

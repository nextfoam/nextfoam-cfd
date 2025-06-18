#include "createControl.H"
#include "createTimeControls.H"

/*---------------------------------------------------------------------------*\
                                Multi Region 
\*---------------------------------------------------------------------------*/
#ifdef MULTIREGION

List<bool> correctPhiFluid(fluidRegions.size());
List<bool> checkMeshCourantNoFluid(fluidRegions.size());
List<bool> moveMeshOuterCorrectorsFluid(fluidRegions.size());

forAll(fluidRegions, i)
{
    correctPhiFluid[i] = 
        pvCoupling.dict
        (
            fluidRegions.get(i)
        ).getOrDefault("correctPhi", fluidRegions[i].dynamic());

    checkMeshCourantNoFluid[i] = 
        pvCoupling.dict
        (
            fluidRegions.get(i)
        ).getOrDefault("checkMeshCourantNo", false);

    moveMeshOuterCorrectorsFluid[i] = 
        pvCoupling.dict
        (
            fluidRegions.get(i)
        ).getOrDefault("moveMeshOuterCorrectors", false);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                                Single Region 
\*---------------------------------------------------------------------------*/
#else

#ifdef NFOAM
pimpleControl& pimple = pvCoupling;
#endif

bool correctPhi
(
    pimple.dict().getOrDefault("correctPhi", mesh.dynamic())
);

bool checkMeshCourantNo
(
    pimple.dict().getOrDefault("checkMeshCourantNo", false)
);

bool moveMeshOuterCorrectors
(
    pimple.dict().getOrDefault("moveMeshOuterCorrectors", false)
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

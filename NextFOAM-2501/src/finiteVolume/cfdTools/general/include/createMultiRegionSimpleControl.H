// This file is added by Gill
multiRegionSimpleControl pvCoupling
(
    runTime,
    reinterpret_cast<PtrList<fvMesh>&>(fluidRegions),
    solidRegions,
    true//foamLog
);

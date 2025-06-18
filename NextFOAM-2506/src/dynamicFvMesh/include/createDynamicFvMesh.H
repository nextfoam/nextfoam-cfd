#if defined(MULTIREGION)

IOobject regionPropertiesIO
(
    "regionProperties",
    runTime.time().constant(),
    runTime.db(),
    IOobjectOption::MUST_READ,
    IOobjectOption::NO_WRITE
);

wordList fluidNames;
wordList solidNames;

if (regionPropertiesIO.typeHeaderOk<IOdictionary>(true))
{
    regionProperties rp(runTime);
    fluidNames = rp["fluid"];
    solidNames = rp["solid"];

    fluidNames.setSize(fluidNames.size());
    solidNames.setSize(solidNames.size());
}
else // Allow multi-region solvers to solve single-region cases.
{
    fluidNames.push_back(Foam::polyMesh::defaultRegion);
    fluidNames.setSize(1);
}

PtrList<dynamicFvMesh> fluidRegions(fluidNames.size());

forAll(fluidNames, i)
{
    Info<< "Create fluid mesh for region " << fluidNames[i]
        << " for time = " << runTime.timeName() << nl << endl;

    fluidRegions.set
    (
        i, 
        dynamicFvMesh::New
        (
            IOobject
            (
                fluidNames[i],
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );
}

PtrList<fvMesh> solidRegions(solidNames.size());

forAll(solidNames, i)
{
    Info<< "Create solid mesh for region " << solidNames[i]
        << " for time = " << runTime.timeName() << nl << endl;

    solidRegions.set
    (
        i, 
        new fvMesh
        (
            IOobject
            (
                solidNames[i],
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );
}


#else
Info<< "Create mesh for time = "
    << runTime.timeName() << nl << endl;

autoPtr<dynamicFvMesh> meshPtr(dynamicFvMesh::New(args, runTime));

dynamicFvMesh& mesh = meshPtr();
#endif

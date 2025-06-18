IOdictionary cellZonesToRegions
(
    IOobject
    (
        "cellZonesToRegions",
        runTime.constant(),
        runTime,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

createDefaultRegions(cellZonesToRegions, cellZoneName);

cellZonesToRegions.regIOobject::write();

// Edit regions section in cellZonesToRegions file for splitMeshRegions
if (port == 0)
{
    fileName str;

    Info << "The default region configuration is as follows:" 
            << nl << endl;

    Info << "regions" 
            << cellZonesToRegions.subDict("regions") << nl << endl;

    Info << "After correctly modifying " 
            << "the dictionary 'regions' of the file "
            << cellZonesToRegions.filePath() 
            << ", press enter."
            << endl;
    getline(cin, str);
}
else
{
    if (port > 65535)
    {
        Info << "Port number " << port << " is invalid..." << nl << endl;
        return -1;
    };

    waitMapping((uint16_t)port);
}

cellZonesToRegions.regIOobject::read();

const dictionary& regionsDict(cellZonesToRegions.subDict("regions"));

cellZoneType.clear();

forAll(regionsDict.toc(), i)
{
    word regionName = regionsDict.toc()[i];

    const dictionary& regionDict(regionsDict.subDict(regionName));

    word regionType = regionDict.get<word>("type");

    List<word> czs = regionDict.get<List<word>>("cellZones");
    List<label> czid(czs.size());

    forAll(czs, j)
    {
        label cellZoneID = pMesh.cellZones().findZoneID(czs[j]);
        czid[j] = cellZoneID;

        label zoneID(cellGroupZoneID[cellZoneID]);

        groupType.erase(zoneID);
        groupType.insert(zoneID, regionType);

        cellZoneType(cellZoneID).append(regionType);
        cellZoneToRegionName(cellZoneID).append(regionName);
    }

    if (czid.size())
    {
        if (regionType == "fluid") 
        {
            fluidRegions.append(czid);
            fluidRegionNames.append(regionName);
        }
        else
        {
            solidRegions.append(czid);
            solidRegionNames.append(regionName);
        }

        regionToZones.append(czid);
    }
}

fluidRegions.shrink();
fluidRegionNames.shrink();
solidRegions.shrink();
solidRegionNames.shrink();
regionToZones.shrink();

if (fluidRegions.size() == 1 && !solidRegions.size())
{
    word name = polyMesh::defaultRegion;
    fluidRegionNames[0] = name;
    cellZoneToRegionName.clear();

    forAll(pMesh.cellZones(), cellZoneI)
    {
        cellZoneToRegionName(cellZoneI).append(name);
    }
}

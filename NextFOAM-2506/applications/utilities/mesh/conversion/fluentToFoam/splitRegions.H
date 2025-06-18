void createDefaultRegions
(
    IOdictionary& cellZonesToRegions,
    const Map<word>& cellZoneName
)
{
    if (cellZonesToRegions.found("fluentCellZones"))
    {
        cellZonesToRegions.subDict("fluentCellZones").clear();
    }

    dictionary cellZonesDict;
    forAll(cellZoneName, i)
    {
        cellZonesDict.add(cellZoneName[i], cellZoneType[i]);
    }
    cellZonesToRegions.add("fluentCellZones", cellZonesDict, true);

    word defaultRegionName("fluid");
    DynamicList<word> fluidCellZones;
    DynamicList<word> solidCellZones;

    forAll(cellZonesDict.toc(), zoneI)
    {
        word name(cellZonesDict.toc()[zoneI]);
        word type(cellZonesDict.get<word>(name));

        if (type == "fluid")
        {
            fluidCellZones.append(name);
        }
        else
        {
            solidCellZones.append(name);
        }
    }

    fluidCellZones.shrink();
    solidCellZones.shrink();

    if (!solidCellZones.size())
    {
        defaultRegionName = "defaultRegion";
    }

    dictionary regions;

    if (fluidCellZones.size())
    {
        dictionary region;
        region.add("cellZones", flatOutput(fluidCellZones), true);
        region.add("type", "fluid", true);

        regions.add(defaultRegionName, region, true);
    }

    if (solidCellZones.size())
    {
        forAll(solidCellZones, i)
        {
            dictionary region;
            word name(solidCellZones[i]);
            wordHashSet zones;
            zones.insert(name);

            region.add("cellZones", flatOutput(zones), true);
            region.add("type", "solid", true);

            regions.add(name, region, true);
        }
    }

    if (cellZonesToRegions.found("regions"))
    {
        cellZonesToRegions.subDict("regions").clear();
    }

    cellZonesToRegions.add("regions", regions, true);
}


// Prepend prefix to selected patches.
void renamePatches
(
    fvMesh& mesh,
    const word& prefix,
    const labelList& patchesToRename
)
{
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    forAll(patchesToRename, i)
    {
        label patchi = patchesToRename[i];
        polyPatch& pp = polyPatches[patchi];

        if (isA<coupledPolyPatch>(pp))
        {
            WarningInFunction
                << "Encountered coupled patch " << pp.name()
                << ". Will only rename the patch itself,"
                << " not any referred patches."
                << " This might have to be done by hand."
                << endl;
        }

        pp.name() = prefix + '_' + pp.name();
    }
    // Recalculate any demand driven data (e.g. group to name lookup)
    polyPatches.updateMesh();
}   


// Select all cells not in the region
labelList getNonRegionCells(const labelList& cellRegion, const label regionI)
{   
    DynamicList<label> nonRegionCells(cellRegion.size());
    forAll(cellRegion, celli)
    {   
        if (cellRegion[celli] != regionI)
        {
            nonRegionCells.append(celli);
        } 
    }
    return nonRegionCells.shrink();
}   


void addToInterface
(
    const polyMesh& mesh,
    const label zoneID,
    const label ownRegion,
    const label neiRegion,
    EdgeMap<Map<label>>& regionsToSize
)
{
    edge interface
    (
        min(ownRegion, neiRegion),
        max(ownRegion, neiRegion)
    );

    auto iter = regionsToSize.find(interface);

    if (iter.good())
    {
        // Check if zone present
        auto zoneIter = iter().find(zoneID);
        if (zoneIter.good())
        {
            // Found zone. Increment count.
            ++(*zoneIter);
        }
        else
        {
            // New or no zone. Insert with count 1.
            iter().insert(zoneID, 1);
        }
    }
    else
    {
        // Create new interface of size 1.
        Map<label> zoneToSize;
        zoneToSize.insert(zoneID, 1);
        regionsToSize.insert(interface, zoneToSize);
    }
}


// Get region-region interface name and sizes.
// Returns interfaces as straight list for looping in identical order.
void getInterfaceSizes
(
    const polyMesh& mesh,
    const bool useFaceZones,
    const labelList& cellRegion,
    const wordList& regionNames,

    edgeList& interfaces,
    List<Pair<word>>& interfaceNames,
    labelList& interfaceSizes,
    labelList& faceToInterface
)
{
    // From region-region to faceZone (or -1) to number of faces.

    EdgeMap<Map<label>> regionsToSize;


    // Internal faces
    // ~~~~~~~~~~~~~~

    forAll(mesh.faceNeighbour(), facei)
    {
        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = cellRegion[mesh.faceNeighbour()[facei]];

        if (ownRegion != neiRegion)
        {
            addToInterface
            (
                mesh,
                (useFaceZones ? mesh.faceZones().whichZone(facei) : -1),
                ownRegion,
                neiRegion,
                regionsToSize
            );
        }
    }


    // Boundary faces
    // ~~~~~~~~~~~~~~

    // Neighbour cellRegion.
    labelList coupledRegion(mesh.nBoundaryFaces());

    forAll(coupledRegion, i)
    {
        label celli = mesh.faceOwner()[i+mesh.nInternalFaces()];
        coupledRegion[i] = cellRegion[celli];
    }
    syncTools::swapBoundaryFaceList(mesh, coupledRegion);

    forAll(coupledRegion, i)
    {
        label facei = i+mesh.nInternalFaces();
        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = coupledRegion[i];

        if (ownRegion != neiRegion)
        {
            addToInterface
            (
                mesh,
                (useFaceZones ? mesh.faceZones().whichZone(facei) : -1),
                ownRegion,
                neiRegion,
                regionsToSize
            );
        }
    }


    if (Pstream::parRun())
    {
        if (Pstream::master())
        {
            // Receive and add to my sizes
            for (const int slave : Pstream::subProcs())
            {
                IPstream fromSlave(Pstream::commsTypes::blocking, slave);

                EdgeMap<Map<label>> slaveSizes(fromSlave);

                forAllConstIters(slaveSizes, slaveIter)
                {
                    const Map<label>& slaveInfo = *slaveIter;

                    auto masterIter = regionsToSize.find(slaveIter.key());

                    if (masterIter.good())
                    {
                        // Same inter-region
                        Map<label>& masterInfo = *masterIter;

                        forAllConstIters(slaveInfo, iter)
                        {
                            const label zoneID = iter.key();
                            const label slaveSize = iter.val();

                            auto zoneIter = masterInfo.find(zoneID);
                            if (zoneIter.good())
                            {
                                *zoneIter += slaveSize;
                            }
                            else
                            {
                                masterInfo.insert(zoneID, slaveSize);
                            }
                        }
                    }
                    else
                    {
                        regionsToSize.insert(slaveIter.key(), slaveInfo);
                    }
                }
            }
        }
        else
        {
            // Send to master
            {
                OPstream toMaster
                (
                    Pstream::commsTypes::blocking,
                    Pstream::masterNo()
                );
                toMaster << regionsToSize;
            }
        }
    }

    // Rework

    Pstream::broadcast(regionsToSize);


    // Now we have the global sizes of all inter-regions.
    // Invert this on master and distribute.
    label nInterfaces = 0;
    forAllConstIters(regionsToSize, iter)
    {
        const Map<label>& info = iter.val();
        nInterfaces += info.size();
    }

    interfaces.setSize(nInterfaces);
    interfaceNames.setSize(nInterfaces);
    interfaceSizes.setSize(nInterfaces);
    EdgeMap<Map<label>> regionsToInterface(nInterfaces);

    nInterfaces = 0;
    forAllConstIters(regionsToSize, iter)
    {
        const edge& e = iter.key();
        const Map<label>& info = iter.val();

        const word& name0 = regionNames[e[0]];
        const word& name1 = regionNames[e[1]];

        forAllConstIters(info, infoIter)
        {
            interfaces[nInterfaces] = iter.key();
            label zoneID = infoIter.key();
            if (zoneID == -1)
            {
                interfaceNames[nInterfaces] = Pair<word>
                (
                    name0 + "_to_" + name1,
                    name1 + "_to_" + name0
                );
            }
            else
            {
                const word& zoneName = mesh.faceZones()[zoneID].name();
                interfaceNames[nInterfaces] = Pair<word>
                (
                    zoneName + "_" + name0 + "_to_" + name1,
                    zoneName + "_" + name1 + "_to_" + name0
                );
            }
            interfaceSizes[nInterfaces] = infoIter();

            if (regionsToInterface.found(e))
            {
                regionsToInterface[e].insert(zoneID, nInterfaces);
            }
            else
            {
                Map<label> zoneAndInterface;
                zoneAndInterface.insert(zoneID, nInterfaces);
                regionsToInterface.insert(e, zoneAndInterface);
            }
            nInterfaces++;
        }
    }


    // Consistent interface information for all processors
    Pstream::broadcasts
    (
        UPstream::worldComm,
        interfaces,
        interfaceNames,
        interfaceSizes,
        regionsToInterface
    );

    // Mark all inter-region faces.
    faceToInterface.setSize(mesh.nFaces(), -1);

    forAll(mesh.faceNeighbour(), facei)
    {
        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = cellRegion[mesh.faceNeighbour()[facei]];

        if (ownRegion != neiRegion)
        {
            label zoneID = -1;
            if (useFaceZones)
            {
                zoneID = mesh.faceZones().whichZone(facei);
            }

            edge interface
            (
                min(ownRegion, neiRegion),
                max(ownRegion, neiRegion)
            );

            faceToInterface[facei] = regionsToInterface[interface][zoneID];
        }
    }
    forAll(coupledRegion, i)
    {
        label facei = i+mesh.nInternalFaces();
        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = coupledRegion[i];

        if (ownRegion != neiRegion)
        {
            label zoneID = -1;
            if (useFaceZones)
            {
                zoneID = mesh.faceZones().whichZone(facei);
            }

            edge interface
            (
                min(ownRegion, neiRegion),
                max(ownRegion, neiRegion)
            );

            faceToInterface[facei] = regionsToInterface[interface][zoneID];
        }
    }
}


// Create mesh for region.
autoPtr<mapPolyMesh> createRegionMesh
(
    const fvMesh& mesh,
    // Region info
    const labelList& cellRegion,
    const label regionI,
    const word& regionName,
    // Interface info
    const labelList& interfacePatches,
    const labelList& faceToInterface,

    autoPtr<fvMesh>& newMesh
)
{
    // Create dummy system/fv*
    fvMeshTools::createDummyFvMeshFiles(mesh, regionName, true);

    // Neighbour cellRegion.
    labelList coupledRegion(mesh.nBoundaryFaces());

    forAll(coupledRegion, i)
    {
        label celli = mesh.faceOwner()[i+mesh.nInternalFaces()];
        coupledRegion[i] = cellRegion[celli];
    }
    syncTools::swapBoundaryFaceList(mesh, coupledRegion);


    // Topology change container. Start off from existing mesh.
    polyTopoChange meshMod(mesh);

    // Cell remover engine
    removeCells cellRemover(mesh);

    // Select all but region cells
    labelList cellsToRemove(getNonRegionCells(cellRegion, regionI));

    // Find out which faces will get exposed. Note that this
    // gets faces in mesh face order. So both regions will get same
    // face in same order (important!)
    labelList exposedFaces = cellRemover.getExposedFaces(cellsToRemove);

    labelList exposedPatchIDs(exposedFaces.size());
    forAll(exposedFaces, i)
    {
        label facei = exposedFaces[i];
        label interfacei = faceToInterface[facei];

        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = -1;

        if (mesh.isInternalFace(facei))
        {
            neiRegion = cellRegion[mesh.faceNeighbour()[facei]];
        }
        else
        {
            neiRegion = coupledRegion[facei-mesh.nInternalFaces()];
        }


        // Check which side is being kept - determines which of the two
        // patches will be used.

        label otherRegion = -1;

        if (ownRegion == regionI && neiRegion != regionI)
        {
            otherRegion = neiRegion;
        }
        else if (ownRegion != regionI && neiRegion == regionI)
        {
            otherRegion = ownRegion;
        }
        else
        {
            FatalErrorInFunction
                << "Exposed face:" << facei
                << " fc:" << mesh.faceCentres()[facei]
                << " has owner region " << ownRegion
                << " and neighbour region " << neiRegion
                << " when handling region:" << regionI
                << exit(FatalError);
        }

        // Find the patch.
        if (regionI < otherRegion)
        {
            exposedPatchIDs[i] = interfacePatches[interfacei];
        }
        else
        {
            exposedPatchIDs[i] = interfacePatches[interfacei]+1;
        }
    }

    // Remove faces
    cellRemover.setRefinement
    (
        cellsToRemove,
        exposedFaces,
        exposedPatchIDs,
        meshMod
    );

    autoPtr<mapPolyMesh> map = meshMod.makeMesh
    (
        newMesh,
        IOobject
        (
            regionName,
            mesh.time().timeName(),
            mesh.time(),
            IOobject::READ_IF_PRESENT,  // read fv* if present
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    return map;
}


void createAndWriteRegion
(
    const fvMesh& mesh,
    const labelList& cellRegion,
    const wordList& regionNames,
    const bool prefixRegion,
    const labelList& faceToInterface,
    const labelList& interfacePatches,
    const label regionI,
    const word& newMeshInstance
)
{
    Info<< "Creating mesh for region " << regionI
        << ' ' << regionNames[regionI] << endl;

    autoPtr<fvMesh> newMesh;
    autoPtr<mapPolyMesh> map = createRegionMesh
    (
        mesh,
        cellRegion,
        regionI,
        regionNames[regionI],
        interfacePatches,
        faceToInterface,
        newMesh
    );


    // Make map of all added patches
    labelHashSet addedPatches(2*interfacePatches.size());
    forAll(interfacePatches, interfacei)
    {
        addedPatches.insert(interfacePatches[interfacei]);
        addedPatches.insert(interfacePatches[interfacei]+1);
    }


    Info<< "Mapping fields" << endl;

    // Map existing fields
    newMesh().updateMesh(map());

    const polyBoundaryMesh& newPatches = newMesh().boundaryMesh();
    newPatches.checkParallelSync(true);

    // Delete empty patches
    // ~~~~~~~~~~~~~~~~~~~~

    // Create reordering list to move patches-to-be-deleted to end
    labelList oldToNew(newPatches.size(), -1);
    DynamicList<label> sharedPatches(newPatches.size());
    label newI = 0;

    Info<< "Deleting empty patches" << endl;

    // Assumes all non-proc boundaries are on all processors!
    forAll(newPatches, patchi)
    {
        const polyPatch& pp = newPatches[patchi];

        if (!isA<processorPolyPatch>(pp))
        {
            if (returnReduceOr(pp.size()))
            {
                oldToNew[patchi] = newI;
                if (!addedPatches.found(patchi))
                {
                    sharedPatches.append(newI);
                }
                newI++;
            }
        }
    }

    // Same for processor patches (but need no reduction)
    forAll(newPatches, patchi)
    {
        const polyPatch& pp = newPatches[patchi];

        if (isA<processorPolyPatch>(pp) && pp.size())
        {
            oldToNew[patchi] = newI++;
        }
    }

    const label nNewPatches = newI;

    // Move all deleteable patches to the end
    forAll(oldToNew, patchi)
    {
        if (oldToNew[patchi] == -1)
        {
            oldToNew[patchi] = newI++;
        }
    }

    //reorderPatches(newMesh(), oldToNew, nNewPatches);
    fvMeshTools::reorderPatches(newMesh(), oldToNew, nNewPatches, true);

    // Rename shared patches with region name
    if (prefixRegion)
    {
        Info<< "Prefixing patches with region name" << endl;

        renamePatches(newMesh(), regionNames[regionI], sharedPatches);
    }


    Info<< "Writing new mesh" << endl;

    newMesh().setInstance(newMeshInstance);
    newMesh().write();
    topoSet::removeFiles(newMesh());
    processorMeshes::removeFiles(newMesh());

    // Write addressing files like decomposePar
    Info<< "Writing addressing to base mesh" << endl;

    labelIOList pointProcAddressing
    (
        IOobject
        (
            "pointRegionAddressing",
            newMesh().facesInstance(),
            polyMesh::meshSubDir,
            newMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        map().pointMap()
    );
    Info<< "Writing map " << pointProcAddressing.name()
        << " from region" << regionI
        << " points back to base mesh." << endl;
    pointProcAddressing.write();

    labelIOList faceProcAddressing
    (
        IOobject
        (
            "faceRegionAddressing",
            newMesh().facesInstance(),
            polyMesh::meshSubDir,
            newMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        newMesh().nFaces()
    );
    forAll(faceProcAddressing, facei)
    {
        // face + turning index. (see decomposePar)
        // Is the face pointing in the same direction?
        label oldFacei = map().faceMap()[facei];

        if
        (
            map().cellMap()[newMesh().faceOwner()[facei]]
         == mesh.faceOwner()[oldFacei]
        )
        {
            faceProcAddressing[facei] = oldFacei+1;
        }
        else
        {
            faceProcAddressing[facei] = -(oldFacei+1);
        }
    }
    Info<< "Writing map " << faceProcAddressing.name()
        << " from region" << regionI
        << " faces back to base mesh." << endl;
    faceProcAddressing.write();

    labelIOList cellProcAddressing
    (
        IOobject
        (
            "cellRegionAddressing",
            newMesh().facesInstance(),
            polyMesh::meshSubDir,
            newMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        map().cellMap()
    );
    Info<< "Writing map " <<cellProcAddressing.name()
        << " from region" << regionI
        << " cells back to base mesh." << endl;
    cellProcAddressing.write();

    labelIOList boundaryProcAddressing
    (
        IOobject
        (
            "boundaryRegionAddressing",
            newMesh().facesInstance(),
            polyMesh::meshSubDir,
            newMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        labelList(nNewPatches, -1)
    );
    forAll(oldToNew, i)
    {
        if (!addedPatches.found(i))
        {
            label newI = oldToNew[i];
            if (newI >= 0 && newI < nNewPatches)
            {
                boundaryProcAddressing[oldToNew[i]] = i;
            }
        }
    }
    Info<< "Writing map " << boundaryProcAddressing.name()
        << " from region" << regionI
        << " boundary back to base mesh." << endl;
    boundaryProcAddressing.write();
}


// Create for every region-region interface with non-zero size two patches.
// First one is for minimumregion to maximumregion.
// Note that patches get created in same order on all processors (if parallel)
// since looping over synchronised 'interfaces'.
labelList addRegionPatches
(
    fvMesh& mesh,
    const wordList& regionNames,
    const edgeList& interfaces,
    const List<Pair<word>>& interfaceNames
)
{
    Info<< nl << "Adding patches" << nl << endl;

    labelList interfacePatches(interfaces.size());

    forAll(interfaces, interI)
    {
        const edge& e = interfaces[interI];
        const Pair<word>& names = interfaceNames[interI];

        //Info<< "For interface " << interI
        //    << " between regions " << e
        //    << " trying to add patches " << names << endl;


        mappedWallPolyPatch patch1
        (
            names[0],
            0,                  // overridden
            0,                  // overridden
            0,                  // overridden
            regionNames[e[1]],  // sampleRegion
            mappedPatchBase::NEARESTPATCHFACE,
            names[1],           // samplePatch
            point::zero,        // offset
            mesh.boundaryMesh()
        );

        interfacePatches[interI] = fvMeshTools::addPatch
        (
            mesh,
            patch1,
            dictionary(),   //optional per field value
            fvPatchFieldBase::calculatedType(),
            true            //validBoundary
        );

        mappedWallPolyPatch patch2
        (
            names[1],
            0,
            0,
            0,
            regionNames[e[0]],  // sampleRegion
            mappedPatchBase::NEARESTPATCHFACE,
            names[0],
            point::zero,        // offset
            mesh.boundaryMesh()
        );
        fvMeshTools::addPatch
        (
            mesh,
            patch2,
            dictionary(),   //optional per field value
            fvPatchFieldBase::calculatedType(),
            true            //validBoundary
        );

        Info<< "For interface between region " << regionNames[e[0]]
            << " and " << regionNames[e[1]] << " added patches" << endl
            << "    " << interfacePatches[interI]
            << "\t" << mesh.boundaryMesh()[interfacePatches[interI]].name()
            << endl
            << "    " << interfacePatches[interI]+1
            << "\t" << mesh.boundaryMesh()[interfacePatches[interI]+1].name()
            << endl;
    }
    return interfacePatches;
}


void getCellRegion
(   
    const polyMesh& mesh,
    const cellZoneMesh& regions,
    labelList& cellRegion
)   
{   
    forAll(regions, regionI)
    {
        const cellZone& cz = regions[regionI];

        forAll(cz, i)
        {
            label celli = cz[i];
            if (cellRegion[celli] == -1)
            {
                cellRegion[celli] = regionI;
            }
            else
            {
                FatalErrorInFunction
                    << "Cell " << celli << " with cell centre "
                    << mesh.cellCentres()[celli]
                    << " is multiple zones. This is not allowed." << endl
                    << "It is in zone " << regions.names()[cellRegion[celli]]
                    << " and in zone " << regions.names()[regionI]
                    << exit(FatalError);
            }
        }
    }

    // Neighbour zoneID.
    //syncTools::swapBoundaryCellList(mesh, cellRegion, neiClusterID);
}


void writeCellToRegion(const fvMesh& mesh, const labelList& cellRegion)
{
    // Write to manual decomposition option
    {
        labelIOList cellToRegion
        (
            IOobject
            (
                "cellToRegion",
                mesh.facesInstance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            cellRegion
        );
        cellToRegion.write();

        Info<< "Writing region per cell file (for manual decomposition) to "
            << cellToRegion.objectPath() << nl << endl;
    }
    // Write for postprocessing
    {
        volScalarField cellToRegion
        (
            IOobject
            (
                "cellToRegion",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(dimless, Zero),
            fvPatchFieldBase::zeroGradientType()
        );
        forAll(cellRegion, celli)
        {
            cellToRegion[celli] = cellRegion[celli];
        }
        cellToRegion.write();

        Info<< "Writing region per cell as volScalarField to "
            << cellToRegion.objectPath() << nl << endl;
    }
}


bool splitMeshRegions
(
    const Time& runTime,
    fvMesh& mesh,
    const cellZoneMesh& regions,
    const List<labelList>& regionToZones,
    const List<word>& fluidRegionNames,
    const List<word>& solidRegionNames
)
{
    const word oldInstance = mesh.pointsInstance();

    mesh.cellZones().checkParallelSync(true);

    // Determine per cell the region it belongs to
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // cellRegion is the labelList with the region per cell.
    labelList cellRegion(mesh.nCells(), -1);

    // Name of region
    wordList regionNames(regions.names());

    label nCellRegions(regions.size());

    getCellRegion(mesh, regions, cellRegion);

    Info<< endl << "Number of regions:" << nCellRegions << nl << endl;

    // Write decomposition to file
    //writeCellToRegion(mesh, cellRegion);


    // Sizes per region
    // ~~~~~~~~~~~~~~~~

    labelList regionSizes(nCellRegions, Zero);

    forAll(cellRegion, celli)
    {
        regionSizes[cellRegion[celli]]++;
    }
    forAll(regionSizes, regionI)
    {
        reduce(regionSizes[regionI], sumOp<label>());
    }

    Info<< "Region\tCells" << nl
        << "------\t-----" << endl;

    forAll(regionSizes, regionI)
    {
        Info<< regionI << '\t' << regionSizes[regionI] << nl;
    }
    Info<< endl;



    // Print region to zone
    Info<< "Region\tZone\tName" << nl
        << "------\t----\t----" << endl;
    forAll(regionToZones, regionI)
    {
        Info<< regionI << '\t' << flatOutput(regionToZones[regionI])
            << '\t'
            << regionNames[regionI] << nl;
    }
    Info<< endl;


    // Since we're going to mess with patches and zones make sure all
    // is synchronised
    mesh.boundaryMesh().checkParallelSync(true);
    mesh.faceZones().checkParallelSync(true);


    // Interfaces
    // ----------
    // per interface:
    // - the two regions on either side
    // - the name
    // - the (global) size
    edgeList interfaces;
    List<Pair<word>> interfaceNames;
    labelList interfaceSizes;
    // per face the interface
    labelList faceToInterface;

    getInterfaceSizes
    (
        mesh,
        false,//true,
        cellRegion,
        regionNames,
        interfaces,
        interfaceNames,
        interfaceSizes,
        faceToInterface
    );

    Info<< "Sizes of interfaces between regions:" << nl << nl
        << "Interface\tRegion\tRegion\tFaces" << nl
        << "---------\t------\t------\t-----" << endl;

    forAll(interfaces, interI)
    {
        const edge& e = interfaces[interI];

        Info<< interI
            << "\t\t" << e[0] << "\t" << e[1]
            << "\t" << interfaceSizes[interI] << nl;
    }
    Info<< endl;

    // Remove any demand-driven fields ('S', 'V' etc)
    mesh.clearOut();

    {
        // Add patches for interfaces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Add all possible patches. Empty ones get filtered later on.
        Info<< nl << "Adding patches" << nl << endl;

        labelList interfacePatches
        (
            addRegionPatches
            (
                mesh,
                regionNames,
                interfaces,
                interfaceNames
            )
        );

        // Create regions
        // ~~~~~~~~~~~~~~
        {
            // Split all
            for (label regionI = 0; regionI < nCellRegions; regionI++)
            {
                Info<< nl
                    << "Region " << regionI << nl
                    << "-------- " << endl;

                createAndWriteRegion
                (
                    mesh,
                    cellRegion,
                    regionNames,
                    false,//prefixRegion,
                    faceToInterface,
                    interfacePatches,
                    regionI,
                    oldInstance
                    //(overwrite ? oldInstance : runTime.timeName())
                );
            }
        }
    }

    // ... Write constant/regionProperties
    fileName regionProperties
    (
        runTime.path()/runTime.constant()/"regionProperties"
    );

    {
        Foam::mkDir(runTime.path()/runTime.constant());

        OFstream osrp(regionProperties);

        IOobject::writeBanner(osrp)
            << "FoamFile\n{\n"
            << "    version     " << osrp.version() << ";\n"
            << "    format      " << osrp.format() << ";\n"
            << "    class       " << "dictionary" << ";\n"
            << "    location    " << "\"constant\"" << ";\n"
            << "    object      " << "regionProperties" << ";\n"
            << "}" << nl;

        IOobject::writeDivider(osrp)
            << nl
            << "regions"
            << nl << "("
            << nl << "    fluid ( ";

        forAll(fluidRegionNames, iter)
        {
            //osrp << iter.key() << " ";
            osrp << fluidRegionNames[iter] << " ";
        }

        osrp  << ")"
              << nl << "    solid ( ";

        forAll(solidRegionNames, iter)
        {
            //osrp << iter.key() << " ";
            osrp << solidRegionNames[iter] << " ";
        }
        osrp  << ")"
              << nl << ");" << nl;

        IOobject::writeEndDivider(osrp);
    }

    return true;
}

label addBafflePatch
(
    fvMesh& mesh,
    const word& patchName,
    const word& coupledPatchName,
    const word& groupName,
    const word& extraGroupName,
    const dictionary& patchDict,
    const word& regionName,
    const bool addExtraGroup
)   
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    if (pbm.findPatchID(patchName) == -1)
    {
        autoPtr<mappedWallPolyPatch> ppPtr
        (
	        new mappedWallPolyPatch
            (
                patchName,
                0,
                0,
                0,
                regionName,
                mappedPatchBase::NEARESTPATCHFACE,
                coupledPatchName,
                point::zero,
                pbm
            )
        );
        auto& pp = *ppPtr;

        pp.inGroups().clear();
        pp.inGroups().append(groupName);

        if (addExtraGroup)
        {   
            pp.inGroups().append(extraGroupName);
        }

        // Add patch, create calculated everywhere
        fvMeshTools::addPatch
        (
            mesh,
            pp,
            dictionary(),   // do not set specialised patchFields
            fvPatchFieldBase::calculatedType(),
            true            // parallel sync'ed addition
        );
    }
    else
    {
        Info<< "Patch '" << patchName
            << "' already exists.  Only "
            << "moving patch faces - type will remain the same"
            << endl;
    }

    return pbm.findPatchID(patchName);
}


// Filter out the empty patches.
void filterPatches(fvMesh& mesh, const wordHashSet& addedPatchNames)
{
    // Remove any zero-sized ones. Assumes
    // - processor patches are already only there if needed
    // - all other patches are available on all processors
    // - but coupled ones might still be needed, even if zero-size
    //   (e.g. processorCyclic)
    // See also logic in createPatch.
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    labelList oldToNew(pbm.size(), -1);
    label newPatchi = 0;
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (!isA<processorPolyPatch>(pp))
        {
            if
            (
                isA<coupledPolyPatch>(pp)
             || returnReduceOr(pp.size())
             || addedPatchNames.contains(pp.name())
            )
            {
                // Coupled (and unknown size) or uncoupled and used
                oldToNew[patchi] = newPatchi++;
            }
        }
    }

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            oldToNew[patchi] = newPatchi++;
        }
    }


    const label nKeepPatches = newPatchi;

    // Shuffle unused ones to end
    if (nKeepPatches != pbm.size())
    {
        Info<< endl
            << "Removing zero-sized patches:" << endl << incrIndent;

        forAll(oldToNew, patchi)
        {
            if (oldToNew[patchi] == -1)
            {
                Info<< indent << pbm[patchi].name()
                    << " type " << pbm[patchi].type()
                    << " at position " << patchi << endl;
                oldToNew[patchi] = newPatchi++;
            }
        }
        Info<< decrIndent;

        fvMeshTools::reorderPatches(mesh, oldToNew, nKeepPatches, true);
        Info<< endl;
    }
}


void modifyOrAddFace
(
    polyTopoChange& meshMod,
    const face& f,
    const label facei,
    const label own,
    const bool flipFaceFlux,
    const label newPatchi,
    const label zoneID,
    const bool zoneFlip,

    bitSet& modifiedFace
)
{
    if (modifiedFace.set(facei))
    {
        // First usage of face. Modify.
        meshMod.setAction
        (
            polyModifyFace
            (
                f,                          // modified face
                facei,                      // label of face
                own,                        // owner
                -1,                         // neighbour
                flipFaceFlux,               // face flip
                newPatchi,                  // patch for face
                false,                      // remove from zone
                zoneID,                     // zone for face
                zoneFlip                    // face flip in zone
            )
        );
    }
    else
    {
        // Second or more usage of face. Add.
        meshMod.setAction
        (
            polyAddFace
            (
                f,                          // modified face
                own,                        // owner
                -1,                         // neighbour
                -1,                         // master point
                -1,                         // master edge
                facei,                      // master face
                flipFaceFlux,               // face flip
                newPatchi,                  // patch for face
                zoneID,                     // zone for face
                zoneFlip                    // face flip in zone
            )
        );
    }
}


// Create faces for fZone faces. Usually newMasterPatches, newSlavePatches
// only size one but can be more for duplicate baffle sets
void createFaces
(
    const bool internalFacesOnly,
    const fvMesh& mesh,
    const faceZone& fZone,
    const labelUList& newMasterPatches,
    const labelUList& newSlavePatches,
    polyTopoChange& meshMod,
    bitSet& modifiedFace,
    label& nModified
)
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    forAll(newMasterPatches, i)
    {
        const label newMasterPatchi = newMasterPatches[i];
        const label newSlavePatchi = newSlavePatches[i];

        // Pass 1. Do selected side of zone
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
        {
            const label zoneFacei = fZone.whichFace(facei);

            if (zoneFacei != -1)
            {
                if (!fZone.flipMap()[zoneFacei])
                {
                    // Use owner side of face
                    modifyOrAddFace
                    (
                        meshMod,
                        mesh.faces()[facei],    // modified face
                        facei,                  // label of face
                        mesh.faceOwner()[facei],// owner
                        false,                  // face flip
                        newMasterPatchi,        // patch for face
                        fZone.index(),          // zone for face
                        false,                  // face flip in zone
                        modifiedFace            // modify or add status
                    );
                }
                else
                {
                    // Use neighbour side of face.
                    // To keep faceZone pointing out of original neighbour
                    // we don't need to set faceFlip since that cell
                    // now becomes the owner
                    modifyOrAddFace
                    (
                        meshMod,
                        mesh.faces()[facei].reverseFace(),  // modified face
                        facei,                      // label of face
                        mesh.faceNeighbour()[facei],// owner
                        true,                       // face flip
                        newMasterPatchi,            // patch for face
                        fZone.index(),              // zone for face
                        false,                      // face flip in zone
                        modifiedFace                // modify or add status
                    );
                }

                ++nModified;
            }
        }


        // Pass 2. Do other side of zone
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
        {
            const label zoneFacei = fZone.whichFace(facei);

            if (zoneFacei != -1)
            {
                if (!fZone.flipMap()[zoneFacei])
                {
                    // Use neighbour side of face
                    modifyOrAddFace
                    (
                        meshMod,
                        mesh.faces()[facei].reverseFace(),  // modified face
                        facei,                          // label of face
                        mesh.faceNeighbour()[facei],    // owner
                        true,                           // face flip
                        newSlavePatchi,                 // patch for face
                        fZone.index(),                  // zone for face
                        true,                           // face flip in zone
                        modifiedFace                    // modify or add
                    );
                }
                else
                {
                    // Use owner side of face
                    modifyOrAddFace
                    (
                        meshMod,
                        mesh.faces()[facei],    // modified face
                        facei,                  // label of face
                        mesh.faceOwner()[facei],// owner
                        false,                  // face flip
                        newSlavePatchi,         // patch for face
                        fZone.index(),          // zone for face
                        true,                   // face flip in zone
                        modifiedFace            // modify or add status
                    );
                }
            }
        }


        // Modify any boundary faces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~

        // Normal boundary:
        // - move to new patch. Might already be back-to-back baffle
        // you want to add cyclic to. Do warn though.
        //
        // Processor boundary:
        // - do not move to cyclic
        // - add normal patches though.

        // For warning once per patch.
        labelHashSet patchWarned;

        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];
            const bool isCoupled = pp.coupled();

            if
            (
                isCoupled
             && (
                    pbm[newMasterPatchi].coupled()
                 || pbm[newSlavePatchi].coupled()
                )
            )
            {
                // Do not allow coupled faces to be moved to different
                // coupled patches.
            }
            else if (isCoupled || !internalFacesOnly)
            {
                forAll(pp, i)
                {
                    label facei = pp.start()+i;

                    label zoneFacei = fZone.whichFace(facei);

                    if (zoneFacei != -1)
                    {
                        if (!isCoupled && patchWarned.insert(patchi))
                        {
                            WarningInFunction
                                << "Found boundary face (patch " << pp.name()
                                << ") in faceZone " << fZone.name()
                                << " to convert to baffle patches "
                                << pbm[newMasterPatchi].name() << "/"
                                << pbm[newSlavePatchi].name() << nl
                                << "    Set internalFacesOnly to true in the"
                                << " createBaffles control dictionary if you"
                                << " don't wish to convert boundary faces."
                                << endl;
                        }

                        modifyOrAddFace
                        (
                            meshMod,
                            mesh.faces()[facei],        // modified face
                            facei,                      // label of face
                            mesh.faceOwner()[facei],    // owner
                            false,                      // face flip
                            fZone.flipMap()[zoneFacei]
                          ? newSlavePatchi
                          : newMasterPatchi,            // patch for face
                            fZone.index(),              // zone for face
                            fZone.flipMap()[zoneFacei], // face flip in zone
                            modifiedFace                // modify or add
                        );

                        ++nModified;
                    }
                }
            }
        }
    }
}


void createBaffles
(
    fvMesh& mesh,
    const List<word>& baffleNames,
    const Map<word>& cellZoneToRegionName
)
{
    const word meshInstance = mesh.pointsInstance();

    dictionary baffles;
    HashSet<word> bafflePatches;

    // Create baffles
    // ~~~~~~~~~~~~~~
    // Is done in multiple steps
    // - create patches with 'calculated' patchFields
    // - move faces into these patches
    // - change the patchFields to the wanted type
    // This order is done so e.g. fixedJump works:
    // - you cannot create patchfields at the same time as patches since
    //   they do an evaluate upon construction
    // - you want to create the patchField only after you have faces
    //   so you don't get the 'create-from-nothing' mapping problem.


    // Pass 1: add patches
    // ~~~~~~~~~~~~~~~~~~~

    forAll(baffleNames, zoneI)
    {
        word groupName = baffleNames[zoneI];

        dictionary dict;
        dict.add("type", "faceZone", true);
        dict.add("zoneName", groupName, true);

        dictionary patchPairs;

        patchPairs.add("type", "mappedWall", true);
        patchPairs.add("sampleMode", "nearestPatchFace", true);
        dict.add("patchPairs", patchPairs, true);
        baffles.add(groupName, dict, true);

        const word masterName = groupName + "_master";
        const word slaveName = groupName + "_slave";
        bafflePatches.insert(masterName);
        bafflePatches.insert(slaveName);

        word extraGroupName;
        bool addExtraGroup(false);
        if (dict.found("group"))
        {
            extraGroupName = dict.get<word>("group");
            addExtraGroup = true;
        }

        const dictionary& patchSource = dict.subDict("patchPairs");

        word groupNameMaster = groupName;
        word groupNameSlave = groupName;


        dictionary patchDictMaster(patchSource);
        patchDictMaster.set("nFaces", 0);
        patchDictMaster.set("startFace", 0);
        //patchDictMaster.set("coupleGroup", groupName);
        patchDictMaster.set("samplePatch", slaveName);

        dictionary patchDictSlave(patchDictMaster);
        patchDictSlave.set("samplePatch", masterName);

        // Note: This is added for the particular case where we want
        // master and slave in different groupNames
        // (ie 3D thermal baffles)

        const bool sameGroup =
            patchSource.getOrDefault("sameGroup", true);

        if (!sameGroup)
        {
            groupNameMaster = groupName + "Group_master";
            groupNameSlave = groupName + "Group_slave";
            patchDictMaster.set("coupleGroup", groupNameMaster);
            patchDictSlave.set("coupleGroup", groupNameSlave);
        }

        const label zoneID = mesh.faceZones().findZoneID(groupName); 
        const faceZone& fZone = mesh.faceZones()[zoneID]; 

        label currentCellZoneID
        (
            mesh.cellZones().whichZone(mesh.faceOwner()[fZone[0]])
        );

        word currentRegionName
        (
            cellZoneToRegionName[currentCellZoneID]
        );

        addBafflePatch
        (
            mesh,
            masterName,
            slaveName,
            groupNameMaster,
            extraGroupName,
            patchDictMaster,
            currentRegionName,
            addExtraGroup
        );

        addBafflePatch
        (
            mesh,
            slaveName,
            masterName,
            groupNameSlave,
            extraGroupName,
            patchDictSlave,
            currentRegionName,
            addExtraGroup
        );
    }

    // Make sure patches and zoneFaces are synchronised across couples
    mesh.boundaryMesh().checkParallelSync(true);
    mesh.faceZones().checkParallelSync(true);



    // Mesh change container
    polyTopoChange meshMod(mesh);

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();


    // Do the actual changes. Note: 
    // - loop in incrementing face order (not necessary if faceZone ordered). 
    //   Preserves any existing ordering on patch faces. 
    // - two passes, do non-flip faces first and flip faces second. This 
    //   guarantees that when e.g. creating a cyclic all faces from one 
    //   side come first and faces from the other side next. 
 
    // Whether first use of face (modify) or consecutive (add) 
    bitSet modifiedFace(mesh.nFaces()); 
    label nModified = 0; 
 
    forAll(baffles.toc(), baffleI)
    { 
        const dictionary& dict = baffles.subDict(baffles.toc()[baffleI]);
        const word& groupName = dict.dictName(); 
 
        const label zoneID = mesh.faceZones().findZoneID(groupName); 
        const faceZone& fZone = mesh.faceZones()[zoneID]; 
 
        DynamicList<label> newMasterPatches; 
        DynamicList<label> newSlavePatches; 
 
        const word masterName(groupName + "_master"); 
        newMasterPatches.push_back(pbm.findPatchID(masterName)); 

        const word slaveName(groupName + "_slave"); 
        newSlavePatches.push_back(pbm.findPatchID(slaveName)); 
 
        createFaces 
        ( 
            true,
            mesh, 
            fZone, 
            newMasterPatches, 
            newSlavePatches, 
            meshMod, 
            modifiedFace, 
            nModified 
        ); 
    }
 
 
    Info<< "Converted " << returnReduce(nModified, sumOp<label>()) 
        << " faces into boundary faces in patches " 
        << bafflePatches.sortedToc() << nl << endl;

    // Change the mesh. Change points directly (no inflation).
    meshMod.changeMesh(mesh, false);

    // Remove any now zero-sized patches
    filterPatches(mesh, bafflePatches);

    mesh.setInstance(meshInstance);
}

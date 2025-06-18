/*---------------------------------------------------------------------------*\
|             setAlphaSurface | Copyright (C) 2019 Dezhi Dai             |
-------------------------------------------------------------------------------
    Cork   | Copyright (C) 2016 Gilbert Bernstein
    libigl | Copyright (C) 2018 Alec Jacobson, Daniele Panozzo and others.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    setAlphaSurface

Description
    Initialize the alpha (VOF) field with an input shape surface mesh
    in .stl format. The previous fluid shape will be retained.

    Reference:
        \verbatim
            Dai, Dezhi (2019).
            setAlphaSurface
            Mendeley Data, V2
            doi 10.17632/wg4sx7sc57.2
            url http://dx.doi.org/10.17632/wg4sx7sc57.2
        \endverbatim

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "triSurface.H"
#include "Fstream.H"
#include "triFace.H"
#include "triFaceList.H"
#include "topoSetSource.H"
#include "cellZoneSet.H"

#include "readSTL.h"
#include "remove_duplicate_vertices.h"
#include "mesh_boolean.h"
#include "readSTL.h"
#include "volume.h"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

triSurface surfaceAdd
(
	const triSurface& surface1,
	const triSurface& surface2
)
{
    const pointField& points1 = surface1.points();
    const pointField& points2 = surface2.points();

    // Final surface
    triSurface combinedSurf;

    // Make new storage
    List<labelledTri> facesAll(surface1.size() + surface2.size());
    vectorField pointsAll(points1.size() + points2.size());


    label pointi = 0;
    // Copy points1 into pointsAll
    for (const auto& pt : points1)
    {
        pointsAll[pointi++] = pt;
    }
    // Add surface2 points
    for (const auto& pt : points2)
    {
        pointsAll[pointi++] = pt;
    }


    label trianglei = 0;

    // Determine map for both regions
    label nNewPatches = 0;
    labelList patch1Map(surface1.patches().size());
    labelList patch2Map(surface2.patches().size());

    patch1Map = identity(surface1.patches().size());
    patch2Map = identity(surface2.patches().size(), patch1Map.size());

    nNewPatches = surface1.patches().size()+surface2.patches().size();

    // Copy triangles1 into trianglesAll
    for (const labelledTri& tri : surface1)
    {
        labelledTri& destTri = facesAll[trianglei++];

        destTri.triFace::operator=(tri);
        destTri.region() = patch1Map[tri.region()];
    }

    // Add (renumbered) surface2 triangles
    for (const labelledTri& tri : surface2)
    {
        labelledTri& destTri = facesAll[trianglei++];
        destTri[0] = tri[0] + points1.size();
        destTri[1] = tri[1] + points1.size();
        destTri[2] = tri[2] + points1.size();
        destTri.region() = patch2Map[tri.region()];
    }

    geometricSurfacePatchList newPatches(nNewPatches);
    forAll(surface1.patches(), patchi)
    {
        newPatches[patch1Map[patchi]] = surface1.patches()[patchi];
    }
    forAll(surface2.patches(), patchi)
    {
        newPatches[patch2Map[patchi]] = surface2.patches()[patchi];
    }

    Info<< "New patches:" << nl;
    forAll(newPatches, patchi)
    {
        Info<< "    " << patchi << '\t' << newPatches[patchi].name() << nl;
    }
    Info<< endl;


    // Construct new surface mesh
    combinedSurf = triSurface(facesAll, newPatches, pointsAll);

    // Merge all common points and do some checks
    combinedSurf.cleanup(true);

    Info<< "Merged surface:" << endl;

    combinedSurf.writeStats(Info);

    Info<< endl;

    return combinedSurf;
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Initialize the alpha (VOF) field with an input shape surface mesh "
        "in .stl format. The previous fluid shape will be retained.\n"
    );

    argList::noFunctionObjects();   // Don't use function objects

    argList::addBoolOption
    (
        "skipTopoSet",
        "skp creating celZones"
    );

    argList::addBoolOption
    (
        "writeBoundaryFields",
        "Write separate boundary field data to time/boundaryFields directory"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const bool skipTopoSet = args.found("skipTopoSet");
    const bool writeBoundaryFields(args.found("writeBoundaryFields"));

    const word dictName("setAlphaSurfaceDict");
    #include "setSystemMeshDictionaryIO.H"

    IOdictionary setAlphaSurfaceDict(dictIO);

    Info<< "Reading " << setAlphaSurfaceDict.name() << nl << endl;

    const word fieldName
    (
        setAlphaSurfaceDict.getOrDefault<word>("field", "alpha.water")
    );
    fileName stlDir(setAlphaSurfaceDict.get<fileName>("STLDir"));
    List<word> stlFileNames(setAlphaSurfaceDict.get<List<word>>("STLNames"));
    point outsidePoint(setAlphaSurfaceDict.get<point>("outsidePoint"));
    pointField outsidePoints(1, outsidePoint);

    fileName caseDir(runTime.rootPath()/runTime.globalCaseName());
    if (isDir(caseDir/stlDir))
    {
        stlDir = caseDir/stlDir;
    }

    fileName unifiedSTL(stlDir/"unified" + ".stl");

    if (Pstream::master())
    {
        Foam::rm(unifiedSTL);

        if (stlFileNames.size() > 1)
        {
            fileName stl1(stlDir/stlFileNames[0]);
            fileName stl2(stlDir/stlFileNames[1]);
        
            const triSurface surf1(stl1, -1);
            const triSurface surf2(stl2, -1);
            triSurface combinedSurf = surfaceAdd(surf2, surf1);

            forAll (stlFileNames, i)
            {
                if ( i > 1)
                {
                    fileName stl3(stlDir/stlFileNames[i]);

                    const triSurface surf3(stl3, -1);
                    combinedSurf = surfaceAdd(surf3, combinedSurf);
                }
            }

            combinedSurf.write(unifiedSTL, false);
        }
    }

    if (stlFileNames.size() == 1)
    {
        unifiedSTL = stlDir/stlFileNames[0];
    }

    List<dictionary> actionEntries;
    List<word> alphaSetCut;
    word alphaSetCell("alphaSetCell");

    forAll (stlFileNames, i)
    {
        dictionary dict;
        word name("alphaSetCut-" + stlFileNames[i]);
        alphaSetCut.append(name);

        dict.add("name", name);
        dict.add("type", "cellSet");
        dict.add("action", "new");
        dict.add("source", "surfaceToCell");

        dictionary sourceDict;
        sourceDict.add("file", stlDir/stlFileNames[i]);
        sourceDict.add("useSurfaceOrientation", "false");
        sourceDict.add("outsidePoints", outsidePoints);
        sourceDict.add("includeCut", "true");
        sourceDict.add("includeInside", "false");
        sourceDict.add("includeOutside", "false");
        sourceDict.add("nearDistance", -1);
        sourceDict.add("curvature", -100);
        dict.add("sourceInfo", sourceDict);

        actionEntries.append(dict);

        dict.clear();
        dict.add("name", name);
        dict.add("type", "cellZoneSet");
        dict.add("action", "new");
        dict.add("source", "setToCellZone");
        sourceDict.clear();
        sourceDict.add("set", name);
        dict.add("sourceInfo", sourceDict);

        actionEntries.append(dict);
    }

    {
        dictionary dict;

        dict.add("name", alphaSetCell);
        dict.add("type", "cellSet");
        dict.add("action", "new");
        dict.add("source", "surfaceToCell");

        dictionary sourceDict;
        sourceDict.add("file", unifiedSTL);
        sourceDict.add("useSurfaceOrientation", "false");
        sourceDict.add("outsidePoints", outsidePoints);
        sourceDict.add("includeCut", "false");
        sourceDict.add("includeInside", "true");
        sourceDict.add("includeOutside", "false");
        sourceDict.add("nearDistance", -1);
        sourceDict.add("curvature", -100);
        dict.add("sourceInfo", sourceDict);

        actionEntries.append(dict);

        dict.clear();
        dict.add("name", alphaSetCell);
        dict.add("type", "cellZoneSet");
        dict.add("action", "new");
        dict.add("source", "setToCellZone");
        sourceDict.clear();
        sourceDict.add("set", alphaSetCell);
        dict.add("sourceInfo", sourceDict);

        actionEntries.append(dict);
    }


    actionEntries.setSize(actionEntries.size());
    alphaSetCut.setSize(alphaSetCut.size());

    if (!skipTopoSet)
    {
        // Execute all actions
        forAll (actionEntries, actionI)
        {
            const dictionary& dict = actionEntries[actionI];
            if (dict.empty())
            {
                continue;
            }
            const word setName(dict.get<word>("name"));
            const word setType(dict.get<word>("type"));

            const topoSetSource::setAction action =
                topoSetSource::actionNames.get("action", dict);

            autoPtr<topoSet> currentSet;

            currentSet = topoSet::New(setType, mesh, setName, 16384);

            const word sourceType(dict.get<word>("source"));

            Info<< "    Applying source " << sourceType << endl;
            autoPtr<topoSetSource> source = topoSetSource::New
            (
                sourceType,
                mesh,
                dict.optionalSubDict("sourceInfo")
            );

            source().applyToSet(action, currentSet());
            // Synchronize for coupled patches.
            currentSet().sync(mesh);
            if (!currentSet().write())
            {
                WarningInFunction
                    << "Failed writing set "
                    << currentSet().objectPath() << endl;
            }
            fileHandler().flush();

            if (currentSet)
            {
                if (currentSet().type() == "cellZoneSet")
                {
                    Info<< "    "
                        << currentSet().type() << ' '
                        << currentSet().name() << " now size "
                        << returnReduce(currentSet().size(), sumOp<label>())
                        << nl << endl;
                }
            }
        }

        Info << nl << "    Wet cellZones have been created..." << nl << endl;

        runTime.printExecutionTime(Info);
    }

    if (stlFileNames.size() > 1)
    {
        Foam::rm(stlDir/"unified.stl");
    }

    Info << nl << "Reading field " << fieldName << nl << endl;
    
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            "0",
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimVolume/dimTime, Zero)
    );

    volScalarField alpha1
    (
        IOobject
        (
            fieldName,
            "0",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    Eigen::MatrixXd temp_V, V, N;
    Eigen::MatrixXi temp_F, VI, VJ, F;

    Info << "Initializing field " << fieldName << nl << endl;
    const scalarField& cellVolumes = mesh.V();
    const cellList& cellFaces = mesh.cells();
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    const labelList& own = mesh.faceOwner();

    forAll(alphaSetCut ,i)
    {
        fileName stlName = stlDir/stlFileNames[i];

        bool success = igl::readSTL(stlName, temp_V, temp_F, N);
        if (success)
        {
            igl::remove_duplicate_vertices
            (
                temp_V, 
                temp_F, 
                1e-6, 
                V, 
                VI, 
                VJ, 
                F
            );
        }

        const label cutCellZoneID = 
            mesh.cellZones().findZoneID(alphaSetCut[i]);

        const labelList& cutCellLabels = mesh.cellZones()[cutCellZoneID];

        forAll(cutCellLabels, celli)
        {
            const label& cellI = cutCellLabels[celli];

            scalar fluidVol(0.0);

            // Localize cell point labels
            const cell& currentCell(cellFaces[cellI]);
            const labelList& globalCellPtLabels(currentCell.labels(faces));
            pointField localPts;
            forAll(globalCellPtLabels, pointI)
            {
                localPts.append(points[globalCellPtLabels[pointI]]);
            }


            // Localize cell face labels
            labelListList localFaces;
            forAll(currentCell, faceI)
            {
                const face& fa = faces[currentCell[faceI]];

                localPts.append(fa.centre(points));

                edgeList edges = fa.edges();

                forAll(edges, edgeI)
                {
                    labelList localTriface;
                    label c0 = localPts.size() - 1;
                    localTriface.append(c0);
                    label c1 = edges[edgeI].start();
                    label c2 = edges[edgeI].end();

                    forAll(globalCellPtLabels, pointI)
                    {
                        if(c1 == globalCellPtLabels[pointI])
                        {
                            c1 = pointI;
                            break;
                        }
                    }

                    forAll(globalCellPtLabels, pointI)
                    {
                        if(c2 == globalCellPtLabels[pointI])
                        {
                            c2 = pointI;
                            break;
                        }
                    }

                    if(cellI == own[currentCell[faceI]])
                    {
                        localTriface.append(c1);
                        localTriface.append(c2);
                    }
                    else
                    {
                        localTriface.append(c2);
                        localTriface.append(c1);
                    }

                    localFaces.append(localTriface);
                }
            }

            Eigen::MatrixXd VC(localPts.size(), 3);
            forAll(localPts, pointI)
            {
                VC(pointI, 0) = localPts[pointI].x();
                VC(pointI, 1) = localPts[pointI].y();
                VC(pointI, 2) = localPts[pointI].z();
            }

            Eigen::MatrixXi FC(localFaces.size(), 3);
            forAll(localFaces, faceI)
            {
                const labelList& fa(localFaces[faceI]);
                FC(faceI, 0) = fa[0];
                FC(faceI, 1) = fa[1];
                FC(faceI, 2) = fa[2];
            }


            // Calculate intersection of the shape surface mesh and cellI
            Eigen::MatrixXd V_bool;
            Eigen::MatrixXi F_bool;

            igl::copyleft::cork::mesh_boolean
            (
                VC, 
                FC, 
                V, 
                F, 
                igl::MESH_BOOLEAN_TYPE_INTERSECT, 
                V_bool, 
                F_bool
            );

            // Calculate intersected volume
            if (F_bool.rows() > 0)
            {
                Eigen::MatrixXd V2(V_bool.rows() + 1, V_bool.cols());
                V2.topRows(V_bool.rows()) = V_bool;
                V2.bottomRows(1).setZero();
                Eigen::MatrixXi T(F_bool.rows(), 4);
                T.leftCols(3) = F_bool;
                T.rightCols(1).setConstant(V_bool.rows());
                Eigen::VectorXd vol;
                igl::volume(V2, T, vol);

                fluidVol = mag(vol.sum());

                scalar alphaFluid(fluidVol / cellVolumes[cellI]);
                alphaFluid = min
                (
                    scalar(1.0), 
                    max(scalar(0.0), alphaFluid)
                );

                alpha1.ref()[cellI] = alphaFluid;
            }
        }
    }

    const label insideCellZoneID = 
        mesh.cellZones().findZoneID(alphaSetCell);

    const labelList& insideCellLabels = mesh.cellZones()[insideCellZoneID];

    forAll(insideCellLabels, celli)
    {
        const label& cellI = insideCellLabels[celli];

        alpha1.ref()[cellI] = 1.0;
    }
    
    alpha1.correctBoundaryConditions();

    Info << "Writing field " << fieldName << nl << endl;

    alpha1.writeOpt(IOobject::AUTO_WRITE);
    alpha1.write();

    if ( writeBoundaryFields && !mesh.thisDb().time().processorCase())
    {
        alpha1.writeBoundaryField(true);
    }

    runTime.printExecutionTime(Info);

    return 0;
}

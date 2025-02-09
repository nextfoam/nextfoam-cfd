/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Test-boundaryMeshEntries

Description
    Find and print "boundary" information

\*---------------------------------------------------------------------------*/

#include "argList.H"

#include "fileName.H"
#include "IOstreams.H"
#include "polyMesh.H"
#include "polyBoundaryMeshEntries.H"
#include "OSspecific.H"
#include "Time.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noFunctionObjects();  // Disallow function objects
    argList::addVerboseOption("additional verbosity");

    #include "setRootCase.H"
    #include "createTime.H"

    fileName coherentInst;
    coherentInst =
    (
        runTime.findInstance
        (
            polyMesh::meshSubDir,
            "coherent",
            IOobject::READ_IF_PRESENT
        )
    );

    // Unfortunately with READ_IF_PRESENT, cannot tell if the file
    // was actually found or not

    Info<< "check: " << (coherentInst/polyMesh::meshSubDir/"coherent") << nl;

    if (!Foam::isFile(coherentInst/polyMesh::meshSubDir/"coherent"))
    {
        coherentInst.clear();
    }

    Info<< "found coherent: " << coherentInst << nl;

    PtrList<entry> entries;

    if (!coherentInst.empty())
    {
        IOdictionary coherent
        (
            IOobject
            (
                "coherent",
                coherentInst,
                polyMesh::meshSubDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            )
        );

        ITstream& is = coherent.lookup("boundary");
        is >> entries;
    }

    Info<< "entries: " << entries << nl;

    Info<< "type:  " << polyBoundaryMeshEntries::types(entries) << nl;
    Info<< "start: " << polyBoundaryMeshEntries::patchStarts(entries) << nl;
    Info<< "size: "  << polyBoundaryMeshEntries::patchSizes(entries) << nl;
    Info<< nl;

    fileName boundaryInst;
    boundaryInst =
    (
        runTime.findInstance
        (
            polyMesh::meshSubDir,
            "boundary",
            IOobject::MUST_READ
        )
    );

    Info<< "found boundary: " << boundaryInst << nl;

    polyBoundaryMeshEntries pbm
    (
        IOobject
        (
            "boundary",
            boundaryInst,
            polyMesh::meshSubDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    Info<< "type:  " << flatOutput(pbm.types()) << nl;
    Info<< "start: " << flatOutput(pbm.patchStarts()) << nl;
    Info<< "size: "  << flatOutput(pbm.patchSizes()) << nl;

    pbm.writeEntry("boundary", Info);

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //

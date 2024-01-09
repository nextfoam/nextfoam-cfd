/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "GeometricField.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "fvPatchField.H"
#include "fvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::multiRegionSolutionControl::storePrevIter() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> vField;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> sField;

    forAll(regions_, regioni)
    {
        const fvMesh* fvmesh(regions_[regioni]);

        HashTable<const vField*>
            vflds(fvmesh->objectRegistry::lookupClass<vField>());

        forAllIters(vflds, iter)
        {
            const vField& vfld = *iter();

            const word& vfldName = vfld.name();

            if
            (
                (
                    (vfldName.find("PrevIter") == std::string::npos)
                 && (vfldName.find("_0") == std::string::npos)
                 && fvmesh->relaxField(vfldName)
                )
             || (vfldName == "U")
            )
            {
                DebugInfo
                    << algorithmName_ << ": storing previous iter for "
                    << vfldName << endl;

                vfld.storePrevIter();
            }
        }

        HashTable<const sField*>
            sflds(fvmesh->objectRegistry::lookupClass<sField>());

        forAllIters(sflds, iter)
        {
            const sField& sfld = *iter();

            const word& sfldName = sfld.name();

            if (sfldName == "phi")
            {
                 DebugInfo
                    << algorithmName_ << ": storing previous iter for "
                    << sfldName << endl;

                sfld.storePrevIter();
            }
        }
    }
}


// ************************************************************************* //

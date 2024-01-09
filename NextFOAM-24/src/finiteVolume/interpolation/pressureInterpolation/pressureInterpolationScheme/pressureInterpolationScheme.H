/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::pressureInterpolationScheme

Description
    Abstract base class for pressure interpolation schemes.

SourceFiles
    pressureInterpolationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef pressureInterpolationScheme_H
#define pressureInterpolationScheme_H

#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                 Class pressureInterpolationScheme Declaration
\*---------------------------------------------------------------------------*/

class pressureInterpolationScheme
:
    public refCount
{
protected:

    // Protected Data

        //- Hold reference to mesh
        const fvMesh& mesh_;

        //- Reference to the fvOptions
        const fv::optionList& fvOptions_;

        //- Has active porous media zone?
        bool hasActivePorousZone_;


private:

    // Private Member Functions

        //- No copy construct
        pressureInterpolationScheme
        (
            const pressureInterpolationScheme&
        ) = delete;

        //- No copy assignment
        void operator=(const pressureInterpolationScheme&) = delete;


public:

    //- Runtime type information
    TypeName("pressureInterpolationScheme");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            pressureInterpolationScheme,
            Mesh,
            (
                const fvMesh& mesh,
                Istream& schemeData
            ),
            (mesh, schemeData)
        );


    // Constructors

        //- Construct from mesh
        pressureInterpolationScheme(const fvMesh& mesh)
        :
            mesh_(mesh),
            fvOptions_
            (
                mesh_.foundObject<fv::optionList>("fvOptions")
              ? mesh_.lookupObject<fv::optionList>("fvOptions")
              : fv::options::New(mesh_)
            ),
            hasActivePorousZone_(fvOptions_.hasActivePorousZone())
        {}


    // Selectors

        //- Return new tmp interpolation scheme
        static tmp<pressureInterpolationScheme> New
        (
            const fvMesh& mesh,
            Istream& schemeData
        );


    //- Destructor
    virtual ~pressureInterpolationScheme() = default;


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return the face-interpolate of the given cell field
        //  with the given weighting factors
        static tmp<surfaceScalarField>
        interpolate
        (
            const volScalarField&,
            const tmp<surfaceScalarField>&
        );

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> weights
        (
            const volScalarField&
        ) const = 0;

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return false;
        }

        //- Return the explicit correction to the face-interpolate
        //  for the given field
        virtual tmp<surfaceScalarField>
        correction(const volScalarField&) const
        {
            return tmp<surfaceScalarField>
            (
                nullptr
            );
        }

        //- Return the face-interpolate of the given cell field
        //  with explicit correction
        virtual tmp<surfaceScalarField>
        interpolate(const volScalarField&) const;

        //- Return the face-interpolate of the given tmp cell field
        //  with explicit correction
        tmp<surfaceScalarField> interpolate(const tmp<volScalarField>&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makePressureInterpolationScheme(SS)                                    \
                                                                               \
defineTypeNameAndDebug(SS, 0);                                                 \
pressureInterpolationScheme::addMeshConstructorToTable<SS>                     \
    add##SS##MeshConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2018 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    Foam::MRFZoneList

Description
    List container for MRF zomes

SourceFiles
    MRFZoneList.C

\*---------------------------------------------------------------------------*/

#ifndef MRFZoneList_H
#define MRFZoneList_H

#include "fvMesh.H"
#include "dictionary.H"
#include "fvMatricesFwd.H"
#include "MRFZone.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class MRFZoneList;
Ostream& operator<<(Ostream& os, const MRFZoneList& models);

/*---------------------------------------------------------------------------*\
                         Class MRFZoneList Declaration
\*---------------------------------------------------------------------------*/

class MRFZoneList
:
    public PtrList<MRFZone>
{
    // Private Member Functions

        //- No copy construct
        MRFZoneList(const MRFZoneList&) = delete;

        //- No copy assignment
        void operator=(const MRFZoneList&) = delete;


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        mutable volVectorField coriolisForce_;   // by Gill

        mutable surfaceScalarField volumeAcceleration_;   // by Gill


public:

    //- Constructor
    MRFZoneList(const fvMesh& mesh, const dictionary& dict);

    //- Destructor
    ~MRFZoneList() = default;


    // Member Functions

        //- Return active status
        bool active(const bool warn = false) const;

        //- Reset the source list
        void reset(const dictionary& dict);

        //- Return the MRF with a given name
        const MRFZone& getFromName(const word& name) const;

        //- Add the frame acceleration
        void addAcceleration
        (
            const volVectorField& U,
            volVectorField& ddtU
        ) const;

        //- Add the frame acceleration contribution to the momentum equation
        void addAcceleration(fvVectorMatrix& UEqn) const;

        //- Add the frame acceleration contribution to the momentum equation
        void addAcceleration
        (
            const volScalarField& rho,
            fvVectorMatrix& UEqn
        ) const;

        //- Return the frame acceleration
        tmp<volVectorField> DDt
        (
            const volVectorField& U
        ) const;

        //- Return the frame acceleration
        tmp<volVectorField> DDt
        (
            const volScalarField& rho,
            const volVectorField& U
        ) const;

        //- Return the MRF absolute flux
        tmp<surfaceScalarField> phi() const;

        tmp<volVectorField> coriolisForce
        (
            const volVectorField& U
        ) const; // by Gill

        tmp<volVectorField> coriolisForce
        (
            const volScalarField& rho,
            const volVectorField& U
        ) const; // by Gill

        tmp<surfaceScalarField> volumeAcceleration() const; // by Gill

        //- Make the given absolute velocity relative within the MRF region
        void makeRelative(volVectorField& U) const;

        //- Make the given absolute flux relative within the MRF region
        void makeRelative(surfaceScalarField& phi) const;

        //- Return the given absolute flux relative within the MRF region
        tmp<surfaceScalarField> relative
        (
            const tmp<surfaceScalarField>& phi
        ) const;

        //- Return the given absolute boundary flux relative within
        //  the MRF region
        tmp<FieldField<fvsPatchField, scalar>> relative
        (
            const tmp<FieldField<fvsPatchField, scalar>>& tphi
        ) const;

        //- Return the given absolute patch flux relative within
        //  the MRF region
        tmp<Field<scalar>> relative
        (
            const tmp<Field<scalar>>& tphi,
            const label patchi
        ) const;

        //- Make the given absolute mass-flux relative within the MRF region
        void makeRelative
        (
            const surfaceScalarField& rho,
            surfaceScalarField& phi
        ) const;

        void addCoriolisFlux
        (
            const surfaceScalarField& rAUf,
            surfaceScalarField& phi
        ) const; // by Gill

        void addCoriolisFlux
        (
            const surfaceScalarField& rhof,
            const surfaceScalarField& rhorAUf,
            surfaceScalarField& phi
        ) const; // by Gill

        //- Make the given relative velocity absolute within the MRF region
        void makeAbsolute(volVectorField& U) const;

        //- Make the given relative flux absolute within the MRF region
        void makeAbsolute(surfaceScalarField& phi) const;

        //- Return the given relative flux absolute within the MRF region
        tmp<surfaceScalarField> absolute
        (
            const tmp<surfaceScalarField>& phi
        ) const;

        //- Make the given relative mass-flux absolute within the MRF region
        void makeAbsolute
        (
            const surfaceScalarField& rho,
            surfaceScalarField& phi
        ) const;

        //- Correct the boundary velocity for the rotation of the MRF region
        void correctBoundaryVelocity(volVectorField& U) const;

        void updateForce // by Gill
        (
            const surfaceScalarField& phi, 
            const volVectorField& U
        ) const;

        void updateForce // by Gill
        (
            const surfaceScalarField& rhof, 
            const surfaceScalarField& phi, 
            const volVectorField& U
        ) const;

        //- Correct the boundary flux for the rotation of the MRF region
        void correctBoundaryFlux
        (
            const volVectorField& U,
            surfaceScalarField& phi
        ) const;

        //- Filter-out the MRF region contribution from the given field
        // setting the corresponding values to zero
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> zeroFilter
        (
            const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tphi
        ) const;

        //- Update MRFZone faces if the mesh topology changes
        void update();


        // I-O

            //- Read dictionary
            bool read(const dictionary& dict);

            //- Write data to Ostream
            bool writeData(Ostream& os) const;

            //- Ostream operator
            friend Ostream& operator<<
            (
                Ostream& os,
                const MRFZoneList& models
            );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MRFZoneListTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

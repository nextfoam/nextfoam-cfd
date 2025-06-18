/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::basicWindSpeedFvPatchVectorField

Group
    grpRASBoundaryConditions grpInletBoundaryConditions

Description
    Inlet conditions with boundary layer profile of atmospheric wind speed
    according to surface roughness prepared for design pressure analysis of
    buildings. (by Building structure standards of the Ministry of
    Land, Infrastructure and Transport of Korea, 2016-317, 82 p)

    \f[
        V_{H} = V_{0} K_{zr} K_{zt} I_{w}
    \f]

    where
    \varatable
        V_{H}   | Designed wind speed                                                   [m/s]
        V_{0}   | Fundamental wind speed                                                [m/s]
        K_{zr}  | Distribution coefficient at reference height H as wind speed altitude [-]
        K_{zt}  | Topography coefficient as rougness factor                             [-]
        I_{w}   | topography coefficient                                                [-]

Usage
    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        // Mandatory entries (unmodifiable)
        type            basicWindSpeed;

        flowDir         uniform (1 0 0); 
        zDir            uniform (0 0 1); 

        V0              28;

        // Roughness factor correspond to 1:A, 2:B, 3:C, 4:D.
        Roughness       1; 
        Kzt             1.0;
        Iw              1.0;           

        // Mandatory/Optional (inherited) entries (unmodifiable)
        ...
    }
    \endverbatim

    The inherited entries are elaborated in:
     - \link fixedValueFvPatchField.H \endlink


SourceFiles
    basicWindSpeedFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef basicWindSpeedFvPatchVectorField_H
#define basicWindSpeedFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "Function1.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class basicWindSpeedFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class basicWindSpeedFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

private:

    // Private Data

        //- Minimum coordinate vector of this patch
        const vector ppMin_;

        //- Reference to the time database
        const Time& time_;

        //- Reference to the patch
        const polyPatch& patch_;

        //- Streamwise flow direction
        autoPtr<Function1<vector>> flowDir_;

        //- Direction of the ground-normal coordinate
        autoPtr<Function1<vector>> zDir_;

        //
        scalar V0_;

        //
        scalar Kzt_;

        //
        int Roughness_;

        //
        scalar Iw_;

public:

    //- Runtime type information
    TypeName("basicWindSpeed");


    // Constructors

        //- Construct from patch and internal field
        basicWindSpeedFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        basicWindSpeedFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- basicWindSpeedFvPatchVectorField onto a new patch
        basicWindSpeedFvPatchVectorField
        (
            const basicWindSpeedFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new basicWindSpeedFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        basicWindSpeedFvPatchVectorField
        (
            const basicWindSpeedFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new basicWindSpeedFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Access
            // Return Wind Velocity;
            vectorField U(const fvPatch& patch) const;

            // Return Flow Direction;
            vectorField flowDir(const fvPatch& patch) const;

            // Return Z direction;
            vector zDir() const;

            // Return Kzr
            scalar Kzr(scalar& z) const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();


        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

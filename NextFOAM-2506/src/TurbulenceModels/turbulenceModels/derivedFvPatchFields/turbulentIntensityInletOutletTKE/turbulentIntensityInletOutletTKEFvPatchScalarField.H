/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016, 2019 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::turbulentIntensityInletOutletTKEFvPatchScalarField

Group
    grpInletBoundaryConditions

Description
    This boundary condition provides a turbulent kinetic energy condition,
    based on user-supplied turbulence intensity, defined as a fraction of the
    mean velocity:

        \f[
            k_p = 1.5 (I |U|)^2
        \f]

    where

    \vartable
        k_p     | kinetic energy at the patch
        I       | turbulence intensity
        U       | velocity field
    \endvartable

    In the event of reverse flow, a zero-gradient condition is applied.

Usage
    \table
        Property     | Description             | Required    | Default value
        turbIntensity    | fraction of mean field [0-1] | yes    |
        U            | velocity field name     | no          | U
        phi          | flux field name         | no          | phi
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type        turbulentIntensityInletOutletTKE;
        turbIntensity   0.05;           // 5% turbulence
        value       uniform 1;      // placeholder
    }
    \endverbatim

See also
    Foam::inletOutletFvPatchField

SourceFiles
    turbulentIntensityInletOutletTKEFvPatchScalarField.C

Author
    JaeHeung Gill, NEXTfoam Co.,Ltd.

\*---------------------------------------------------------------------------*/

#ifndef turbulentIntensityInletOutletTKEFvPatchScalarField_H
#define turbulentIntensityInletOutletTKEFvPatchScalarField_H

#include "inletOutletFvPatchFields.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class turbulentIntensityInletOutletTKEFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentIntensityInletOutletTKEFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private data

        //- Turbulent intensity as fraction of mean velocity
        autoPtr<PatchFunction1<scalar>> turbIntensity_;

        //- Name of the velocity field
        word UName_;


public:

    //- Runtime type information
    TypeName("turbulentIntensityInletOutletTKE");


    // Constructors

        //- Construct from patch and internal field
        turbulentIntensityInletOutletTKEFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentIntensityInletOutletTKEFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentIntensityInletOutletTKEFvPatchScalarField onto a new patch
        turbulentIntensityInletOutletTKEFvPatchScalarField
        (
            const turbulentIntensityInletOutletTKEFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentIntensityInletOutletTKEFvPatchScalarField
        (
            const turbulentIntensityInletOutletTKEFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentIntensityInletOutletTKEFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        turbulentIntensityInletOutletTKEFvPatchScalarField
        (
            const turbulentIntensityInletOutletTKEFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentIntensityInletOutletTKEFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

            //- Return the name of the velocity field
            const word& UName() const
            {
                return UName_;
            }

            //- Return reference to the name of the velocity field
            //  to allow adjustment
            word& UName()
            {
                return UName_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

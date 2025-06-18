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
    Foam::inletOutletTotalTemperatureFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition provides a total temperature for inflow boundaries
    and zero gradients for outflow boundaries. A user-specified 'T0' value is 
    applied as the static temperature in incompressible cases for inflow 
    boundaries.

Usage
    \table
        Property     | Description             | Required    | Default value
        U            | velocity field name     | no          | U
        phi          | flux field name         | no          | phi
        psi          | compressibility field name |  no      | thermo:psi
        gamma        | heat capacity ration (Cp/Cv) | yes    |
        T0           | total temperature [K]  | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            inletOutletTotalTemperature;
        U               U;
        phi             phi;
        psi             psi;
        gamma           gamma;
        T0              uniform 0;
        value           uniform 0;
    }
    \endverbatim

Author
    Modified JaeHeung Gill, NEXTfoam Co.,Ltd.

See also
    Foam::inletOutletFvPatchField

SourceFiles
    inletOutletTotalTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef inletOutletTotalTemperatureFvPatchScalarField_H
#define inletOutletTotalTemperatureFvPatchScalarField_H

#include "inletOutletFvPatchFields.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class inletOutletTotalTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class inletOutletTotalTemperatureFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private data

        //- Name of the velocity field
        word UName_;

        //- Name of the compressibility field used to calculate the wave speed
        word psiName_;

        //- Heat capacity ratio
        scalar gamma_;

        //- Total temperature
        autoPtr<PatchFunction1<scalar>> T0_;


public:

    //- Runtime type information
    TypeName("inletOutletTotalTemperature");


    // Constructors

        //- Construct from patch and internal field
        inletOutletTotalTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        inletOutletTotalTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  inletOutletTotalTemperatureFvPatchScalarField onto a new patch
        inletOutletTotalTemperatureFvPatchScalarField
        (
            const inletOutletTotalTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        inletOutletTotalTemperatureFvPatchScalarField
        (
            const inletOutletTotalTemperatureFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        inletOutletTotalTemperatureFvPatchScalarField
        (
            const inletOutletTotalTemperatureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return fvPatchField<scalar>::Clone(*this);
        }

        //- Clone with an internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return fvPatchField<scalar>::Clone(*this, iF);
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

            //- Return the heat capacity ratio
            scalar gamma() const
            {
                return gamma_;
            }

            //- Return reference to the heat capacity ratio to allow adjustment
            scalar& gamma()
            {
                return gamma_;
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

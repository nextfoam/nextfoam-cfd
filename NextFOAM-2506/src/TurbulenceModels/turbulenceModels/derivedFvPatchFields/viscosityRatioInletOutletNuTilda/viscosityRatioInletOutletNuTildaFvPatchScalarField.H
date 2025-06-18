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
    Foam::viscosityRatioInletOutletNuTildaFvPatchScalarField

Group
    grpCmpRASBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition provides a modified turbulent viscosity,
    \f$\tilde{\nu}\f$(nuTilda) inlet/outlet condition based on 
    a specified turbulent viscosity ratio. 
    The patch values are calculated using:

        \f[
            \tilde{\nu}_p = \nu\lambda
        \f]
   
    where

    \vartable
        \tilde{\nu}_p   | patch nuTilda values
        \lambda         | turbulent viscosity ratio (\frac{\tilde{\nu}}{\nu})
    \endvartable

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        viscosityRatio | turbulence viscosity ratio  | yes         |
        phi          | flux field name         | no          | phi
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            viscosityRatioInletOutletNuTilda;
        viscosityRatio  5;
        value           uniform 200;   // placeholder
    }
    \endverbatim

Note
    In the event of reverse flow, a zero-gradient condition is applied

SeeAlso
    Foam::inletOutletFvPatchField

SourceFiles
    viscosityRatioInletOutletNuTildaFvPatchScalarField.C

Author
    JaeHeung Gill, NEXTfoam Co.,Ltd.

\*---------------------------------------------------------------------------*/

#ifndef viscosityRatioInletOutletNuTildaFvPatchField_H
#define viscosityRatioInletOutletNuTildaFvPatchField_H

#include "inletOutletFvPatchFields.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
 Class viscosityRatioInletOutletNuTildaFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class viscosityRatioInletOutletNuTildaFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private data

        //- turbulent viscosity ratio
        autoPtr<PatchFunction1<scalar>> viscosityRatio_;


public:

    //- Runtime type information
    TypeName("viscosityRatioInletOutletNuTilda");


    // Constructors

        //- Construct from patch and internal field
        viscosityRatioInletOutletNuTildaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        viscosityRatioInletOutletNuTildaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  viscosityRatioInletOutletNuTildaFvPatchScalarField onto a new patch
        viscosityRatioInletOutletNuTildaFvPatchScalarField
        (
            const viscosityRatioInletOutletNuTildaFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        viscosityRatioInletOutletNuTildaFvPatchScalarField
        (
            const viscosityRatioInletOutletNuTildaFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new viscosityRatioInletOutletNuTildaFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        viscosityRatioInletOutletNuTildaFvPatchScalarField
        (
            const viscosityRatioInletOutletNuTildaFvPatchScalarField&,
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
                new viscosityRatioInletOutletNuTildaFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

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

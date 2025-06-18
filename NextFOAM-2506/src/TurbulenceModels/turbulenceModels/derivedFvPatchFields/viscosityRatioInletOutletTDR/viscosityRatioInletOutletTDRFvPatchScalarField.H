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
    Foam::viscosityRatioInletOutletTDRFvPatchScalarField

Group
    grpCmpRASBoundaryConditions grpInletBoundaryConditions

Description
    This boundary condition provides a turbulence dissipation rate, 
    \f$\epsilon\f$(epsilon) or a turbulence specific dissipation rate,
    \f$\omega\f$(omega) inlet/outlet condition based on a turbulent viscosity 
    ratio. The patch values are calculated using:

        \f[
            \epsilon_p = C_{\mu}\frac{k^2}{\nu}\lambda^{-1}
        \f]

        \f[
            \omega_p = \frac{k}{\nu}\lambda^{-1}
        \f]
  
    where

    \vartable
        \epsilon_p  | patch epsilon values
        \omega_p    | patch omega values
        C_{\mu}     | Model coefficient, set to 0.09
        \lambda     | turbulent viscosity ratio (\frac{\nu_t}{\nu})
    \endvartable

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        viscosityRatio | turbulent viscosity ratio  | yes    |
        phi          | flux field name         | no          | phi
        k            | turbulence kinetic energy field name | no | k
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            viscosityRatioInletOutletTDR;
        viscosityRatio  5;
        value           uniform 200;   // placeholder
    }
    \endverbatim

Note
    In the event of reverse flow, a zero-gradient condition is applied

SeeAlso
    Foam::inletOutletFvPatchField

SourceFiles
    viscosityRatioInletOutletTDRFvPatchScalarField.C

Author
    JaeHeung Gill, NEXTfoam Co.,Ltd.

\*---------------------------------------------------------------------------*/

#ifndef viscosityRatioInletOutletTDRFvPatchField_H
#define viscosityRatioInletOutletTDRFvPatchField_H

#include "inletOutletFvPatchFields.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
 Class viscosityRatioInletOutletTDRFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class viscosityRatioInletOutletTDRFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private data

        //- turbulent viscosity ratio
        autoPtr<PatchFunction1<scalar>> viscosityRatio_;

        //- Name of the turbulent kinetic energy field
        word kName_;


public:

    //- Runtime type information
    TypeName("viscosityRatioInletOutletTDR");


    // Constructors

        //- Construct from patch and internal field
        viscosityRatioInletOutletTDRFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        viscosityRatioInletOutletTDRFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  viscosityRatioInletOutletTDRFvPatchScalarField onto a new patch
        viscosityRatioInletOutletTDRFvPatchScalarField
        (
            const viscosityRatioInletOutletTDRFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        viscosityRatioInletOutletTDRFvPatchScalarField
        (
            const viscosityRatioInletOutletTDRFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new viscosityRatioInletOutletTDRFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        viscosityRatioInletOutletTDRFvPatchScalarField
        (
            const viscosityRatioInletOutletTDRFvPatchScalarField&,
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
                new viscosityRatioInletOutletTDRFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access
        
            //- Return the name of the velocity field
            const word& kName() const
            {
                return kName_;
            }

            //- Return reference to the name of the velocity field
            //  to allow adjustment
            word& kName()
            {
                return kName_;
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

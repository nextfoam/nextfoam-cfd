/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
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
    Foam::subsonicOutletFvPatchField

Group
    grpOutletBoundaryConditions

Description
    Non-reflecting outlet boundary condition primarily for use with internal 
    flow.

Author
    JaeHeung Gill, NEXTFoam Co.,Ltd.

Usage
    \table
        Property     | Description               | Required | Default
        pExit        | outlet static pressure    | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    outlet
    {
        type            subsonicOutlet;
        pExit           100000;
        value           uniform 100000;
    }
    \endverbatim

Note

See also
    Foam::fixedValueFvPatchField


SourceFiles
    subsonicOutletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef subsonicOutletFvPatchField_H
#define subsonicOutletFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "turbulentFluidThermoModel.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class subsonicOutletFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class subsonicOutletFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        //- exit static pressure

        autoPtr<PatchFunction1<scalar>> pExit_;

        //- Current time index
        label curTimeIndex_;

    // Private member functions

        //- calculate and return current boundary field
        Field<Type> bcValue
        (
            const fvPatchField<scalar>& pp,
            const fvPatchField<vector>& pU,
            const fluidThermo& thermo
        ) const;


public:

    //- Runtime type information
    TypeName("subsonicOutlet");


    // Constructors

        //- Construct from patch and internal field
        subsonicOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        subsonicOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given subsonicOutletFvPatchField
        //  onto a new patch
        subsonicOutletFvPatchField
        (
            const subsonicOutletFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        subsonicOutletFvPatchField
        (
            const subsonicOutletFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new subsonicOutletFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        subsonicOutletFvPatchField
        (
            const subsonicOutletFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new subsonicOutletFvPatchField<Type>(*this, iF)
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
                const fvPatchField<Type>&,
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

#ifdef NoRepository
#   include "subsonicOutletFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

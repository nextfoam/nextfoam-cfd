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
    Foam::farfieldRiemannFvPatchField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition provides a non-reflecting free-stream condition 
    for compressible external flow based on Riemann invariants.

Author
    JaeHeung Gill, NEXTfoam Co.,Ltd.

Usage
    \table
        Property     | Description                      | Required | Default
        MInf         | free stream Mach number          | yes |
        pInf         | free stream static pressure      | yes |
        TInf         | free stream static temperature   | yes |
        flowDir      | flow direction                   | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        type            farfieldRiemann;
        MInf            0.729;
        pInf            108987.657;
        TInf            255.556;
        flowDir         (0.999187375576736 0.0403061842410544 0);
        value           uniform 108987.657;
    }
    \endverbatim

Note

See also
    Foam::fixedValueFvPatchField


SourceFiles
    farfieldRiemannFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef farfieldRiemannFvPatchField_H
#define farfieldRiemannFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class farfieldRiemannFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class farfieldRiemannFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        //- Mach number
        scalar MInf_;

        //- farfield static pressure
        scalar pInf_;

        //- farfield static temperature
        scalar TInf_;

        //- flow direction
        vector flowDir_;


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
    TypeName("farfieldRiemann");


    // Constructors

        //- Construct from patch and internal field
        farfieldRiemannFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        farfieldRiemannFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>& fld
        );

        //- Construct from patch, internal field and dictionary
        farfieldRiemannFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given farfieldRiemannFvPatchField
        //  onto a new patch
        farfieldRiemannFvPatchField
        (
            const farfieldRiemannFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        farfieldRiemannFvPatchField
        (
            const farfieldRiemannFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new farfieldRiemannFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        farfieldRiemannFvPatchField
        (
            const farfieldRiemannFvPatchField<Type>&,
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
                new farfieldRiemannFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return Mach number
            const scalar& MInf() const
            {
                return MInf_;
            }

            scalar& MInf()
            {
                return MInf_;
            }

            //- Return farfield static pressure
            const scalar& pInf() const
            {
                return pInf_;
            }

            scalar& pInf()
            {
                return pInf_;
            }

            //- Return farfield static temperature
            const scalar& TInf() const
            {
                return TInf_;
            }

            scalar& TInf()
            {
                return TInf_;
            }

            //- Return flow direction
            const vector& flowDir() const
            {
                return flowDir_;
            }

            vector& flowDir()
            {
                return flowDir_;
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
#   include "farfieldRiemannFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

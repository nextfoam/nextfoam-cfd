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
    Foam::omegaBlendedWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a wall function for the specific
    dissipation rate (i.e. \c omega) and the turbulent kinetic energy
    production contribution (i.e. \c G) for low- and high-Reynolds number
    applications. This boundary condition should be used in conjunction 
    with nutSpaldingWallFunction, then it automatically uses the wall 
    function coefficients specified in nutSpaldingWallFunction. This 
    boundary condition uses a blending approach that is tuned to have a 
    consistent wall shear stress and non-dimensional velocity profile(u+) 
    along y+ in the boundary layer on the flat plate.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries
        type            omegaBlendedWallFunction;

        // Optional entries
        beta1           0.075;
    }
    \endverbatim

    \table
      Property  | Description                     | Type   | Reqd | Deflt
      type      | Type name: omegaBlendedWallFunction    | word   | yes  | -
      beta1     | Model coefficient               | scalar | no   | 0.075
    \endtable

    The inherited entries are elaborated in:
      - \link wallFunctionCoefficients.H \endlink

Author
    JaeHeung Gill, NEXTfoam Co.,Ltd.

SourceFiles
    omegaBlendedWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef omegaBlendedWallFunctionFvPatchScalarField_H
#define omegaBlendedWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchField.H"
#include "nutSpaldingWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class turbulenceModel;

/*---------------------------------------------------------------------------*\
             Class omegaBlendedWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class omegaBlendedWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
protected:

    // Protected Data

        //- Tolerance used in weighted calculations
        static scalar tolerance_;

        //- Initialised flag
        bool initialised_;

        //- Master patch ID
        label master_;

        //- beta1 coefficient
        scalar beta1_;

        //- Local copy of turbulence G field
        scalarField G_;

        //- Local copy of turbulence omega field
        scalarField omega_;

        //- List of averaging corner weights
        List<List<scalar>> cornerWeights_;


    // Protected Member Functions

        //- Set the master patch - master is responsible for updating all
        //- wall function patches
        virtual void setMaster();

        //- Create the averaging weights for cells which are bounded by
        //- multiple wall function faces
        virtual void createAveragingWeights();

        //- Helper function to return const access to an nut patch
        const nutSpaldingWallFunctionFvPatchScalarField& nutPatch
        (
            const label patchi
        );

        //- Helper function to return non-const access to an omega patch
        virtual omegaBlendedWallFunctionFvPatchScalarField& omegaPatch
        (
            const label patchi
        );

        //- Main driver to calculate the turbulence fields
        virtual void calculateTurbulenceFields
        (
            const turbulenceModel& turbulence,
            scalarField& G0,
            scalarField& omega0
        );

        //- Calculate the omega and G
        virtual void calculate
        (
            const turbulenceModel& turbulence,
            const List<scalar>& cornerWeights,
            const fvPatch& patch,
            scalarField& G,
            scalarField& omega
        );

        //- Return non-const access to the master patch ID
        virtual label& master()
        {
            return master_;
        }

        //- Write local wall function variables
        void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("omegaBlendedWallFunction");


    // Constructors

        //- Construct from patch and internal field
        omegaBlendedWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        omegaBlendedWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- omegaBlendedWallFunctionFvPatchScalarField
        //- onto a new patch
        omegaBlendedWallFunctionFvPatchScalarField
        (
            const omegaBlendedWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        omegaBlendedWallFunctionFvPatchScalarField
        (
            const omegaBlendedWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new omegaBlendedWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        omegaBlendedWallFunctionFvPatchScalarField
        (
            const omegaBlendedWallFunctionFvPatchScalarField&,
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
                new omegaBlendedWallFunctionFvPatchScalarField(*this, iF)
            );
        }

    //- Destructor
    virtual ~omegaBlendedWallFunctionFvPatchScalarField() = default;


    // Member Functions

        // Access

            //- Return non-const access to the master's G field
            scalarField& G(bool init = false);

            //- Return non-const access to the master's omega field
            scalarField& omega(bool init = false);


        // Evaluation

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Update the coefficients associated with the patch field
            virtual void updateWeightedCoeffs(const scalarField& weights);

            //- Manipulate matrix
            virtual void manipulateMatrix(fvMatrix<scalar>& matrix);

            //- Manipulate matrix with given weights
            virtual void manipulateMatrix
            (
                fvMatrix<scalar>& matrix,
                const scalarField& weights
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

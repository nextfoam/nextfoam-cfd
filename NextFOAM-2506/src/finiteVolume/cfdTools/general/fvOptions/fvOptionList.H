/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2020,2023 PCOpt/NTUA
    Copyright (C) 2020,2023 FOSS GP
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
    Foam::fv::optionList

Description
    List of finite volume options

SourceFile
    optionList.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvOptionList_H
#define Foam_fvOptionList_H

#include "fvOption.H"
#include "PtrList.H"
#include "GeometricField.H"
#include "geometricOneField.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

namespace fv
{
    class optionList;
}

Ostream& operator<<(Ostream& os, const fv::optionList& options);

namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class optionList Declaration
\*---------------------------------------------------------------------------*/

class optionList
:
    public PtrList<fv::option>
{
protected:

    // Protected Data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Time index to check that all defined sources have been applied
        label checkTimeIndex_;

        //- Markers for porous media zone cells
        mutable autoPtr<volScalarField> porousZoneCellMarker_;

        //- Markers for porous media interface zone cells
        mutable autoPtr<labelList> porousInterfaceCells_;

        //- Markers for active porous media zone cells
        mutable autoPtr<labelList> activePorousCells_;

        //- Markers for porous media interface zone faces
        mutable autoPtr<surfaceScalarField> porousInterfaceFaces_;

        //- Markers for porous media interior zone faces
        mutable autoPtr<surfaceScalarField> porousInteriorFaces_;

        //- Markers for porous media zone faces
        mutable autoPtr<surfaceScalarField> porousFaces_;

        //- Linear interpolation weighting factors
        mutable autoPtr<surfaceScalarField> linearWeights_;

        //- Pressure gradient due to porous media resistance
        mutable autoPtr<volVectorField> porousGradpPtr_;

        //- Has active porous media zone?
        bool hasActivePorousZone_;


    // Protected Member Functions

        //- Return "options" sub-dictionary (if present) or return dict
        static const dictionary& optionsDict(const dictionary& dict);

        //- Read options dictionary
        bool readOptions(const dictionary& dict);

        //- Check that all sources have been applied
        void checkApplied() const;

        //- tangential momentum source for porousBafflePressure
        template<class RhoFieldType>
        void addTangentialUSource
        (
            fvMatrix<vector>& mtx,
            const RhoFieldType& rho,
            const GeometricField<vector, fvPatchField, volMesh>& field,
            const word& fieldName
        );

        template<class RhoFieldType>
        void addTangentialUSource
        (
            fvMatrix<scalar>& mtx,
            const RhoFieldType& rho,
            const GeometricField<scalar, fvPatchField, volMesh>& field,
            const word& fieldName
        )
        {}

        template<class RhoFieldType>
        void addTangentialUSource
        (
            fvMatrix<symmTensor>& mtx,
            const RhoFieldType& rho,
            const GeometricField<symmTensor, fvPatchField, volMesh>& field,
            const word& fieldName
        )
        {}

        //- Return source for equation with specified name and dimensions
        template<class Type>
        tmp<fvMatrix<Type>> source
        (
            GeometricField<Type, fvPatchField, volMesh>& field,
            const word& fieldName,
            const dimensionSet& ds
        );

        //- No copy construct
        optionList(const optionList&) = delete;

        //- No copy assignment
        void operator=(const optionList&) = delete;


public:

    //- Runtime type information
    TypeName("optionList");


    // Constructors

        //- Default construct from mesh
        explicit optionList(const fvMesh& mesh);

        //- Construct from mesh and dictionary
        optionList(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    virtual ~optionList();


    // Member Functions

        //- Reset the source list
        void reset(const dictionary& dict);

        //- Get porous media zone cell markers
        autoPtr<volScalarField> getPorousZoneCellMarker() const;

        //- Return porous media zone cell marker
        const volScalarField& porousZoneCellMarker() const;

        //- Get porous media interface zone cell markers
        autoPtr<labelList> getPorousInterfaceCells() const;

         //- Return reference to porous media interface zone cell markers
        const labelList& porousInterfaceCells() const;

        //- Get active porous media zone cell markers
        autoPtr<labelList> getActivePorousCells() const;

         //- Return reference to active porous media zone cell markers
        const labelList& activePorousCells() const;
       
        //- Get porous media zone interface markers
        tmp<surfaceScalarField> getPorousInterfaceFaces() const;

        //- Return reference to porous zone interface markers
        const surfaceScalarField& porousInterfaceFaces() const;

        //- Get porous media zone interior and non-coupled boundaries markers
        tmp<surfaceScalarField> getPorousInteriorFaces() const;

        //- Return reference to porous zone interior and non-coupled boundaries
        //markers
        const surfaceScalarField& porousInteriorFaces() const;

        //- Get porous media face markers
        tmp<surfaceScalarField> getPorousFaces() const;

        //- Return reference to porous media face markers
        const surfaceScalarField& porousFaces() const;

        //- Get linear interpolation weighting factors
        tmp<surfaceScalarField> getLinearWeights() const;

        //- Return linear interpolation weighting factors
        const surfaceScalarField& linearWeights() const;

        //- Correct interpolation weighting factors for porous media faces
        void correctPorousFaceWeights
        (
                surfaceScalarField& weights,
                bool flip = false
        ) const;

        //- Return whether there is something to apply to the field
        bool appliesToField(const word& fieldName) const;

        //- Find whether there is active porous media zone
        bool findActivePorous() const;

        //- Return whether there is active porous media zone
        bool hasActivePorousZone() const
        {
            return hasActivePorousZone_;
        }

        //- Update porous media solid properties regarding porosity
        virtual void adjustTransport
        (
            volSymmTensorField& kappaEff,
            volSymmTensorField& alphaEff,
            const volScalarField& Cp
        );


        // Sources

            //- Return source for equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                GeometricField<Type, fvPatchField, volMesh>& field
            );

            //- Return source for equation with specified name
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                GeometricField<Type, fvPatchField, volMesh>& field,
                const word& fieldName
            );

            //- Return source for equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field
            );

            //- Return source for equation with specified name
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field,
                const word& fieldName
            );

            //- Return source for equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field
            );

            //- Return source for equation with specified name
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field,
                const word& fieldName
            );

            //- Return source for equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const volScalarField& alpha,
                const geometricOneField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field
            );

            //- Return source for equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const geometricOneField& alpha,
                const volScalarField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field
            );

            //- Return source for equation
            template<class Type>
            tmp<fvMatrix<Type>> operator()
            (
                const geometricOneField& alpha,
                const geometricOneField& rho,
                GeometricField<Type, fvPatchField, volMesh>& field
            );

            //- Return source for equation with second time derivative
            template<class Type>
            tmp<fvMatrix<Type>> d2dt2
            (
                GeometricField<Type, fvPatchField, volMesh>& field
            );

            //- Return source for equation with second time derivative
            template<class Type>
            tmp<fvMatrix<Type>> d2dt2
            (
                GeometricField<Type, fvPatchField, volMesh>& field,
                const word& fieldName
            );


        // Constraints

            //- Apply constraints to equation
            template<class Type>
            void constrain(fvMatrix<Type>& eqn);


        // Correction

            //- Apply correction to field
            template
            <
                class Type, template<class> class PatchField, class GeoMesh
            >
            void correct(GeometricField<Type, PatchField, GeoMesh>& field);


        //- Post process sensitivity field related to the fvOption
        template<class Type>
        void postProcessSens
        (
            Field<Type>& sensField,
            const word& fieldName = word::null,
            const word& designVariablesName = word::null
        );

        //- Post process auxiliary sensitivity field related to the fvOption
        template<class Type, template<class> class PatchField, class GeoMesh>
        void postProcessAuxSens
        (
            const GeometricField<Type, PatchField, GeoMesh>& primal,
            const GeometricField<Type, PatchField, GeoMesh>& adjoint,
            scalarField& sensField,
            const word& fieldName = word::null
        );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);

            //- Write data to Ostream
            virtual bool writeData(Ostream& os) const;

            //- Ostream operator
            friend Ostream& operator<<
            (
                Ostream& os,
                const optionList& options
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvOptionListTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

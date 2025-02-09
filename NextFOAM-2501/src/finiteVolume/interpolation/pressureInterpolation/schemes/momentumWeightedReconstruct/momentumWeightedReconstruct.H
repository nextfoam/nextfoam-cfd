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
    Foam::momentumWeightedReconstruct

Group
    grpFvPressureInterpolationSchemes

Description
    Pressure interpolation scheme class which uses discretized momentum 
    equation's diagonal coefficients as weighting factors.

Author
    JaeHeung Gill, NEXTfoam Co.,Ltd.

SourceFiles
    momentumWeightedReconstruct.C

\*---------------------------------------------------------------------------*/

#ifndef momentumWeightedReconstruct_H
#define momentumWeightedReconstruct_H

#include "pressureInterpolationScheme.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class momentumWeightedReconstruct Declaration
\*---------------------------------------------------------------------------*/

class momentumWeightedReconstruct
:
    public pressureInterpolationScheme
{
    // Private Data

        mutable bool correct_;

        mutable surfaceScalarField weights_;


    // Private Member Functions

        //- No copy assignment
        void operator=(const momentumWeightedReconstruct&) = delete;


public:

    //- Runtime type information
    TypeName("momentumWeightedReconstruct");


    // Constructors

        //- Construct from mesh
        momentumWeightedReconstruct(const fvMesh& mesh)
        :
            pressureInterpolationScheme(mesh),
            correct_(true),
            weights_(calcWeights())
        {}

        //- Construct from Istream
        momentumWeightedReconstruct(const fvMesh& mesh, Istream&)
        :
            pressureInterpolationScheme(mesh),
            correct_(true),
            weights_(calcWeights())
        {}

        //- Construct from faceFlux and Istream
        momentumWeightedReconstruct
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream&
        )
        :
            pressureInterpolationScheme(mesh),
            correct_(true),
            weights_(calcWeights())
        {}


    // Member Functions

        //- Calculate the interpolation weighting factors
        tmp<surfaceScalarField> calcWeights() const
        {
            dictionary controlledVariables;
            
            if (mesh_.foundObject<IOdictionary>("controlledVariables"))
            {
                controlledVariables = 
                    mesh_.lookupObject<IOdictionary>("controlledVariables");
            }

            if (controlledVariables.found("boundedVariables"))
            {
                wordHashSet boundedVariables
                (
                    controlledVariables.get<wordHashSet>("boundedVariables")
                );

                if 
                (
                    boundedVariables.found("p_rgh")
                 || boundedVariables.found("p")
                )
                {
                    correct_ = false;
                }
            }

            bool calcWeight(true);

            IOobject phiHeader
            (
                "phi",
                Foam::Time::timeName(mesh_.time().startTime().value()),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            );

            if 
            ( 
                !phiHeader.typeHeaderOk<surfaceScalarField>(true)
             && mesh_.time().timeIndex() == (mesh_.time().startTimeIndex() + 1)
            )
            {
                calcWeight = false;
                correct_ = false;
            }

            tmp<surfaceScalarField> tw
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        "momentumWeightedReconstructWeights",
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    mesh_.surfaceInterpolation::weights()
                )
            );

            if (calcWeight)
            {
                tmp<volScalarField> tAU
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "UAp",
                            mesh_.time().timeName(),
                            mesh_
                        ),
                        mesh_.lookupObject<volScalarField>("UAp")
                    )
                );

                const scalarField& AU = tAU.ref().primitiveField();

                scalarField& w = tw.ref().primitiveFieldRef();

                const labelUList& own = mesh_.owner();
                const labelUList& nei = mesh_.neighbour();

                if (correct_)
                {
                    forAll(own, facei)
                    {
                        scalar AOwn = AU[own[facei]];
                        scalar ANei = AU[nei[facei]];

                        w[facei] = AOwn/(AOwn + ANei);
                    }
                }
                else
                {
                    forAll(own, facei)
                    {
                        scalar AOwn = AU[own[facei]];
                        scalar ANei = AU[nei[facei]];

                        w[facei] = ANei/(AOwn + ANei);
                    }
                }

                surfaceScalarField::Boundary& bw = tw.ref().boundaryFieldRef();

                forAll(bw, patchi)
                {
                    const fvPatchScalarField& pAU = 
                        tAU().boundaryField()[patchi];

                    if (pAU.coupled())
                    {
                        tmp<scalarField> tAUP(pAU.patchInternalField());
                        const scalarField& AUP = tAUP();

                        tmp<scalarField> tAUN(pAU.patchNeighbourField());
                        const scalarField& AUN = tAUN();

                        if (correct_)
                        {
                            forAll(pAU, facei)
                            {
                                scalar AOwn = AUP[facei];
                                scalar ANei = AUN[facei];

                                bw[patchi][facei] = AOwn/(AOwn + ANei);
                            }
                        }
                        else
                        {
                            forAll(pAU, facei)
                            {
                                scalar AOwn = AUP[facei];
                                scalar ANei = AUN[facei];

                                bw[patchi][facei] = ANei/(AOwn + ANei);
                            }
                        }
                    }
                }
            }

            return tw;
        }

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights(const volScalarField&) const
        {
            if (hasActivePorousZone_)
            {
                fvOptions_.correctPorousFaceWeights(weights_, true);
            }

            return weights_;
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return correct_;
        }

        //- Return the explicit correction to the face-interpolate
        virtual tmp<surfaceScalarField> correction
        (
            const volScalarField&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

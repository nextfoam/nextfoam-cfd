/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::RASModels::realizableKEtwoLayer

Group
    grpRASTurbulence

Description
    Realizable k-epsilon turbulence model with two-layer approach for 
    incompressible and compressible flows.

    References:
    \verbatim
        Shih, T. H., Liou, W. W., Shabbir, A., Yang, Z., & Zhu, J. (1994).
        A new k-epsilon eddy viscosity model for high Reynolds number
        turbulent flows: Model development and validation.
        NASA STI/Recon Technical Report N, 95, 11442.

        Shih, T. H., Liou, W. W., Shabbir, A., Yang, Z., & Zhu, J. (1995).
        A New k-epsilon Eddy Viscosity Model for High Reynolds Number
        Turbulent Flows.
        Computers and Fluids, 24(3), 227-238.

        Chen, H. C., Patel, V. C. (1988).
        Near-Wall Turbulence Models for Complex Flows Including Separation. 
        AIAA Journal, 26(6), 641–648.

        Jongen, T. (1992).
        Simulation and Modeling of Turbulent Incompressible Flows. 
        PhD thesis. EPF Lausanne, Lausanne, Switzerland.

        Wolfshtein, M. (1969).
        The Velocity and Temperature Distribution of One-Dimensional Flow with 
        Turbulence Augmentation and Pressure Gradient. 
        Int. J. Heat Mass Transfer. 12. 301–318.
    \endverbatim

    The default model coefficients are
    \verbatim
        realizableKEtwoLayerCoeffs
        {
            A0          4.04;
            C2          1.9;
            sigmak      1.0;
            sigmaEps    1.2;
            Cmu         0.09;
            kappa       0.41;
            ReyStar     60;
            deltaRey    10;
            Anu         70;
        }
    \endverbatim

Author
    JaeHeung Gill, NEXTfoam Co.,Ltd.

SourceFiles
    realizableKEtwoLayer.C

\*---------------------------------------------------------------------------*/

#ifndef realizableKEtwoLayer_H
#define realizableKEtwoLayer_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                        Class realizableKE Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class realizableKEtwoLayer
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{

protected:

    // Protected data

        // Model coefficients

            // for realizableKE
            dimensionedScalar A0_;
            dimensionedScalar C2_;
            dimensionedScalar sigmak_;
            dimensionedScalar sigmaEps_;

            // for realizableKEtwoLayer
            dimensionedScalar Cmu_;
            dimensionedScalar Cmu75_;
            dimensionedScalar kappa_;
            dimensionedScalar ReyStar_;
            dimensionedScalar deltaRey_;
            dimensionedScalar Anu_;

            // for turbulent source due to buoyancy
            dimensionedScalar Prt_;


        // Fields

            volScalarField k_;
            volScalarField epsilon_;

            volScalarField rCmu_;
            volScalarField::Internal lambda_;
            volScalarField::Internal Rey_;

            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField& y_;


   // Protected Member Functions

        tmp<volScalarField> rCmu
        (
            const volTensorField& gradU,
            const volScalarField& S2,
            const volScalarField& magS
        );

        virtual void correctNut
        (
            const volTensorField& gradU,
            const volScalarField& S2,
            const volScalarField& magS
        );

        void blendAndRelax
        (
            fvScalarMatrix& eqn, 
            const volScalarField::Internal& lambda,
            const volScalarField::Internal& es
        );

        virtual void correctNut();
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> epsilonSource() const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("realizableKEtwoLayer");

    // Constructors

        //- Construct from components
        realizableKEtwoLayer
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~realizableKEtwoLayer() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    (this->nut_/sigmak_ + this->nu())
                )
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DepsilonEff",
                    (this->nut_/sigmaEps_ + this->nu())
                )
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "realizableKEtwoLayer.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

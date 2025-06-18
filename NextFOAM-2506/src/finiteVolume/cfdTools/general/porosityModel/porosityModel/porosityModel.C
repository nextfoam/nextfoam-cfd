/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "porosityModel.H"
#include "zeroGradientFvPatchFields.H"
#include "solutionControl.H"
#include "multiRegionSolutionControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(porosityModel, 0);
    defineRunTimeSelectionTable(porosityModel, mesh);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::porosityModel::adjustNegativeResistance(dimensionedVector& resist)
{
    scalar maxCmpt = cmptMax(resist.value());

    if (maxCmpt < 0)
    {
        FatalErrorInFunction
            << "Cannot have all resistances set to negative, resistance = "
            << resist
            << exit(FatalError);
    }
    else
    {
        vector& val = resist.value();
        for (label cmpt = 0; cmpt < vector::nComponents; cmpt++)
        {
            if (val[cmpt] < 0)
            {
                val[cmpt] *= -maxCmpt;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModel::porosityModel
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const wordRe& cellZoneName
)
:
    regIOobject
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    name_(name),
    mesh_(mesh),
    dict_(dict),
    coeffs_(dict.optionalSubDict(modelType + "Coeffs")),
    active_(true),
    zoneName_(cellZoneName),
    cellZoneIDs_(),
    csysPtr_
    (
        coordinateSystem::New(mesh, coeffs_, coordinateSystem::typeName)
    ), 
    porosity_
    (
        dimensionedScalar::getOrAddToDict("porosity", coeffs_, dimless, 1.0)
    ),
    porosityBf_
    (
        mesh_.boundary(),
        volScalarField::Internal::null(),
        fvPatchFieldBase::calculatedType()
    ),
    steady_(true),
    solveEnergy_(false),
    dualCell_
    (
        coeffs_.getOrDefault<bool>("dualCell", false)
    )
{
    if (zoneName_.empty())
    {
        dict.readIfPresent("active", active_);
        dict_.readEntry("cellZone", zoneName_);
    }

    cellZoneIDs_ = mesh_.cellZones().indices(zoneName_);

    Info<< "    creating porous zone: " << zoneName_ << endl;

    if (returnReduceAnd(cellZoneIDs_.empty()) && Pstream::master())
    {
        FatalErrorInFunction
            << "Cannot find porous cellZone " << zoneName_ << endl
            << "Valid zones : "
            << flatOutput(mesh_.cellZones().names()) << nl
            << "Valid groups: "
            << flatOutput(mesh_.cellZones().groupNames()) << nl
            << exit(FatalError);
    }

    Info<< incrIndent << indent << csys() << decrIndent << endl;

    const pointField& points = mesh_.points();
    const cellList& cells = mesh_.cells();
    const faceList& faces = mesh_.faces();

    for (const label zonei : cellZoneIDs_)
    {
        const cellZone& cZone = mesh_.cellZones()[zonei];

        boundBox bb;

        for (const label celli : cZone)
        {
            const cell& c = cells[celli];
            const pointField cellPoints(c.points(faces, points));

            for (const point& pt : cellPoints)
            {
                bb.add(csys().localPosition(pt));
            }
        }

        bb.reduce();

        Info<< "    local bounds: " << bb.span() << nl << endl;
    }

    checkSolverState();

    updatePorosity();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModel::checkSolverState()
{
    dualCell_ = coeffs_.getOrDefault<bool>("dualCell", false);

    if (mesh_.foundObject<surfaceScalarField>("faceRho"))
    {
        if (mesh_.foundObject<solutionControl>("solutionControl"))
        {
            const solutionControl& control
            (
                mesh_.lookupObject<solutionControl>("solutionControl")
            );

            solveEnergy_ = control.solveEnergy();

            if (control.algorithmName() == "SIMPLE")
            {
                steady_ = true;
            }
            else
            {
                steady_ = false;
            }
        }
        else if 
        (
            mesh_.time().foundObject<multiRegionSolutionControl>("SIMPLE")
        )
        {
            steady_ = true;

            const dictionary& control
            (
                mesh_.solutionDict().subDict("SIMPLE")
            );

            solveEnergy_ = control.getOrDefault<bool>("solveEnergy", true);
        }
        else if 
        (
            mesh_.time().foundObject<multiRegionSolutionControl>("PIMPLE")
        )
        {
            steady_ = false;

            const dictionary& control
            (
                mesh_.solutionDict().subDict("PIMPLE")
            );

            solveEnergy_ = control.getOrDefault<bool>("solveEnergy", true);
        }
    }
}


void Foam::porosityModel::updatePorosity()
{
    if (solveEnergy_ && porosity_.value() < 1.0)
    {
        auto tPorosityv
        (
            tmp<volScalarField>::New
            (
                IOobject
                (
                    "porosityv_" + zoneName_,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("unitPorosityv", dimless, scalar(1)),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        volScalarField& porosityv = tPorosityv.ref();

        for (const label zonei : cellZoneIDs_)
        {
            const cellZone& cZone = mesh_.cellZones()[zonei];

            for (const label celli : cZone)
            {
                porosityv[celli] = porosity_.value();   
            }
        }

        porosityv.correctBoundaryConditions();

        porosityBf_ == porosityv.boundaryField();

        tPorosityv.clear();

        if (dualCell_)
        {
            const word nbrRegionName = coeffs_.get<word>("solidRegionName");

            const fvMesh& nbrRegion = 
                mesh_.time().lookupObject<fvMesh>(nbrRegionName);

            volScalarField& betavSolid = const_cast<volScalarField&>
            (
                nbrRegion.lookupObject<volScalarField>("betavSolid")
            );

            betavSolid = 1.0 - porosity_;

            betavSolid.correctBoundaryConditions();

            coeffs_.add
            (
                "solidProperties", 
                dictionary
                (
                    nbrRegion.lookupObject<IOdictionary>
                    (
                        "thermophysicalProperties"
                    )
                ),
                true
            );

            //- Verify that the Dual Cell Heat Exchange model dictionaries 
            //  are set
            word hxDictName("HX_" + mesh_.name() + "-to-" + nbrRegionName);
            word nbrHxDictName("HX_" + nbrRegionName + "-to-" + mesh_.name());

            if 
            (
                !mesh_.lookupObject<IOdictionary>
                (
                    "fvOptions"
                ).toc().found(hxDictName)
             || !nbrRegion.lookupObject<IOdictionary>
                (
                    "fvOptions"
                ).toc().found(nbrHxDictName)
            )
            {
                FatalErrorInFunction
                    << "The inter-region heat transfer source must be " 
                    << "specified in the 'fvOptions' file as "  
                    << "'" << hxDictName << "' for region " << mesh_.name() 
                    << " and "
                    << "'" << nbrHxDictName << "' for region " 
                    << nbrRegionName 
                    << abort(FatalError);
            }
        }
    }
}


void Foam::porosityModel::adjustTransport
(
    volSymmTensorField& kappaEff,
    volSymmTensorField& alphaEff,
    const volScalarField& Cp
)
{
    if (porosity_.value() < 1.0)
    {
        scalar gamma(porosity_.value());
        scalar beta(1.0 - gamma);

        const dictionary& solidProperties(coeffs_.subDict("solidProperties"));

        dimensioned<symmTensor> kps("kappas", kappaEff.dimensions(), Zero);

        const word kappaType
        (
            solidProperties.subDict("thermoType").getWord("transport")
        );

        const dictionary& transport
        (
            solidProperties.subDict("mixture").subDict("transport")
        );

        if (kappaType == "constIso")
        {
            kps.value() = symmTensor::I*transport.get<scalar>("kappa");
        }
        else if (kappaType == "constAnIso")
        {
            alphaEff.rename("Anialpha");
            kappaEff.rename("Anikappa");

            mesh_.solution().enableCache("grad(T)");

            autoPtr<coordinateSystem> coordSys
            (
                coordinateSystem::New
                (
                    mesh_,
                    solidProperties,
                    coordinateSystem::typeName
                )
            );

            kps.value() = coordSys().transformPrincipal
            (
                transport.get<vector>("kappa")
            );
        }

        for (const label zonei : cellZoneIDs_)
        {
            const cellZone& cZone = mesh_.cellZones()[zonei];

            for (const label celli : cZone)
            {
                symmTensor& kappa = kappaEff.primitiveFieldRef()[celli];
                symmTensor& alpha = alphaEff.primitiveFieldRef()[celli];

                kappa *= gamma;
                alpha *= gamma;

                if (!dualCell_)
                {
                    kappa += beta*kps.value();
                    alpha += beta*kps.value()/Cp[celli];
                }
            }

            auto& kappaBf = kappaEff.boundaryFieldRef();
            auto& alphaBf = alphaEff.boundaryFieldRef();
            const auto& cpBf = Cp.boundaryField();

            forAll(kappaBf, patchi)
            {
                kappaBf[patchi] *= porosityBf_[patchi];
                alphaBf[patchi] *= porosityBf_[patchi];

                if (!dualCell_)
                {
                    kappaBf[patchi] += 
                        (1.0 - porosityBf_[patchi])*kps.value();
                    alphaBf[patchi] += 
                        (1.0 - porosityBf_[patchi])*kps.value()/cpBf[patchi];
                }
            }
        }
    }
}


void Foam::porosityModel::transformModelData()
{
    if (!mesh_.upToDatePoints(*this))
    {
        calcTransformModelData();

        // set model up-to-date wrt points
        mesh_.setUpToDatePoints(*this);
    }
}


Foam::tmp<Foam::vectorField> Foam::porosityModel::porosityModel::force
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu
)
{
    transformModelData();

    auto tforce = tmp<vectorField>::New(U.size(), Zero);

    if (!cellZoneIDs_.empty())
    {
        this->calcForce(U, rho, mu, tforce.ref());
    }

    return tforce;
}


void Foam::porosityModel::addResistance(fvVectorMatrix& UEqn)
{
    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    this->correct(UEqn);
}


void Foam::porosityModel::addResistance
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
)
{
    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    this->correct(UEqn, rho, mu);
}


void Foam::porosityModel::addResistance
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU,
    bool correctAUprocBC
)
{
    if (cellZoneIDs_.empty())
    {
        return;
    }

    transformModelData();
    this->correct(UEqn, AU);

    if (correctAUprocBC)
    {
        // Correct the boundary conditions of the tensorial diagonal to ensure
        // processor boundaries are correctly handled when AU^-1 is interpolated
        // for the pressure equation.
        AU.correctBoundaryConditions();
    }
}


bool Foam::porosityModel::writeData(Ostream& os) const
{
    return true;
}


bool Foam::porosityModel::read(const dictionary& dict)
{
    dict.readIfPresent("active", active_);

    coeffs_ = dict.optionalSubDict(type() + "Coeffs");

    dict.readEntry("cellZone", zoneName_);
    cellZoneIDs_ = mesh_.cellZones().indices(zoneName_);

    porosity_ = dimensionedScalar::getOrAddToDict
    (
        "porosity", 
        coeffs_, 
        dimless, 
        1.0
    );

    checkSolverState();

    updatePorosity();

    return true;
}


// ************************************************************************* //

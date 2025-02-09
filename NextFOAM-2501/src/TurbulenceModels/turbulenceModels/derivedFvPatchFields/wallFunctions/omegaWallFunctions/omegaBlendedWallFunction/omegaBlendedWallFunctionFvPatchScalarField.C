/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016, 2019 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "omegaBlendedWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::omegaBlendedWallFunctionFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::omegaBlendedWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const auto& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<omegaBlendedWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            omegaBlendedWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            opf.master() = master;
        }
    }
}


void Foam::omegaBlendedWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const auto& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    const fvMesh& mesh = omega.mesh();

    if (initialised_ && !mesh.changing())
    {
        return;
    }

    volScalarField weights
    (
        IOobject
        (
            "weights",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    );

    DynamicList<label> omegaPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<omegaBlendedWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            omegaPatches.append(patchi);

            const labelUList& faceCells = bf[patchi].patch().faceCells();
            for (const auto& celli : faceCells)
            {
                ++weights[celli];
            }
        }
    }

    cornerWeights_.setSize(bf.size());
    for (const auto& patchi : omegaPatches)
    {
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    G_.setSize(internalField().size(), 0.0);
    omega_.setSize(internalField().size(), 0.0);

    initialised_ = true;
}


const Foam::nutSpaldingWallFunctionFvPatchScalarField&
Foam::omegaBlendedWallFunctionFvPatchScalarField::nutPatch
(
    const label patchi
)
{
    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const volScalarField::Boundary& bf = turbModel.nut()().boundaryField();

    return refCast<const nutSpaldingWallFunctionFvPatchScalarField>(bf[patchi]);
}


Foam::omegaBlendedWallFunctionFvPatchScalarField&
Foam::omegaBlendedWallFunctionFvPatchScalarField::omegaPatch
(
    const label patchi
)
{
    const auto& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    const auto& opf =
        refCast<const omegaBlendedWallFunctionFvPatchScalarField>(bf[patchi]);

    return const_cast<omegaBlendedWallFunctionFvPatchScalarField&>(opf);
}


void Foam::omegaBlendedWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbModel,
    scalarField& G0,
    scalarField& omega0
)
{
    // accumulate all of the G and omega contributions
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            omegaBlendedWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            opf.calculate(turbModel, w, opf.patch(), G0, omega0);
        }
    }

    // apply zero-gradient condition for omega
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            omegaBlendedWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            opf == scalarField(omega0, opf.patch().faceCells());
        }
    }
}


void Foam::omegaBlendedWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& omega0
)
{
    const label patchi = patch.index();

    const scalar Cmu25 = pow025(nutPatch(patchi).wallCoeffs().Cmu());
    const scalar kappa = nutPatch(patchi).wallCoeffs().kappa();
    const scalar E = nutPatch(patchi).wallCoeffs().E();
    const scalar yPlusLam = nutPatch(patchi).wallCoeffs().yPlusLam();

    const scalarField& y = turbModel.y()[patchi];

    const labelUList& faceCells = patch.faceCells();

    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const scalarField& nutw = tnutw();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradUw(mag(Uw.snGrad()));
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    const scalarField uTau(sqrt((nutw + nuw)*magGradUw));

    forAll(faceCells, facei)
    {
        const label& celli = faceCells[facei];
        const scalar& w = cornerWeights[facei];

        const scalar yPlus = uTau[facei]*y[facei]/nuw[facei];

        scalar omegaVis = w*6*nuw[facei]/(beta1_*sqr(y[facei]));
        const scalar omegaLog = w*sqrt(k[celli])/(Cmu25*kappa*y[facei]);

        scalar yPlusSub = 1.0;
        scalar scaleOmegaVis(0.5);

        if (yPlus < yPlusSub)
        {
            omega0[celli] += w*omegaVis*scaleOmegaVis;

        }
        else
        {
            omegaVis *= 
                ((1.0 - scaleOmegaVis)*tanh(yPlus - 1) + scaleOmegaVis);
            scalar op = tanh(pow((yPlus - yPlusSub)/10, 0.81));
            scalar bf = 0.5;

            omega0[celli] +=
                pow(pow(omegaVis*(1 - op), bf) + pow(omegaLog*op, bf), 1/bf);
        }

        scalar kUp = min(kappa*magUp[facei]/max(uTau[facei], VSMALL), 40);
        scalar dypdup = 1.0 + kappa*(exp(kUp) - 1.0 - kUp - 0.5*sqr(kUp))/E;
        scalar dudy = sqr(uTau[facei])/mag(dypdup)/nuw[facei];

        G0[celli] += w*nutw[facei]*magGradUw[facei]*dudy;
    }
}


void Foam::omegaBlendedWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    os.writeEntryIfDifferent<scalar>("beta1", 0.075, beta1_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::omegaBlendedWallFunctionFvPatchScalarField::omegaBlendedWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    initialised_(false),
    master_(-1),
    beta1_(0.075),
    G_(),
    omega_(),
    cornerWeights_()
{}


Foam::omegaBlendedWallFunctionFvPatchScalarField::omegaBlendedWallFunctionFvPatchScalarField
(
    const omegaBlendedWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    initialised_(false),
    master_(-1),
    beta1_(ptf.beta1_),
    G_(),
    omega_(),
    cornerWeights_()
{}


Foam::omegaBlendedWallFunctionFvPatchScalarField::omegaBlendedWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    initialised_(false),
    master_(-1),
    beta1_(dict.getOrDefault<scalar>("beta1", 0.075)),
    G_(),
    omega_(),
    cornerWeights_()
{
    // Apply zero-gradient condition on start-up
    this->extrapolateInternal();
}


Foam::omegaBlendedWallFunctionFvPatchScalarField::omegaBlendedWallFunctionFvPatchScalarField
(
    const omegaBlendedWallFunctionFvPatchScalarField& owfpsf
)
:
    fixedValueFvPatchField<scalar>(owfpsf),
    initialised_(false),
    master_(-1),
    beta1_(owfpsf.beta1_),
    G_(),
    omega_(),
    cornerWeights_()
{}


Foam::omegaBlendedWallFunctionFvPatchScalarField::omegaBlendedWallFunctionFvPatchScalarField
(
    const omegaBlendedWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(owfpsf, iF),
    initialised_(false),
    master_(-1),
    beta1_(owfpsf.beta1_),
    G_(),
    omega_(),
    cornerWeights_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField& Foam::omegaBlendedWallFunctionFvPatchScalarField::G
(
    bool init
)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            G_ = 0.0;
        }

        return G_;
    }

    return omegaPatch(master_).G();
}


Foam::scalarField& Foam::omegaBlendedWallFunctionFvPatchScalarField::omega
(
    bool init
)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            omega_ = 0.0;
        }

        return omega_;
    }

    return omegaPatch(master_).omega(init);
}


void Foam::omegaBlendedWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), omega(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& omega0 = this->omega();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G = db().lookupObjectRef<FieldType>(turbModel.GName());

    FieldType& omega = const_cast<FieldType&>(internalField());

    forAll(*this, facei)
    {
        const label celli = patch().faceCells()[facei];

        G[celli] = G0[celli];
        omega[celli] = omega0[celli];
    }

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::omegaBlendedWallFunctionFvPatchScalarField::updateWeightedCoeffs
(
    const scalarField& weights
)
{
    if (updated())
    {
        return;
    }

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), omega(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& omega0 = this->omega();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G = db().lookupObjectRef<FieldType>(turbModel.GName());

    FieldType& omega = const_cast<FieldType&>(internalField());

    scalarField& omegaf = *this;

    // only set the values if the weights are > tolerance
    forAll(weights, facei)
    {
        const scalar w = weights[facei];

        if (w > tolerance_)
        {
            const label celli = patch().faceCells()[facei];

            G[celli] = (1.0 - w)*G[celli] + w*G0[celli];
            omega[celli] = (1.0 - w)*omega[celli] + w*omega0[celli];
            omegaf[facei] = omega[celli];
        }
    }

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::omegaBlendedWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    matrix.setValues(patch().faceCells(), patchInternalField());

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void Foam::omegaBlendedWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const Field<scalar>& weights
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    DynamicList<label> constraintCells(weights.size());
    DynamicList<scalar> constraintValues(weights.size());
    const labelUList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& fld = internalField();

    forAll(weights, facei)
    {
        // only set the values if the weights are > tolerance
        if (tolerance_ < weights[facei])
        {
            const label celli = faceCells[facei];

            constraintCells.append(celli);
            constraintValues.append(fld[celli]);
        }
    }

    if (debug)
    {
        Pout<< "Patch: " << patch().name()
            << ": number of constrained cells = " << constraintCells.size()
            << " out of " << patch().size()
            << endl;
    }

    matrix.setValues(constraintCells, constraintValues);

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void Foam::omegaBlendedWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        omegaBlendedWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //

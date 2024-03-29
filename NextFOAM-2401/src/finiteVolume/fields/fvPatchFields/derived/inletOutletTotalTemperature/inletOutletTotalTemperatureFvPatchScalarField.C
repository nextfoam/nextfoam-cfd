/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "inletOutletTotalTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inletOutletTotalTemperatureFvPatchScalarField::
inletOutletTotalTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    UName_("U"),
    psiName_("thermo::psi"),
    gamma_(1.0),
    T0_(nullptr)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


Foam::inletOutletTotalTemperatureFvPatchScalarField::
inletOutletTotalTemperatureFvPatchScalarField
(
    const inletOutletTotalTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_),
    T0_(ptf.T0_.clone(p.patch()))
{}


Foam::inletOutletTotalTemperatureFvPatchScalarField::
inletOutletTotalTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    UName_(dict.getOrDefault<word>("U", "U")),
    psiName_(dict.getOrDefault<word>("psi", "thermo:psi")),
    gamma_(dict.get<scalar>("gamma")),
    T0_(PatchFunction1<scalar>::New(p.patch(), "T0", dict))
{
    fvPatchFieldBase::readDict(dict);

    this->phiName_ = dict.getOrDefault<word>("phi", "phi");

    this->refValue() = Zero;

    if (!this->readValueEntry(dict))
    {
        fvPatchField<scalar>::operator=
        (
            T0_->value(this->db().time().timeOutputValue())
        );
    }

    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


Foam::inletOutletTotalTemperatureFvPatchScalarField::
inletOutletTotalTemperatureFvPatchScalarField
(
    const inletOutletTotalTemperatureFvPatchScalarField& tppsf
)
:
    inletOutletFvPatchScalarField(tppsf),
    UName_(tppsf.UName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    T0_(tppsf.T0_.clone(this->patch().patch()))
{}


Foam::inletOutletTotalTemperatureFvPatchScalarField::
inletOutletTotalTemperatureFvPatchScalarField
(
    const inletOutletTotalTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    T0_(tppsf.T0_.clone(this->patch().patch()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::inletOutletTotalTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    inletOutletFvPatchScalarField::autoMap(m);

    if (T0_)
    {
        T0_->autoMap(m);
    }
}


void Foam::inletOutletTotalTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    inletOutletFvPatchScalarField::rmap(ptf, addr);

    const inletOutletTotalTemperatureFvPatchScalarField& tiptf =
        refCast<const inletOutletTotalTemperatureFvPatchScalarField>(ptf);

    if (T0_)
    {
        T0_->rmap(tiptf.T0_(), addr);
    }
}


void Foam::inletOutletTotalTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const auto& Up =
        patch().lookupPatchField<volVectorField>(UName_);

    const auto& phip =
        patch().lookupPatchField<surfaceScalarField>(this->phiName_);

    const auto& psip =
        patch().lookupPatchField<volScalarField>(psiName_);

    scalar gM1ByG = (gamma_ - 1.0)/gamma_;

    this->refValue() =
        T0_->value(this->db().time().timeOutputValue())
        /(1.0 + 0.5*psip*gM1ByG*(neg(phip))*magSqr(Up));
    this->valueFraction() = neg(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}


void Foam::inletOutletTotalTemperatureFvPatchScalarField::write(Ostream& os)
const
{
    fvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("phi", "phi", this->phiName_);
    os.writeEntryIfDifferent<word>("psi", "thermo::psi", psiName_);
    os.writeEntry("gamma", gamma_);

    if (T0_)
    {
        T0_->writeData(os);
    }

    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        inletOutletTotalTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //

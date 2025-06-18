/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "subsonicOutletFvPatchField.H"
#include "thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Field<Type> subsonicOutletFvPatchField<Type>::bcValue
(
            const fvPatchField<scalar>& pp,
            const fvPatchField<vector>& pU,
            const fluidThermo& thermo
) const
{
    const scalar t = this->db().time().timeOutputValue();

    const label patchi = this->patch().index();

    const tmp<volScalarField> tgamma = thermo.gamma();
    const scalarField gammai
    (
        tgamma().boundaryField()[patchi].patchInternalField()
    );

    const scalar molWeight
    (
        readScalar
        (
            thermo.subDict("mixture").subDict("specie").lookup("molWeight")
        )
    );

    const scalar R(RR/molWeight);

    const vectorField n(this->patch().nf());

    const scalarField pi(pp.patchInternalField());
    const vectorField Ui(pU.patchInternalField());
    const scalarField Ti
    (
        thermo.T().boundaryField()[this->patch().index()].patchInternalField()
    );

    scalarField ci(sqrt(gammai*R*Ti));
    scalarField Mi((Ui & n)/ci);

    Field<Type> bcValue(this->patch().size(), pTraits<Type>::zero);

    if(this->internalField().name() == "p")
    {
        scalarField pb(pExit_->value(t));

        forAll(Mi, facei)
        {
            if(Mi[facei] >= 1.0)
            {
                pb[facei] = pi[facei];
            }
        }

        bcValue = pb*pTraits<Type>::one;
    }
    else if(this->internalField().name() == "T")
    {
        bcValue = Ti*pTraits<Type>::one;
    }

    return bcValue;
}


template<>
Field<vector> subsonicOutletFvPatchField<vector>::bcValue
(
            const fvPatchField<scalar>& pp,
            const fvPatchField<vector>& pU,
            const fluidThermo& thermo
) const
{
    const vectorField Ui(pU.patchInternalField());

    Field<vector> bcValue(this->patch().size(), pTraits<vector>::zero);

    if(this->internalField().name() == "U")
    {
        bcValue = Ui;
    }

    return bcValue;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
subsonicOutletFvPatchField<Type>::subsonicOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    pExit_(nullptr),
    curTimeIndex_(-1)
{}


template<class Type>
subsonicOutletFvPatchField<Type>::subsonicOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    pExit_
    (
        PatchFunction1<scalar>::New(p.patch(), "pExit", dict)
    ),
    curTimeIndex_(-1)
{
    if (dict.found("value"))
    {
        fixedValueFvPatchField<Type>::operator==
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        this->evaluate();
    }
}


template<class Type>
subsonicOutletFvPatchField<Type>::subsonicOutletFvPatchField
(
    const subsonicOutletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    pExit_(ptf.pExit_.clone(p.patch())),
    curTimeIndex_(-1)
{}


template<class Type>
subsonicOutletFvPatchField<Type>::subsonicOutletFvPatchField
(
    const subsonicOutletFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    pExit_(ptf.pExit_.clone(this->patch().patch())),
    curTimeIndex_(-1)
{}


template<class Type>
subsonicOutletFvPatchField<Type>::subsonicOutletFvPatchField
(
    const subsonicOutletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    pExit_(ptf.pExit_.clone(this->patch().patch())),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void subsonicOutletFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);

    if (pExit_)
    {
        pExit_->autoMap(m);
    }
}


template<class Type>
void subsonicOutletFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const subsonicOutletFvPatchField<Type>& tiptf =
        refCast<const subsonicOutletFvPatchField<Type> >(ptf);

    if (pExit_)
    {
        pExit_->rmap(tiptf.pExit_(), addr);
    }
}


template<class Type>
void subsonicOutletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        Field<Type>& patchField = *this;

        const fvPatch& patch = this->patch();

        static word thermoName(fluidThermo::dictName);

        patchField = this->bcValue
        (
            patch.template lookupPatchField<volScalarField, scalar>("p"),
            patch.template lookupPatchField<volVectorField, vector>("U"),
            this->db().objectRegistry::template lookupObject<fluidThermo>
            (
                thermoName
            )
        );

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void subsonicOutletFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    if (pExit_)
    {
        pExit_->writeData(os);
    }

    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

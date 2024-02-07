/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 OpenCFD Ltd.
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

#include "exprValue.H"
#include "ITstream.H"
#include "Switch.H"
#include <cstring>  // For memcpy, memset

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
static void fillTokens(const Type& val, tokenList& toks)
{
    const direction nCmpt = pTraits<Type>::nComponents;
    const direction nParen = 2*(pTraits<Type>::rank || (nCmpt > 1) ? 1 : 0);

    toks.resize_nocopy(nCmpt + nParen);

    auto iter = toks.begin();

    if (nParen)
    {
        *iter = token::BEGIN_LIST;
        ++iter;
    }

    for (direction cmpt = 0; cmpt < nCmpt; ++cmpt)
    {
        *iter = component(val, cmpt);
        ++iter;
    }

    if (nParen)
    {
        *iter = token::END_LIST;
        ++iter;
    }
}


//- Specialized for bool
template<>
void fillTokens<bool>(const bool& val, tokenList& toks)
{
    toks.resize_nocopy(1);
    toks.front() = token::boolean(val);
}


//- Specialized for label
template<>
void fillTokens<label>(const label& val, tokenList& toks)
{
    toks.resize_nocopy(1);
    toks.front() = val;
}


//- Specialized for scalar
template<>
void fillTokens<scalar>(const scalar& val, tokenList& toks)
{
    toks.resize_nocopy(1);
    toks.front() = val;
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::expressions::valueTypeCode
Foam::expressions::exprValue::peekType(const ITstream& is)
{
    expressions::valueTypeCode whichCode(expressions::valueTypeCode::INVALID);

    const token& tok0 = is.peek();

    if (tok0.isPunctuation(token::BEGIN_LIST))
    {
        // Expecting "( content )" - eg, (x y z), (xx xy ...)

        // First component starts after the opening '('.
        // Can use the current index if the '(' actually came from
        // the putBack.

        const label firstCmpti = (is.tokenIndex() + (is.hasPutback() ? 0 : 1));

        // Search for closing ')', require all components to be numbers
        for (label endCmpti = firstCmpti; endCmpti < is.size(); ++endCmpti)
        {
            const token& tok = is[endCmpti];

            ///InfoErr
            ///    << "check token: " << (endCmpti-firstCmpti) << " : "
            ///    << is[endCmpti].name() << nl;

            if (tok.isPunctuation(token::END_LIST))
            {
                // Select based on the number of components
                // cf. pTraits<Type>::nComponents

                switch (endCmpti - firstCmpti)
                {
                    case 0:  // Explicitly provided '()' - ie, none
                        whichCode = expressions::valueTypeCode::NONE;
                        break;

                    case 1:  // pTraits<sphericalTensor>::nComponents
                        whichCode = exprTypeTraits<sphericalTensor>::value;
                        break;

                    // FUTURE?
                    // case 2:  // pTraits<complex>::nComponents
                    //     whichCode = exprTypeTraits<complex>::value;
                    //     break;

                    case 3:  // pTraits<vector>::nComponents
                        whichCode = exprTypeTraits<vector>::value;
                        break;

                    case 6:  // pTraits<symmTensor>::nComponents
                        whichCode = exprTypeTraits<symmTensor>::value;
                        break;

                    case 9:  // pTraits<tensor>::nComponents
                        whichCode = exprTypeTraits<tensor>::value;
                        break;
                }

                // Closing ')' terminates peeking
                break;
            }
            else if (!tok.isNumber())
            {
                // All components should be numeric
                break;
            }
        }
    }
    else if (tok0.good())
    {
        /// InfoErr<< "check token: " << tok0.info() << nl;

        if (tok0.isScalar())
        {
            whichCode = exprTypeTraits<scalar>::value;
        }
        else if (tok0.isLabel())
        {
            whichCode = exprTypeTraits<label>::value;
        }
        else if (Switch(tok0).good())
        {
            whichCode = exprTypeTraits<bool>::value;
        }

        // Treat anything else as 'invalid', which also implicitly
        // includes the token "bad"
        // else if (tok0.isWord("bad"))
        // {
        //     whichCode = expressions::valueTypeCode::INVALID;
        // }
    }

    return whichCode;
}


bool Foam::expressions::exprValue::read
(
    const std::string& str,
    exprValue& val
)
{
    ITstream is(str);

    // No trailing non-whitespace!
    return (val.readTokens(is) && !is.nRemainingTokens());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::expressions::exprValue::clear()
{
    std::memset(static_cast<void*>(this), '\0', sizeof(*this));
    // Redundant:  typeCode_ = expressions::valueTypeCode::NONE;
}


void Foam::expressions::exprValue::deepCopy(const exprValue& rhs)
{
    if (this != &rhs)
    {
        // Self-assignment is a no-op
        std::memcpy(static_cast<void*>(this), &rhs, sizeof(*this));
    }
}


Foam::tokenList Foam::expressions::exprValue::tokens(bool prune) const
{
    // Handling for NONE, INVALID:
    //   - NONE    => pair of ( ) brackets
    //   - INVALID => "bad" as a word
    //
    // With prune:
    //   - no output for either

    tokenList toks;

    if (!prune)
    {
        if (typeCode_ == expressions::valueTypeCode::NONE)
        {
            toks.resize(2);
            toks.front() = token::BEGIN_LIST;
            toks.back() = token::END_LIST;
            return toks;
        }
        else if (typeCode_ == expressions::valueTypeCode::INVALID)
        {
            toks.emplace_back(word("bad"));
            return toks;
        }
    }

    switch (typeCode_)
    {
        #undef  doLocalCode
        #define doLocalCode(Type, UnusedParam)                        \
                                                                      \
        case expressions::valueTypeCode::type_##Type :                \
        {                                                             \
            const Type* dataPtr = data_.get<Type>();                  \
            if (dataPtr)                                              \
            {                                                         \
                fillTokens<Type>(*dataPtr, toks);                     \
            }                                                         \
            break;                                                    \
        }

        FOR_ALL_EXPR_VALUE_TYPES(doLocalCode);
        #undef doLocalCode

        // exprValue may only be a subset of valueTypeCode types
        default: break;
    }

    return toks;
}


void Foam::expressions::exprValue::write(Ostream& os, bool prune) const
{
    // Handling for NONE, INVALID:
    //   - NONE    => pair of ( ) brackets
    //   - INVALID => "bad" as a word
    //
    // With prune:
    //   - no output for either

    if (!prune)
    {
        if (typeCode_ == expressions::valueTypeCode::NONE)
        {
            os << token::BEGIN_LIST << token::END_LIST;
            return;
        }
        else if (typeCode_ == expressions::valueTypeCode::INVALID)
        {
            os << word("bad");
            return;
        }
    }

    switch (typeCode_)
    {
        #undef  doLocalCode
        #define doLocalCode(Type, UnusedParam)                        \
                                                                      \
        case expressions::valueTypeCode::type_##Type :                \
        {                                                             \
            const Type* dataPtr = data_.get<Type>();                  \
            if (dataPtr)                                              \
            {                                                         \
                os << *dataPtr;                                       \
            }                                                         \
            break;                                                    \
        }

        FOR_ALL_EXPR_VALUE_TYPES(doLocalCode);
        #undef doLocalCode

        // exprValue may only be a subset of valueTypeCode types
        default: break;
    }
}


bool Foam::expressions::exprValue::read(Istream& is)
{
    ITstream* stream = dynamic_cast<ITstream*>(&is);

    // Reading via tokens - simple for now
    // Expect either a single token (scalar, label, word etc)
    // or ( ... ) content

    ITstream toks;

    if (!stream)
    {
        token tok(is);

        is.fatalCheck(FUNCTION_NAME);

        if (tok.isPunctuation(token::BEGIN_LIST))
        {
            // Expecting "( content )" - eg, (x y z), (xx xy ...)
            do
            {
                toks.add_tokens(tok);

                is >> tok;
                is.fatalCheck(FUNCTION_NAME);
            }
            while (!tok.isPunctuation(token::END_LIST));

            if (tok.isPunctuation(token::END_LIST))
            {
                toks.add_tokens(tok);
            }
        }
        else if (tok.good())
        {
            toks.add_tokens(tok);
        }

        // Truncate to number tokens read
        toks.resize(toks.tokenIndex());
        toks.seek(0);

        stream = &toks;
    }

    return readTokens(*stream);
}


bool Foam::expressions::exprValue::readTokens(ITstream& is)
{
    clear();  // type: none, value: zero

    const valueTypeCode whichCode(exprValue::peekType(is));

    if (whichCode == expressions::valueTypeCode::NONE)
    {
        typeCode_ = whichCode;
        is.skip(2);  // Skip tokens: '( )'
        return true;
    }

    // This one should be rare or even impossible
    if (whichCode == expressions::valueTypeCode::INVALID)
    {
        typeCode_ = whichCode;

        if (is.bad())
        {
            return false;
        }

        const token& tok0 = is.peek();

        if (tok0.isWord("bad"))
        {
            is.skip(1);  // Skip token: "bad"
            return true;
        }
    }

    switch (whichCode)
    {
        #undef  doLocalCode
        #define doLocalCode(Type, UnusedParam)                        \
                                                                      \
        case expressions::valueTypeCode::type_##Type :                \
        {                                                             \
            data_.set<Type>(pTraits<Type>(is));                       \
            typeCode_ = whichCode;                                    \
            return true;                                              \
        }

        FOR_ALL_EXPR_VALUE_TYPES(doLocalCode);
        #undef doLocalCode

        // exprValue may only be a subset of valueTypeCode types
        default: break;
    }

    return false;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::expressions::exprValue::operator==(const exprValue& rhs) const
{
    if (typeCode_ != rhs.typeCode_)
    {
        // Types must match
        return false;
    }
    else if (this == &rhs)
    {
        return true;
    }

    switch (typeCode_)
    {
        #undef  doLocalCode
        #define doLocalCode(Type, UnusedParam)                        \
                                                                      \
        case expressions::valueTypeCode::type_##Type :                \
        {                                                             \
            const Type* a = data_.get<Type>();                        \
            const Type* b = rhs.data_.get<Type>();                    \
            return (a && b && (*a == *b));                            \
            break;                                                    \
        }

        FOR_ALL_EXPR_VALUE_TYPES(doLocalCode);
        #undef doLocalCode

        // exprValue may only be a subset of valueTypeCode types
        default: break;
    }

    return false;
}


bool Foam::expressions::exprValue::operator<(const exprValue& rhs) const
{
    // Not yet sortable
    return false;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>
(
    Istream& is,
    expressions::exprValue& val
)
{
    val.read(is);
    return is;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const expressions::exprValue& val
)
{
    val.write(os, false);  // no pruning
    return os;
}


template<>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<expressions::exprValue>& iproxy
)
{
    const auto& val = *iproxy;

    if (val.typeCode() == expressions::valueTypeCode::NONE)
    {
        os << "none";
    }
    else if (val.typeCode() == expressions::valueTypeCode::INVALID)
    {
        os << "bad";
    }
    else
    {
        os << val.valueTypeName() << ": ";
        val.write(os);  // pruning is immaterial - !good() already handled
    }

    return os;
}


// ************************************************************************* //

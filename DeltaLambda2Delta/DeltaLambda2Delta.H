/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 Custom Development
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
    Foam::LESModels::DeltaLambda2Delta

Description
    Flow-dependent length scale based on the eigenvector of the symmetric
    tensor M = S^2 + Omega^2 associated with the lambda-2 vortex criterion
    (Jeong & Hussain, 1995).

    The local vortex rotation axis n_R is taken as the eigenvector
    corresponding to the largest eigenvalue mu_1 of M:

    \verbatim
        M = symm(S & S + Omega & Omega)
        [mu_1, mu_2, mu_3]  = eigenValues(M)  with mu_1 >= mu_2 >= mu_3
        n_R = eigenvector of M for mu_1  (first row of eigenVectors(M))
    \endverbatim

    Note: OpenFOAM's eigenValues(symmTensor) returns eigenvalues in
    ascending order (mu_3 <= mu_2 <= mu_1), so the largest eigenvalue
    is component z (index 2) and its eigenvector is the third row of
    eigenVectors(M).

    Because M is always real symmetric, n_R is always well-defined —
    no fallback is required.  For maximum physical relevance, restrict
    results to cells where mu_2 < 0 (confirmed lambda-2 vortex core);
    outside those cells the axis direction is geometrically meaningful
    but not necessarily physically significant.

    Delta_R is the maximum projected cell width in the swirling plane:

    \verbatim
        Delta_R = deltaCoeff * max_i( 2 * mag(n_R ^ (fc_i - cc)) )
    \endverbatim

    Reference:
    \verbatim
        Jeong, J., & Hussain, F. (1995).
        On the identification of a vortex.
        Journal of Fluid Mechanics, 285, 69-94.
        DOI:10.1017/S0022112095000462
    \endverbatim

    Proposal source: propose_formula/proposal_04_lambda2.md
    Style reference: LESdeltas/DeltaOmegaTildeDelta/DeltaOmegaTildeDelta.H
    Additional reference for eigenvalue usage: src/functionObjects/field/Lambda2/Lambda2.C

SourceFiles
    DeltaLambda2Delta.C

\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
  Proposal  : Eigenvector of S²+Ω² Tensor (λ2 Criterion Axis)
  Source    : propose_formula/proposal_04_lambda2.md
  Reference : LESdeltas/DeltaOmegaTildeDelta/DeltaOmegaTildeDelta.H
              src/functionObjects/field/Lambda2/Lambda2.C
  Notes     : OpenFOAM's eigenValues(symmTensor) returns ascending order;
              largest eigenvalue is component z (index 2).
              The eigenvector for mu_1 is the third row of eigenVectors(M).
              No fallback needed — M is always real symmetric.
\*---------------------------------------------------------------------------*/

#ifndef LESModels_DeltaLambda2Delta_H
#define LESModels_DeltaLambda2Delta_H

#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                    Class DeltaLambda2Delta Declaration
\*---------------------------------------------------------------------------*/

class DeltaLambda2Delta
:
    public LESdelta
{
    // Private Data

        //- Model coefficient (default 1.035, matching DeltaOmegaTildeDelta)
        scalar deltaCoeff_;


    // Private Member Functions

        //- Calculate the delta values
        void calcDelta();

        //- No copy construct
        DeltaLambda2Delta(const DeltaLambda2Delta&) = delete;

        //- No copy assignment
        void operator=(const DeltaLambda2Delta&) = delete;


public:

    //- Runtime type information
    TypeName("DeltaLambda2Delta");


    // Constructors

        //- Construct from name, turbulenceModel and dictionary
        DeltaLambda2Delta
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary&
        );


    //- Destructor
    virtual ~DeltaLambda2Delta() = default;


    // Member Functions

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&);

        //- Correct delta values
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

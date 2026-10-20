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
    Foam::LESModels::DeltaSwirlDelta

Description
    Flow-dependent length scale based on the swirling-strength (lambda_ci)
    vortex identification criterion.

    When the velocity gradient tensor has a complex-conjugate eigenvalue pair
    (lambda_ci > 0), the corresponding real eigenvector identifies the local
    vortex rotation axis n_R free from pure-shear contamination.  Delta_R is
    then the maximum projected cell width in the swirling plane:

    \verbatim
        Delta_R = deltaCoeff * max_i( 2 * mag(n_R ^ (fc_i - cc)) )
    \endverbatim

    where fc_i are face-centre vectors and cc is the cell centre.  In the
    absence of swirl (lambda_ci < SMALL), the method falls back to the
    maximum face-centre distance from the cell centre.

    Reference:
    \verbatim
        Zhou, J., Adrian, R. J., Balachandar, S., & Kendall, T. M. (1999).
        Mechanisms for generating coherent packets of hairpin vortices in
        channel flow.  Journal of Fluid Mechanics, 387, 353-396.
        DOI:10.1017/S002211209900467X

        Chakraborty, P., Balachandar, S., & Adrian, R. J. (2005).
        On the relationships between local vortex identification schemes.
        Journal of Fluid Mechanics, 535, 189-214.
        DOI:10.1017/S0022112005004726
    \endverbatim

    Proposal source: propose_formula/proposal_02_swirling_strength.md
    Style reference: DeltaOmegaTildeDelta (LESdeltas/DeltaOmegaTildeDelta)

SourceFiles
    DeltaSwirlDelta.C

\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
  Proposal  : Swirling Strength / lambda_ci Rotation Axis
  Source    : propose_formula/proposal_02_swirling_strength.md
  Reference : LESdeltas/DeltaOmegaTildeDelta/DeltaOmegaTildeDelta.H
  Notes     : Uses complex eigendecomposition of grad(U); falls back to
              max face distance when lambda_ci < SMALL (no swirl detected).
\*---------------------------------------------------------------------------*/

#ifndef LESModels_DeltaSwirlDelta_H
#define LESModels_DeltaSwirlDelta_H

#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                      Class DeltaSwirlDelta Declaration
\*---------------------------------------------------------------------------*/

class DeltaSwirlDelta
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
        DeltaSwirlDelta(const DeltaSwirlDelta&) = delete;

        //- No copy assignment
        void operator=(const DeltaSwirlDelta&) = delete;


public:

    //- Runtime type information
    TypeName("DeltaSwirlDelta");


    // Constructors

        //- Construct from name, turbulenceModel and dictionary
        DeltaSwirlDelta
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary&
        );


    //- Destructor
    virtual ~DeltaSwirlDelta() = default;


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

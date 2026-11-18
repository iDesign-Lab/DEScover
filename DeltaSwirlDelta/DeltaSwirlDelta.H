/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Upstream CFD GmbH
    Copyright (C) 2022 OpenCFD Ltd.
    Copyright (C) 2026 Jiawei Chen, Xiang I. A. Yang, Ping He
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
    Delta formulation that accounts for the orientation of the local swirl
    axis given by the swirling-strength (lambda_ci) criterion. It is the
    DeltaOmegaTilde length scale of Shur et al. (2015) with the vorticity
    direction replaced by the swirl axis n_R:

    \verbatim
        delta = deltaCoeff*max_f(2*mag(n_R ^ (fc_f - cc)))
    \endverbatim

    where fc_f are the face centres and cc is the cell centre of the cell.

    n_R is the real (right) eigenvector of the velocity-gradient tensor
    A_ij = du_i/dx_j in cells where A has a complex-conjugate eigenvalue
    pair lambda_cr +/- i*lambda_ci (Zhou et al. 1999; Chakraborty et al.
    2005). Note that fvc::grad(U) in OpenFOAM is the transpose of A.

    Where no swirl exists (lambda_ci <= 1e-6*mag(A)), the length scale falls
    back to deltaCoeff*max_f(2*mag(fc_f - cc)).

    Usage
    \verbatim
    hmaxOuter       DeltaSwirlDelta;

    DeltaSwirlDeltaCoeffs
    {
        deltaCoeff      1.035;      // default 1.035
    }
    \endverbatim

    References:
    \verbatim
        Shur, M. L., Spalart, P. R., Strelets, M. K., & Travin, A. K. (2015).
        An enhanced version of DES with rapid transition
        from RANS to LES in separated flows.
        Flow, turbulence and combustion, 95(4), 709-737.
        DOI:10.1007/s10494-015-9618-0

        Zhou, J., Adrian, R. J., Balachandar, S., & Kendall, T. M. (1999).
        Mechanisms for generating coherent packets of hairpin vortices in
        channel flow. Journal of Fluid Mechanics, 387, 353-396.
        DOI:10.1017/S002211209900467X

        Chakraborty, P., Balachandar, S., & Adrian, R. J. (2005).
        On the relationships between local vortex identification schemes.
        Journal of Fluid Mechanics, 535, 189-214.
        DOI:10.1017/S0022112005004726
    \endverbatim

SourceFiles
    DeltaSwirlDelta.C

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

        //- Model coefficient (default 1.035, as in DeltaOmegaTildeDelta)
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

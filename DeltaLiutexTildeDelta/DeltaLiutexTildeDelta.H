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
    Foam::LESModels::DeltaLiutexTildeDelta

Description
    Delta formulation that accounts for the orientation of the Liutex
    (Rortex) vector, i.e. the local axis of rigid rotation.

    \verbatim
        delta = deltaCoeff*max_f(2*mag(n_R ^ (fc_f - cc)))
    \endverbatim

    where fc_f are the face centres and cc is the cell centre of the cell.

    n_R is the real (right) eigenvector of the velocity-gradient tensor
    A_ij = du_i/dx_j, sign-aligned with the vorticity, in cells where A has
    one real eigenvalue lambda_r and a complex-conjugate pair
    lambda_cr +/- i*lambda_ci (Liu et al. 2018; Gao & Liu 2018). The
    eigenvalues are obtained from the characteristic equation with
    Cardano's formula and n_R as the null vector of A - lambda_r*I.
    Note that fvc::grad(U) in OpenFOAM is the transpose of A.

    Where no rigid rotation exists (three real eigenvalues, zero Liutex
    magnitude, or lambda_ci <= 1e-6*mag(A)), n_R = 0 and hence delta = 0,
    in the same way as DeltaOmegaTilde gives delta = 0 for zero vorticity.
    The length scale is intended to be used as the outer length scale
    (hmaxOuter) of the IDDESDelta provided with this code, which bounds it
    from below.

    Usage
    \verbatim
    hmaxOuter       DeltaLiutexTilde;

    DeltaLiutexTildeCoeffs
    {
        deltaCoeff      1.035;      // default 1.035
    }
    \endverbatim

    References:
    \verbatim
        Liu, C., Gao, Y., Tian, S., & Dong, X. (2018).
        Rortex - A new vortex vector definition and vorticity tensor and
        vector decompositions.
        Physics of Fluids, 30(3), 035103.
        DOI:10.1063/1.5023001

        Gao, Y., & Liu, C. (2018).
        Rortex and comparison with eigenvalue-based vortex identification
        criteria.
        Physics of Fluids, 30(8), 085107.
        DOI:10.1063/1.5040112
    \endverbatim

SourceFiles
    DeltaLiutexTildeDelta.C

\*---------------------------------------------------------------------------*/

#ifndef LESModels_DeltaLiutexTildeDelta_H
#define LESModels_DeltaLiutexTildeDelta_H

#include "LESdelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                    Class DeltaLiutexTildeDelta Declaration
\*---------------------------------------------------------------------------*/

class DeltaLiutexTildeDelta
:
    public LESdelta
{
    // Private Data

        //- Run-time selectable delta for hmax
        //  Defaults to the maxDeltaxyz model if not supplied
        autoPtr<LESdelta> hmaxPtr_;

        //- Model coefficient
        scalar deltaCoeff_;

        //- Flag to indicate whether hmax requires updating
        bool requireUpdate_;


    // Private Member Functions

        //- Calculate the delta values
        void calcDelta();

        //- No copy construct
        DeltaLiutexTildeDelta(const DeltaLiutexTildeDelta&) = delete;

        //- No copy assignment
        void operator=(const DeltaLiutexTildeDelta&) = delete;


public:

    //- Runtime type information
    TypeName("DeltaLiutexTilde");


    // Constructors

        //- Construct from name, turbulenceModel and dictionary
        DeltaLiutexTildeDelta
        (
            const word& name,
            const turbulenceModel& turbulence,
            const dictionary&
        );


    //- Destructor
    virtual ~DeltaLiutexTildeDelta() = default;


    // Member Functions

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&);

        //- Return the hmax delta field
        const volScalarField& hmax() const
        {
            return hmaxPtr_();
        }

        // Correct values
        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

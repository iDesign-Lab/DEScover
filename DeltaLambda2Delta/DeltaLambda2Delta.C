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

\*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*\
  Proposal  : Eigenvector of S²+Ω² Tensor (λ2 Criterion Axis)
  Source    : propose_formula/proposal_04_lambda2.md
  Reference : LESdeltas/DeltaOmegaTildeDelta/DeltaOmegaTildeDelta.C
              src/functionObjects/field/Lambda2/Lambda2.C
  Notes     : Algorithm:
              1. Decompose grad(U) into strain S = symm(gradU) and
                 rotation Omega = skew(gradU) (OpenFOAM convention).
              2. Form M = symm(S & S + Omega & Omega)  (symmTensor).
              3. eigenValues(M) returns ascending real vector
                 [mu_3, mu_2, mu_1]; largest is component z.
              4. eigenVectors(M, eVals) returns tensor whose rows are
                 the eigenvectors; row z corresponds to mu_1.
              5. n_R = row z of eigenVectors(M)  (already unit length).
              6. Delta_R = max_faces( 2*|n_R ^ (fc-cc)| ).
              No fallback is needed — M is always real symmetric.
\*---------------------------------------------------------------------------*/

#include "DeltaLambda2Delta.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(DeltaLambda2Delta, 0);
    addToRunTimeSelectionTable(LESdelta, DeltaLambda2Delta, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::DeltaLambda2Delta::calcDelta()
{
    const fvMesh& mesh = turbulenceModel_.mesh();

    const volVectorField& U = turbulenceModel_.U();

    // Velocity gradient tensor
    tmp<volTensorField> tgradU = fvc::grad(U);
    const volTensorField& gradU = tgradU.cref();

    const cellList& cells          = mesh.cells();
    const vectorField& cellCentres = mesh.cellCentres();
    const vectorField& faceCentres = mesh.faceCentres();

    forAll(cells, celli)
    {
        const tensor& gradUc = gradU[celli];

        // ----------------------------------------------------------------
        // Step 1: Strain-rate and rotation-rate tensors
        //   symm(gradU) = S  (symmetric part)
        //   skew(gradU) = -Omega  (OpenFOAM skew = antisymmetric part with
        //                         skew_ij = (gradU_ij - gradU_ji)/2, which
        //                         equals -Omega_ij in the classical sign)
        //   For M = S^2 + Omega^2 we need Omega^2; note that
        //   skew(A) & skew(A) == Omega & Omega because
        //   (-Omega)(-Omega) = Omega^2.
        // ----------------------------------------------------------------
        const symmTensor S     = symm(gradUc);
        const tensor     W     = skew(gradUc);   // W = -Omega (sign-neutral
                                                  // for squared product)

        // M = symm( S*S + W*W )
        // Both S*S and W*W produce symmetric results; taking symm() makes
        // floating-point symmetry exact.
        const symmTensor M = symm
        (
            (S & S) + (W & W)
        );

        // ----------------------------------------------------------------
        // Step 2: Eigenvalues of M (ascending order: mu_3 <= mu_2 <= mu_1)
        //   In OpenFOAM: eigenValues(symmTensor) -> vector [z=min ... x=max]
        //   Wait — OpenFOAM returns ascending order, component layout:
        //     eVals[0] = smallest  (vector::X)
        //     eVals[1] = middle    (vector::Y)  <- lambda_2
        //     eVals[2] = largest   (vector::Z)  <- mu_1 (rotation axis)
        // ----------------------------------------------------------------
        const vector eVals = eigenValues(M);

        // ----------------------------------------------------------------
        // Step 3: Eigenvectors — eigenVectors(symmTensor, vector) returns
        //   a tensor whose rows are unit eigenvectors in the same order as
        //   the eigenvalues supplied.
        //   Row z (eVecs.z()) corresponds to eVals[2] = mu_1.
        // ----------------------------------------------------------------
        const tensor eVecs = eigenVectors(M, eVals);

        // Unit vector along the principal rotation axis (mu_1)
        // eVecs.z() is the third row of the eigenvector tensor
        const vector nR = eVecs.z();

        // ----------------------------------------------------------------
        // Step 4: Delta_R = max over faces of 2*|n_R ^ (fc - cc)|
        // ----------------------------------------------------------------
        const labelList& cFaces = cells[celli];
        const point& cc = cellCentres[celli];

        scalar deltaMaxTmp = 0;
        for (const label facei : cFaces)
        {
            const scalar tmp =
                2.0*mag(nR ^ (faceCentres[facei] - cc));
            if (tmp > deltaMaxTmp)
            {
                deltaMaxTmp = tmp;
            }
        }

        delta_[celli] = deltaCoeff_*deltaMaxTmp;
    }

    const label nD = mesh.nGeometricD();

    if (nD == 2)
    {
        WarningInFunction
            << "Case is 2D, LES is not strictly applicable" << nl
            << endl;
    }
    else if (nD != 3)
    {
        FatalErrorInFunction
            << "Case must be either 2D or 3D" << exit(FatalError);
    }

    delta_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::DeltaLambda2Delta::DeltaLambda2Delta
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict
)
:
    LESdelta(name, turbulence),
    deltaCoeff_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "deltaCoeff",
            1.035
        )
    )
{
    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::DeltaLambda2Delta::read(const dictionary& dict)
{
    const dictionary& coeffsDict = dict.optionalSubDict(type() + "Coeffs");
    coeffsDict.readIfPresent<scalar>("deltaCoeff", deltaCoeff_);
    calcDelta();
}


void Foam::LESModels::DeltaLambda2Delta::correct()
{
    calcDelta();
}


// ************************************************************************* //

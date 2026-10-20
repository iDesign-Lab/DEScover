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

/*---------------------------------------------------------------------------*\
  Proposal  : Swirling Strength / lambda_ci Rotation Axis
  Source    : propose_formula/proposal_02_swirling_strength.md
  Reference : LESdeltas/DeltaOmegaTildeDelta/DeltaOmegaTildeDelta.C
  Notes     : Cell-by-cell eigendecomposition of the (non-symmetric) velocity
              gradient tensor.  The real eigenvector of grad(U) is taken as
              the local swirl axis n_R when the imaginary eigenvalue
              lambda_ci > SMALL.  Falls back to maximum face distance
              otherwise.
\*---------------------------------------------------------------------------*/

#include "DeltaSwirlDelta.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(DeltaSwirlDelta, 0);
    addToRunTimeSelectionTable(LESdelta, DeltaSwirlDelta, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::DeltaSwirlDelta::calcDelta()
{
    const fvMesh& mesh = turbulenceModel_.mesh();

    // Velocity gradient tensor field
    const volVectorField& U = turbulenceModel_.U();
    tmp<volTensorField> tgradU = fvc::grad(U);
    const volTensorField& gradU = tgradU.cref();

    const cellList& cells         = mesh.cells();
    const vectorField& cellCentres = mesh.cellCentres();
    const vectorField& faceCentres = mesh.faceCentres();

    forAll(cells, celli)
    {
        const tensor& gradUc = gradU[celli];

        // -----------------------------------------------------------------
        // Step 1: Complex eigenvalues of grad(U)
        //   eigenValues returns Vector<complex> with three entries.
        //   For one-real + two-complex-conjugate case, one entry has
        //   near-zero imaginary part (the real eigenvalue) and the other
        //   two have ±lambda_ci imaginary parts.
        // -----------------------------------------------------------------
        const Vector<complex> eVals = eigenValues(gradUc);

        const scalar lciX = mag(eVals.x().imag());
        const scalar lciY = mag(eVals.y().imag());
        const scalar lciZ = mag(eVals.z().imag());
        const scalar lambdaCi = max(lciX, max(lciY, lciZ));

        const labelList& cFaces = cells[celli];
        const point& cc = cellCentres[celli];

        // Fallback: max face half-width (no preferred direction)
        scalar deltaMaxTmp = 0;
        for (const label facei : cFaces)
        {
            const scalar tmp =
                2.0*mag(faceCentres[facei] - cc);
            if (tmp > deltaMaxTmp)
            {
                deltaMaxTmp = tmp;
            }
        }

        if (lambdaCi < SMALL)
        {
            // No meaningful swirl detected — use isotropic fallback
            delta_[celli] = deltaCoeff_*deltaMaxTmp;
            continue;
        }

        // -----------------------------------------------------------------
        // Step 2: Identify index of the real eigenvalue (smallest |imag|)
        // -----------------------------------------------------------------
        label realIdx = 0;
        scalar minImag = lciX;
        if (lciY < minImag)
        {
            minImag = lciY;
            realIdx = 1;
        }
        if (lciZ < minImag)
        {
            realIdx = 2;
        }

        // -----------------------------------------------------------------
        // Step 3: Extract real eigenvector corresponding to lambda_r
        //   eigenVectors returns Tensor<complex>; each row is the
        //   eigenvector for the eigenvalue at the same index.
        // -----------------------------------------------------------------
        const Tensor<complex> eVecs = eigenVectors(gradUc, eVals);

        vector nR(Zero);
        switch (realIdx)
        {
            case 0:
            {
                const Vector<complex>& row = eVecs.x();
                nR = vector
                (
                    row.x().real(),
                    row.y().real(),
                    row.z().real()
                );
                break;
            }
            case 1:
            {
                const Vector<complex>& row = eVecs.y();
                nR = vector
                (
                    row.x().real(),
                    row.y().real(),
                    row.z().real()
                );
                break;
            }
            default:  // case 2
            {
                const Vector<complex>& row = eVecs.z();
                nR = vector
                (
                    row.x().real(),
                    row.y().real(),
                    row.z().real()
                );
                break;
            }
        }

        const scalar nRmag = mag(nR);
        if (nRmag < SMALL)
        {
            // Numerically degenerate eigenvector — fall back
            delta_[celli] = deltaCoeff_*deltaMaxTmp;
            continue;
        }
        nR /= nRmag;

        // -----------------------------------------------------------------
        // Step 4: Delta_R = max over faces of 2*|n_R ^ (fc - cc)|
        //   This gives the maximum cell width in the swirling plane.
        // -----------------------------------------------------------------
        deltaMaxTmp = 0;
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

Foam::LESModels::DeltaSwirlDelta::DeltaSwirlDelta
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

void Foam::LESModels::DeltaSwirlDelta::read(const dictionary& dict)
{
    const dictionary& coeffsDict = dict.optionalSubDict(type() + "Coeffs");
    coeffsDict.readIfPresent<scalar>("deltaCoeff", deltaCoeff_);
    calcDelta();
}


void Foam::LESModels::DeltaSwirlDelta::correct()
{
    calcDelta();
}


// ************************************************************************* //

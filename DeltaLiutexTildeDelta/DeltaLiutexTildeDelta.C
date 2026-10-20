/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Upstream CFD GmbH
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "DeltaLiutexTildeDelta.H"
#include "fvcCurl.H"
#include "maxDeltaxyz.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(DeltaLiutexTildeDelta, 0);
    addToRunTimeSelectionTable(LESdelta, DeltaLiutexTildeDelta, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::DeltaLiutexTildeDelta::calcDelta()
{
    const fvMesh& mesh = turbulenceModel_.mesh();
    const volVectorField& U0 = turbulenceModel_.U();
    tmp<volVectorField> tvorticity = fvc::curl(U0);
    const volVectorField& vorticity = tvorticity.cref();

    // Start of Liutex calculation
    tmp<volTensorField> tgradUTensor = fvc::grad(U0);
    const volTensorField& gradUTensor = tgradUTensor.cref();

    volScalarField imageEigenValueField
    (
        IOobject
        (
            "imageEigenValueField",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless/dimTime, 0.0)
    );    
    volVectorField realEigenVectorField
    (
        IOobject
        (
            "realEigenVectorField",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("realEigenVectorField", dimless, Zero)
    );  
    //通过对每一个cell操作，得到
    forAll(gradUTensor, cellI)
    {
        scalarSquareMatrix gradUCell(3, 0);  // 初始化一个 3x3 的零矩阵
        tensor gradU = gradUTensor[cellI];  // 获取一个单元格的梯度
        gradUCell(0, 0) = gradU[0];        // 逐元素填充矩阵
        gradUCell(0, 1) = gradU[1];
        gradUCell(0, 2) = gradU[2];
        gradUCell(1, 0) = gradU[3];
        gradUCell(1, 1) = gradU[4];
        gradUCell(1, 2) = gradU[5];
        gradUCell(2, 0) = gradU[6];
        gradUCell(2, 1) = gradU[7];
        gradUCell(2, 2) = gradU[8];

        scalar A11 = gradU[0];        // 逐元素填充矩阵
        scalar A12 = gradU[1];
        scalar A13 = gradU[2];
        scalar A21 = gradU[3];
        scalar A22 = gradU[4];
        scalar A23 = gradU[5];
        scalar A31 = gradU[6];
        scalar A32 = gradU[7];
        scalar A33 = gradU[8];

        scalar I1 = A11 + A22 + A33;
        scalar I2 = (A11*A22 + A11*A33 + A22*A33) - (A12*A21 + A13*A31 + A23*A32);
        scalar I3 = A11*(A22*A33 -A23*A32)-A12*(A21*A33 - A23*A31) + A13*(A21*A32 - A22*A31);

        scalar p = (3.0*I2 -I1*I1)/3.0;
        scalar q = (2.0*I1*I1*I1 - 9.0*I1*I2 + 27.0*I3)/27.0;
       
        scalar delta = sqr(q/2.0) + pow(p/3.0,3.0);
        vector realEigenVector(0, 0, 0); // 默认旋转矢量
        scalar imageEigenValue = 0;
        if (delta > 0) 
        {
            scalar sqrtDelta = sqrt(delta);
            scalar S1 = cbrt(-q/2.0 + sqrtDelta);
            scalar S2 = cbrt(-q/2.0 - sqrtDelta);
            scalar realEigenValue = S1 + S2 + I1/3.0;
            imageEigenValue = mag(sqrt(3.0)/2.0 * (S1 - S2));

            // 计算实特征向量
            // 构造矩阵 B = A - λI
            scalar B11 = A11 - realEigenValue;
            scalar B12 = A12;
            scalar B13 = A13;
            scalar B21 = A21;
            scalar B22 = A22 - realEigenValue;
            scalar B23 = A23;
            scalar B31 = A31;
            scalar B32 = A32;
            scalar B33 = A33 - realEigenValue;
            // 假设第三个分量为1.0，求解其他分量
            // 求解 2x2 线性方程组
            scalar det = B11*B22 - B12*B21;
            // 检查行列式是否接近零，如果是，尝试其他假设
            // 求解 v1 和 v2
            scalar v1 = (B13*B22 - B12*B23) / det;
            scalar v2 = (B11*B23 - B13*B21) / det;
            scalar v3 = 1.0;
            // 归一化
            scalar norm = sqrt(v1*v1 + v2*v2 + v3*v3);
            
            realEigenVector = vector(v1/norm, v2/norm, v3/norm);
                        
        }
        realEigenVectorField[cellI] = realEigenVector;  
        imageEigenValueField[cellI] = imageEigenValue;   
    }
    // 确保realEigenVectorField与W成锐角
    const tmp<volScalarField> WdotEigtmp = vorticity & realEigenVectorField;
    const volScalarField& WdotEig = WdotEigtmp();
    // 创建一个临时场来存储修正后的特征向量场
    volVectorField correctedEigVectorField = realEigenVectorField;
    // 对于WdotEig < 0的情况（钝角），反转eigenvector的方向
    forAll(WdotEig, cellI)
    {
        if (WdotEig[cellI] < 0)
        {
            correctedEigVectorField[cellI] = -realEigenVectorField[cellI];
        }
    }
    // 现在使用修正后的特征向量场进行计算
    const tmp<volScalarField> Wrtmp = mag(vorticity & correctedEigVectorField); 
    const volScalarField& Wr = Wrtmp();
    // 计算 Wrlambda 但确保结果非负
    const tmp<volScalarField> rawWrlambdatmp = sqr(Wr) - 4.0 * sqr(imageEigenValueField);
    const volScalarField& rawWrlambda = rawWrlambdatmp();
    // 创建一个新的场，将负值截断为零
    volScalarField Wrlambda
    (
        IOobject
        (
            "Wrlambda",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", rawWrlambda.dimensions(), 0.0)
    );
    // 遍历每个单元格，确保值非负
    forAll(Wrlambda, cellI)
    {
        Wrlambda[cellI] = max(rawWrlambda[cellI], 0.0);
    }
    tmp<volVectorField> tLiutex = (Wr - sqrt(Wrlambda)) * correctedEigVectorField;
    const volVectorField& Liutex = tLiutex.cref();
    const volVectorField nvecLiutex
    (
        Liutex
      / (
            max
            (
                mag(Liutex),
                dimensionedScalar(dimless/dimTime, SMALL)
            )
        )
    );
    tLiutex.clear();   
    // End of Liutex calculation
    const cellList& cells = mesh.cells();
    const vectorField& cellCentres = mesh.cellCentres();
    const vectorField& faceCentres = mesh.faceCentres();

    forAll(cells, celli)
    {
        const labelList& cFaces = cells[celli];
        const point& cc = cellCentres[celli];
        const vector& nv = nvecLiutex[celli];

        scalar deltaMaxTmp = 0;

        for (const label facei : cFaces)
        {
            const point& fc = faceCentres[facei];
            const scalar tmp = 2.0*mag(nv ^ (fc - cc));

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

    // Handle coupled boundaries
    delta_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::DeltaLiutexTildeDelta::DeltaLiutexTildeDelta
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict
)
:
    LESdelta(name, turbulence),
    hmaxPtr_(nullptr),
    deltaCoeff_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "deltaCoeff",
            1.035
        )
    ),
    requireUpdate_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<bool>
        (
            "requireUpdate",
            true
        )
    )
{
    if (dict.optionalSubDict(type() + "Coeffs").found("hmax"))
    {
        // User-defined hmax
        hmaxPtr_ =
            LESdelta::New
            (
                IOobject::groupName("hmax", turbulence.U().group()),
                turbulence,
                dict.optionalSubDict("hmaxCoeffs"),
                "hmax"
            );
    }
    else
    {
        Info<< "Employing " << maxDeltaxyz::typeName << " for hmax" << endl;

        hmaxPtr_.reset
        (
            new maxDeltaxyz
            (
                IOobject::groupName("hmax", turbulence.U().group()),
                turbulence,
                dict.optionalSubDict("hmaxCoeffs")
            )
        );
    }

    calcDelta();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::DeltaLiutexTildeDelta::read(const dictionary& dict)
{
    const dictionary& coeffsDict = dict.optionalSubDict(type() + "Coeffs");

    coeffsDict.readIfPresent<scalar>("deltaCoeff", deltaCoeff_);
    coeffsDict.readIfPresent<bool>("requireUpdate", requireUpdate_);

    calcDelta();
}


void Foam::LESModels::DeltaLiutexTildeDelta::correct()
{
    if (turbulenceModel_.mesh().changing() && requireUpdate_)
    {
        hmaxPtr_->correct();
    }

    calcDelta();
}


// ************************************************************************* //

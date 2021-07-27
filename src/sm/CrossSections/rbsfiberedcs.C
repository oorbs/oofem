/*
 *
 *                      ____ ____  ____  ____ _____
 *                     / __   __ \/ __ \/ __ ) ___/
 *                    / / /  / / / /_/ / __  \__ \
 *                   / /_/  /_/ / _, _/ /_/ /__/ /
 *                   \____ ____/_/ |_/_____/____/
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "sm/CrossSections/rbsfiberedcs.h"
#include "sm/CrossSections/fiberedcs.h"
#include "sm/Elements/structuralelement.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "gausspoint.h"
#include "material.h"
#include "floatarray.h"
#include "verbose.h"
#include "contextioerr.h"
#include "classfactory.h"

namespace oofem {
REGISTER_CrossSection(RBSFiberedCrossSection);

void RBSFiberedCrossSection::initializeFrom( int number, IntArray fiberMaterials,
                                             FloatArray fiberThicks, FloatArray fiberWidths, int numberOfFibers,
                                             double thick, double width, FloatArray fiberYcoords, FloatArray fiberZcoords )
{
    this->number         = number ? number : this->number;
    this->fiberMaterials = fiberMaterials.giveSize() ? fiberMaterials : this->fiberMaterials;
    this->fiberThicks    = fiberThicks.giveSize() ? fiberThicks : this->fiberThicks;
    this->fiberWidths    = fiberWidths.giveSize() ? fiberWidths : this->fiberWidths;
    this->fiberYcoords   = fiberYcoords.giveSize() ? fiberYcoords : this->fiberYcoords;
    this->fiberZcoords   = fiberZcoords.giveSize() ? fiberZcoords : this->fiberZcoords;
    this->thick          = thick ? thick : this->thick;
    this->width          = width ? width : this->width;

    int num = this->fiberMaterials.giveSize();
    if ( num != this->fiberThicks.giveSize()    ||
         num != this->fiberWidths.giveSize()    ||
         num != this->fiberYcoords.giveSize()   ||
         num != this->fiberZcoords.giveSize() ) {
        OOFEM_ERROR("%s Array size mismatch ", _IFT_FiberedCrossSection_fibermaterials);
    }

    if ( num <= 0 ) {
        OOFEM_ERROR("%s number of fibers == 0 is not allowed", _IFT_FiberedCrossSection_fibermaterials);
    }

    this->area = fiberThicks.dotProduct(fiberWidths);
}


FloatArrayF<6> RBSFiberedCrossSection::giveGeneralizedStress_Beam3d( const FloatArrayF<6> &strain,
    GaussPoint *gp, TimeStep *tStep ) const
{
    FloatArray fiberStrain;
    auto element = static_cast< StructuralElement * >( gp->giveElement() );
    auto interface = static_cast< FiberedCrossSectionInterface * >( element->giveInterface(FiberedCrossSectionInterfaceType) );

    if ( interface == nullptr ) {
        OOFEM_ERROR("element with no fiber support encountered");
    }

    FloatArrayF<6> answer;
    double shearCoef = 2.;

    for ( int i = 1; i <= this->fiberMaterials.giveSize(); i++ ) {
        auto fiberGp = this->giveSlaveGaussPoint(gp, i - 1);
        auto fiberMat = static_cast< StructuralMaterial * >( domain->giveMaterial( fiberMaterials.at(i) ) );

        double fiberThick  = this->fiberThicks.at(i);
        double fiberWidth  = this->fiberWidths.at(i);
        double fiberYCoord = fiberGp->giveNaturalCoordinate(1);
        double fiberZCoord = fiberGp->giveNaturalCoordinate(2);

        double fiberArea = fiberWidth * fiberThick;
        double shearArea = shearCoef * fiberArea;

        interface->FiberedCrossSectionInterface_computeStrainVectorInFiber(fiberStrain, strain, fiberGp, tStep);

        auto reducedFiberStress = fiberMat->giveRealStressVector_Fiber(fiberStrain, fiberGp, tStep);

        // perform integration
        //
        // strainVector3d:      {eps_x,    eps_y,    eps_z, gamma_yz,      gamma_zx, gamma_xy}
        // strainVectorShell:   {eps_x, gamma_xz, gamma_xy, \der{phi_x}{x}, kappa_y, kappa_z }
        //

        // 1) membrane terms N, Qz, Qy
        answer.at(1) += reducedFiberStress.at(1) * fiberArea;
        answer.at(2) += reducedFiberStress.at(2) * shearArea;
        answer.at(3) += reducedFiberStress.at(3) * shearArea;
        // 2) bending terms Tx, My, Mz
        answer.at(4) += reducedFiberStress.at(2) * fiberArea * fiberYCoord
                      - reducedFiberStress.at(3) * fiberArea * fiberZCoord;
        answer.at(5) += reducedFiberStress.at(1) * fiberArea * fiberZCoord;
        answer.at(6) -= reducedFiberStress.at(1) * fiberArea * fiberYCoord;
    }

    // why only the first fiber's material is used to update the master Gauss point?
    auto status = static_cast< StructuralMaterialStatus * >( domain->giveMaterial( fiberMaterials.at(1) )->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);

    return answer;
}


FloatMatrixF<6,6> RBSFiberedCrossSection::give3dBeamStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
//
{
    double Ip = 0.0, A = 0.0, Ik, G = 0.0;
    double shearCoef = 2.;

    FloatMatrixF<6,6> beamStiffness;

    // perform integration over layers
    for ( int i = 1; i <= this->fiberMaterials.giveSize(); i++ ) {
        auto fiberGp = giveSlaveGaussPoint(gp, i - 1);
        auto mat = dynamic_cast< StructuralMaterial * >( domain->giveMaterial( fiberMaterials.at( fiberGp->giveNumber() ) ) );
        auto fiberMatrix = mat->giveFiberStiffMtrx(rMode, fiberGp, tStep);

        double fiberThick  = this->fiberThicks.at(i);
        double fiberWidth  = this->fiberWidths.at(i);
        double fiberZCoord = fiberZcoords.at(i);
        double fiberYCoord = fiberYcoords.at(i);
        double fiberYCoord2 = fiberYCoord * fiberYCoord;
        double fiberZCoord2 = fiberZCoord * fiberZCoord;

        double fiberArea = fiberWidth * fiberThick;
        double shearArea = shearCoef * fiberArea;

        // perform integration

        // 1) membrane terms N, Qz, Qy

        beamStiffness.at(1, 1) += fiberMatrix.at(1, 1) * fiberArea;

        beamStiffness.at(2, 2) += fiberMatrix.at(2, 2) * shearArea;

        beamStiffness.at(3, 3) += fiberMatrix.at(3, 3) * shearArea;

        // 2) bending terms Tx, My, Mz

        Ip += fiberArea * fiberZCoord2 + fiberArea * fiberYCoord2;
        A  += fiberArea;
        G  += fiberMatrix.at(2, 2) * fiberArea;
        //GA = fiberMatrix.at(2, 2) * fiberWidth * fiberThick;

        beamStiffness.at(5, 5) += fiberMatrix.at(1, 1) * fiberArea * fiberZCoord2;
        beamStiffness.at(6, 6) += fiberMatrix.at(1, 1) * fiberArea * fiberYCoord2;
        //beamStiffness.at(4, 4) += GA * fiberZCoord2;
        //beamStiffness.at(4, 4) += GA * fiberYCoord2;
    }

    G /= A;
    Ik = A * A * A * A / ( 40.0 * Ip );
    beamStiffness.at(4, 4) = G * Ik;
    return beamStiffness;
}

} // end namespace oofem

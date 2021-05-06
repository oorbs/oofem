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

#ifndef rbsfiberedcs_h
#define rbsfiberedcs_h

#include "sm/CrossSections/fiberedcs.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "sm/Materials/structuralmaterial.h"
#include "element.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "interface.h"

///@name Input fields for RBSFiberedCrossSection
//@{
#define _IFT_RBSFiberedCrossSection_Name "rbsfiberedcs"
//#define _IFT_FiberedCrossSection_nfibers "nfibers"
//#define _IFT_FiberedCrossSection_fibermaterials "fibermaterials"
//#define _IFT_FiberedCrossSection_thicks "thicks"
//#define _IFT_FiberedCrossSection_widths "widths"
//#define _IFT_FiberedCrossSection_fiberycentrecoords "fiberycentrecoords"
//#define _IFT_FiberedCrossSection_fiberzcentrecoords "fiberzcentrecoords"
//#define _IFT_FiberedCrossSection_thick "thick"
//#define _IFT_FiberedCrossSection_width "width"
//@}

namespace oofem {
class GaussPoint;
//class FiberedCrossSectionModelInterface;

/**
 * This class extends the fibered cross section to make it suitable to be used in
 * RBSM simulations
 */
class RBSFiberedCrossSection : public FiberedCrossSection
{
protected:

public:
    RBSFiberedCrossSection(int n, Domain * d) : FiberedCrossSection(n, d) {}
    void initializeFrom( int number = 0, IntArray fiberMaterials = {},
        FloatArray fiberThicks = {}, FloatArray fiberWidths = {},
        int numberOfFibers = 0, double thick = 0.0, double width = 0.0,
        FloatArray fiberYcoords = {}, FloatArray fiberZcoords = {} ) override;


    FloatArrayF<6> giveGeneralizedStress_Beam3d(const FloatArrayF<6> &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF<6,6> give3dBeamStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    // identification and auxiliary functions
    const char *giveInputRecordName() const override { return _IFT_RBSFiberedCrossSection_Name; }
    const char *giveClassName() const override { return "RBSFiberedCrossSection"; }

};
} // end namespace oofem
#endif // rbsfiberedcs_h

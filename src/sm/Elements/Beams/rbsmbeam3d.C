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

#include "sm/Elements/Beams/rbsmbeam3d.h"
#include "sm/Materials/structuralms.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "engngm.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "fei3dlinelin.h"
#include "classfactory.h"
#include "elementinternaldofman.h"
#include "masterdof.h"
#include "bctracker.h"

#include "bodyload.h"
#include "boundaryload.h"

#include "sm/Elements/Beams/beam3d.h"

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element( RBSMBeam3d );


RBSMBeam3d::RBSMBeam3d( int n, Domain *aDomain ) : LIBeam3d( n, aDomain )
{
    numberOfGaussPoints = 1;

    referenceNode = 0;
    referenceAngle = 0;
    this->zaxis.clear();
}

RBSMBeam3d::~RBSMBeam3d()
{
}

void RBSMBeam3d::initialize()
{
    // move initialization to here from rbsmtetra.C "instantiate beam element"
}


void RBSMBeam3d::initializeFrom( InputRecord &ir )
{
    //LIBeam3d::initializeFrom( ir );
    StructuralElement :: initializeFrom(ir);
    referenceNode = 0;
    referenceAngle = 0;
    this->zaxis.clear();
    if ( ir.hasField(_IFT_Beam3d_zaxis) ) {
        IR_GIVE_FIELD(ir, this->zaxis, _IFT_Beam3d_zaxis);
    } else if ( ir.hasField(_IFT_Beam3d_refnode) ) {
        IR_GIVE_FIELD(ir, referenceNode, _IFT_Beam3d_refnode);
        if ( referenceNode == 0 ) {
            OOFEM_WARNING("wrong reference node specified. Using default orientation.");
        }
    } else if ( ir.hasField(_IFT_Beam3d_refangle) ) {
        IR_GIVE_FIELD(ir, referenceAngle, _IFT_Beam3d_refangle);
    } else {
        throw ValueInputException(ir, _IFT_Beam3d_zaxis, "axis, reference node, or angle not set");
    }
}


int RBSMBeam3d :: giveLocalCoordinateSystem(FloatMatrix &answer)
{
FloatArray lx, ly, lz, help(3);
Node *nodeA, *nodeB;
nodeA = this->giveNode(1);
nodeB = this->giveNode(2);

lx.beDifferenceOf( nodeB->giveCoordinates(), nodeA->giveCoordinates() );
lx.normalize();

if ( this->referenceNode ) {
Node *refNode = this->giveDomain()->giveNode(this->referenceNode);
help.beDifferenceOf( refNode->giveCoordinates(), nodeA->giveCoordinates() );

lz.beVectorProductOf(lx, help);
lz.normalize();
} else if ( this->zaxis.giveSize() > 0 ) {
lz = this->zaxis;
lz.add(lz.dotProduct(lx), lx);
lz.normalize();
} else {
FloatMatrix rot(3, 3);
double theta = referenceAngle * M_PI / 180.0;

rot.at(1, 1) = cos(theta) + pow(lx.at(1), 2) * ( 1 - cos(theta) );
rot.at(1, 2) = lx.at(1) * lx.at(2) * ( 1 - cos(theta) ) - lx.at(3) * sin(theta);
rot.at(1, 3) = lx.at(1) * lx.at(3) * ( 1 - cos(theta) ) + lx.at(2) * sin(theta);

rot.at(2, 1) = lx.at(2) * lx.at(1) * ( 1 - cos(theta) ) + lx.at(3) * sin(theta);
rot.at(2, 2) = cos(theta) + pow(lx.at(2), 2) * ( 1 - cos(theta) );
rot.at(2, 3) = lx.at(2) * lx.at(3) * ( 1 - cos(theta) ) - lx.at(1) * sin(theta);

rot.at(3, 1) = lx.at(3) * lx.at(1) * ( 1 - cos(theta) ) - lx.at(2) * sin(theta);
rot.at(3, 2) = lx.at(3) * lx.at(2) * ( 1 - cos(theta) ) + lx.at(1) * sin(theta);
rot.at(3, 3) = cos(theta) + pow(lx.at(3), 2) * ( 1 - cos(theta) );

help.at(3) = 1.0;         // up-vector
// here is ly is used as a temp var
if ( fabs( lx.dotProduct(help) ) > 0.999 ) { // Check if it is vertical
ly = {
    0., 1., 0.
};
} else {
ly.beVectorProductOf(lx, help);
}
lz.beProductOf(rot, ly);
lz.normalize();
}

ly.beVectorProductOf(lz, lx);
ly.normalize();

answer.resize(3, 3);
answer.zero();
for ( int i = 1; i <= 3; i++ ) {
answer.at(1, i) = lx.at(i);
answer.at(2, i) = ly.at(i);
answer.at(3, i) = lz.at(i);
}

return 1;
}


void RBSMBeam3d::FiberedCrossSectionInterface_computeStrainVectorInFiber(
    FloatArray &answer, const FloatArray &masterGpStrain,
    GaussPoint *slaveGp, TimeStep *tStep )
{
    double layerYCoord, layerZCoord;

    layerZCoord = slaveGp->giveNaturalCoordinate( 2 );
    layerYCoord = slaveGp->giveNaturalCoordinate( 1 );

    answer.resize( 3 ); // {Exx,GMzx,GMxy}

    answer.at( 1 ) = masterGpStrain.at( 1 )
                     + masterGpStrain.at( 5 ) * layerZCoord
                     - masterGpStrain.at( 6 ) * layerYCoord;
    answer.at( 2 ) = masterGpStrain.at( 2 )
                     + masterGpStrain.at( 4 ) * layerYCoord
        ;
    answer.at( 3 ) = masterGpStrain.at( 3 )
                     - masterGpStrain.at( 4 ) * layerZCoord
        ;
}


void RBSMBeam3d::setNumberOfGaussPoints( int nip )
{
    numberOfGaussPoints = nip;
}

// put this in an interface:
// (v-2.5) calculate confined stresses to apply Poisson's effect:
void RBSMBeam3d::RBSMTetraInterface_computeStressVector( FloatArray &answer, TimeStep *tStep )
//, int useUpdatedGpRecord )
{
    /// todo: separate strain calculation from stress calculation
    FloatMatrix b;
    FloatArray u, stress, strain;

    // This function can be quite costly to do inside the loops when one has many slave dofs.
    this->computeVectorOf(VM_Total, tStep, u);
    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    for ( GaussPoint *gp : *this->giveDefaultIntegrationRulePtr() ) {
        this->computeBmatrixAt( gp, b );

        if ( !this->isActivated( tStep ) ) {
            strain.resize( StructuralMaterial ::giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
            strain.zero();
        }
        strain.beProductOf( b, u );
        // change to compute confined stress components: *** <COMPUTE CONFINED STRESS COMPONENTS> ***
        // this->computeStressVector( stress, strain, gp, tStep );
        stress = this->giveStructuralCrossSection()->giveRealStress_3d( strain, gp, tStep );

        // updates gp stress and strain record according to current increment of displacement
        if ( stress.giveSize() == 0 ) {
            break;
        }
        answer = stress;
        return;
    }
    OOFEM_ERROR( "Could not obtain stress of the RBS-Mindlin beam element")
}

} // end namespace oofem
// <temporary/>
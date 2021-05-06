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

#include "sm/Elements/Beams/rbsbeam3d.h"
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

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element( RBSBeam3d );

//FEI3dLineLin Beam3d :: interp;

RBSBeam3d::RBSBeam3d( int n, Domain *aDomain ) : Beam3d( n, aDomain )
{
    numberOfGaussPoints = 1;
}

RBSBeam3d::~RBSBeam3d()
{
}

void RBSBeam3d::initializeFrom( InputRecord &ir )
{
    Beam3d::initializeFrom( ir );
}

//FEInterpolation *RBSBeam3d ::giveInterpolation() const { return &interp; }


void RBSBeam3d::computeBmatrixAt( GaussPoint *gp, FloatMatrix &answer, int li, int ui )
{
    double l, ksi, kappay, kappaz, c1y, c1z;
    TimeStep *tStep = this->domain->giveEngngModel()->giveCurrentStep();

    l      = this->computeLength();
    ksi    = 0.5 + 0.5 * gp->giveNaturalCoordinate( 1 );
    kappay = this->giveKappayCoeff( tStep );
    kappaz = this->giveKappazCoeff( tStep );
    c1y    = 1. + 2. * kappay;
    c1z    = 1. + 2. * kappaz;

    // eeps = {\eps_x, \gamma_xz, \gamma_xy, \der{phi_x}{x}, \kappa_y, \kappa_z}^T
    answer.resize( 6, 12 );
    answer.zero();

    answer.at( 1, 1 ) = -1. / l;
    answer.at( 1, 7 ) = 1. / l;

    answer.at( 2, 3 )  = ( -2. * kappay ) / ( l * c1y );
    answer.at( 2, 5 )  = kappay / ( c1y );
    answer.at( 2, 9 )  = 2. * kappay / ( l * c1y );
    answer.at( 2, 11 ) = kappay / ( c1y );

    answer.at( 3, 2 )  = ( -2. * kappaz ) / ( l * c1z );
    answer.at( 3, 6 )  = -kappaz / ( c1z );
    answer.at( 3, 8 )  = 2. * kappaz / ( l * c1z );
    answer.at( 3, 12 ) = -kappaz / ( c1z );


    answer.at( 4, 4 )  = -1. / l;
    answer.at( 4, 10 ) = 1. / l;

    answer.at( 5, 3 )  = ( 6. - 12. * ksi ) / ( l * l * c1y );
    answer.at( 5, 5 )  = ( -2. * ( 2. + kappay ) + 6. * ksi ) / ( l * c1y );
    answer.at( 5, 9 )  = ( -6. + 12. * ksi ) / ( l * l * c1y );
    answer.at( 5, 11 ) = ( -2. * ( 1. - kappay ) + 6. * ksi ) / ( l * c1y );

    answer.at( 6, 2 )  = -1.0 * ( 6. - 12. * ksi ) / ( l * l * c1z );
    answer.at( 6, 6 )  = ( -2. * ( 2. + kappaz ) + 6. * ksi ) / ( l * c1z );
    answer.at( 6, 8 )  = -1.0 * ( -6. + 12. * ksi ) / ( l * l * c1z );
    answer.at( 6, 12 ) = ( -2. * ( 1. - kappaz ) + 6. * ksi ) / ( l * c1z );
}

void RBSBeam3d ::computeNmatrixAt( const FloatArray &iLocCoord, FloatMatrix &answer )
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp. Used for numerical calculation of consistent mass
// matrix. Must contain only interpolation for displacement terms,
// not for any rotations. (Inertia forces do not work on rotations).
// r = {u1,v1,w1,fi_x1,fi_y1,fi_z1,u2,v2,w2,fi_x2,fi_y2,fi_21}^T
{
    double l, ksi, ksi2, ksi3, kappay, kappaz, c1y, c1z;
    TimeStep *tStep = this->domain->giveEngngModel()->giveCurrentStep();

    l      = this->computeLength();
    ksi    = 0.5 + 0.5 * iLocCoord.at( 1 );
    kappay = this->giveKappayCoeff( tStep );
    kappaz = this->giveKappazCoeff( tStep );
    c1y    = 1. + 2. * kappay;
    c1z    = 1. + 2. * kappaz;
    ksi2   = ksi * ksi;
    ksi3   = ksi2 * ksi;

    answer.resize( 6, 12 );
    answer.zero();

    answer.at( 1, 1 )  = 1. - ksi;
    answer.at( 1, 7 )  = ksi;
    answer.at( 2, 2 )  = ( c1z - 2. * kappaz * ksi - 3. * ksi2 + 2. * ksi3 ) / c1z;
    answer.at( 2, 6 )  = -l * ( -( 1. + kappaz ) * ksi + ( 2. + kappaz ) * ksi2 - ksi3 ) / c1z;
    answer.at( 2, 8 )  = ( 2. * kappaz * ksi + 3. * ksi2 - 2. * ksi3 ) / c1z;
    answer.at( 2, 12 ) = -l * ( kappaz * ksi + ( 1. - kappaz ) * ksi2 - ksi3 ) / c1z;
    answer.at( 3, 3 )  = ( c1y - 2. * kappay * ksi - 3. * ksi2 + 2. * ksi3 ) / c1y;
    answer.at( 3, 5 )  = l * ( -( 1. + kappay ) * ksi + ( 2. + kappay ) * ksi2 - ksi3 ) / c1y;
    answer.at( 3, 9 )  = ( 2. * kappay * ksi + 3. * ksi2 - 2. * ksi3 ) / c1y;
    answer.at( 3, 11 ) = l * ( kappay * ksi + ( 1. - kappay ) * ksi2 - ksi3 ) / c1y;

    // rotations excluded
    answer.at( 4, 4 )  = 1. - ksi;
    answer.at( 4, 10 ) = ksi;
    answer.at( 5, 3 )  = ( 6. * ksi - 6. * ksi2 ) / ( l * c1y );
    answer.at( 5, 5 )  = ( c1y - 2. * ( 2. + kappay ) * ksi + 3. * ksi2 ) / c1y;
    answer.at( 5, 9 )  = -( 6. * ksi - 6. * ksi2 ) / ( l * c1y );
    answer.at( 5, 11 ) = ( -2. * ( 1. - kappay ) * ksi + 3. * ksi2 ) / c1y;
    answer.at( 6, 2 )  = -( 6. * ksi - 6. * ksi2 ) / ( l * c1z );
    answer.at( 6, 6 )  = ( c1z - 2. * ( 2. + kappaz ) * ksi + 3. * ksi2 ) / c1z;
    answer.at( 6, 8 )  = ( 6. * ksi - 6. * ksi2 ) / ( l * c1z );
    answer.at( 6, 12 ) = ( -2. * ( 1. - kappaz ) * ksi + 3. * ksi2 ) / c1z;
}


//void RBSBeam3d ::computeStiffnessMatrix( FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep ){}

//double RBSBeam3d ::computeVolumeAround( GaussPoint *gp ) {}


void RBSBeam3d ::computeConstitutiveMatrixAt( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
{
    answer = this->giveStructuralCrossSection()->give3dBeamStiffMtrx( rMode, gp, tStep );
}


void RBSBeam3d ::computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep )
{
    answer = this->giveStructuralCrossSection()->giveGeneralizedStress_Beam3d( strain, gp, tStep );
}


void RBSBeam3d::FiberedCrossSectionInterface_computeStrainVectorInFiber(
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


void RBSBeam3d ::computeKappaCoeffs( TimeStep *tStep )
{
    // kappa_y = (6*E*Iy)/(k*G*A*l^2)

    FloatMatrix d;
    double l = this->computeLength();

    this->computeConstitutiveMatrixAt( d, ElasticStiffness, integrationRulesArray[0]->getIntegrationPoint( 0 ), tStep );

    //  kappay = 6. * d.at(5, 5) / ( d.at(3, 3) * l * l );
    //  kappaz = 6. * d.at(6, 6) / ( d.at(2, 2) * l * l );
    if ( d.at( 3, 3 ) != 0. ) {
        kappay = 6. * d.at( 5, 5 ) / ( d.at( 3, 3 ) * l * l );
    } else {
        kappay = 0.;
    }
    if ( d.at( 2, 2 ) != 0. ) {
        kappaz = 6. * d.at( 6, 6 ) / ( d.at( 2, 2 ) * l * l );
    } else {
        kappaz = 0.;
    }
}

double RBSBeam3d ::giveKappayCoeff( TimeStep *tStep )
{
    if ( kappay < 0.0 ) {
        this->computeKappaCoeffs( tStep );
    }

    return kappay;
}

double RBSBeam3d ::giveKappazCoeff( TimeStep *tStep )
{
    if ( kappaz < 0.0 ) {
        this->computeKappaCoeffs( tStep );
    }

    return kappaz;
}


void RBSBeam3d::giveInternalForcesVectorAtPoint( FloatArray &answer, TimeStep *tStep, FloatArray &coords )
{
    // computes exact global end-forces vector
    FloatArray loadEndForces, iF;
    IntArray leftIndx = {
        1, 2, 3, 4, 5, 6
    };
    this->giveEndForcesVector( iF, tStep );

    answer.beSubArrayOf( iF, leftIndx );
    Node *nodeA;

    nodeA     = this->giveNode( 1 );
    double dx = nodeA->giveCoordinate( 1 ) - coords.at( 1 );
    double dy = nodeA->giveCoordinate( 2 ) - coords.at( 2 );
    double dz = nodeA->giveCoordinate( 3 ) - coords.at( 3 );
    double ds = sqrt( dx * dx + dy * dy + dz * dz );

    answer.at( 5 ) += iF.at( 3 ) * ds;
    answer.at( 6 ) -= iF.at( 2 ) * ds;


    // loop over body load array first
    int nBodyLoads = this->giveBodyLoadArray()->giveSize();
    FloatArray help;

    for ( int i = 1; i <= nBodyLoads; i++ ) {
        int id           = bodyLoadArray.at( i );
        Load *load       = domain->giveLoad( id );
        bcGeomType ltype = load->giveBCGeoType();
        if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ForceLoadBVT ) ) {
            this->computeInternalForcesFromBodyLoadVectorAtPoint( help, load, tStep, VM_Total, coords, ds ); // this one is local
            answer.add( help );
        } else {
            if ( load->giveBCValType() != TemperatureBVT && load->giveBCValType() != EigenstrainBVT ) {
                // temperature and eigenstrain is handled separately at computeLoadVectorAt subroutine
                OOFEM_ERROR( "body load %d is of unsupported type (%d)", id, ltype );
            }
        }
    }

    // loop over boundary load array
    int nBoundaryLoads = this->giveBoundaryLoadArray()->giveSize() / 2;
    for ( int i = 1; i <= nBoundaryLoads; i++ ) {
        int n      = boundaryLoadArray.at( 1 + ( i - 1 ) * 2 );
        int id     = boundaryLoadArray.at( i * 2 );
        Load *load = domain->giveLoad( n );
        BoundaryLoad *bLoad;
        if ( ( bLoad = dynamic_cast<BoundaryLoad *>( load ) ) ) {
            bcGeomType ltype = load->giveBCGeoType();
            if ( ltype == EdgeLoadBGT ) {
                this->computeInternalForcesFromBoundaryEdgeLoadVectorAtPoint( help, bLoad, id,
                    ExternalForcesVector, VM_Total, tStep, coords, ds, false );
                answer.add( help );
            } else {
                OOFEM_ERROR( "boundary load %d is of unsupported type (%d)", id, ltype );
            }
        }
    }
    // add exact end forces due to non-nodal loading
    //    this->computeForceLoadVectorAt(loadEndForces, tStep, VM_Total, coords); // will compute only contribution of loads applied directly on receiver (not using sets)
    if ( loadEndForces.giveSize() ) {
        answer.subtract( loadEndForces );
    }

    // add exact end forces due to non-nodal loading applied indirectly (via sets)
    BCTracker *bct                   = this->domain->giveBCTracker();
    BCTracker ::entryListType bcList = bct->getElementRecords( this->number );

    for ( BCTracker ::entryListType ::iterator it = bcList.begin(); it != bcList.end(); ++it ) {
        GeneralBoundaryCondition *bc = this->domain->giveBc( ( *it ).bcNumber );
        BodyLoad *bodyLoad;
        BoundaryLoad *boundaryLoad;
        if ( bc->isImposed( tStep ) ) {
            if ( ( bodyLoad = dynamic_cast<BodyLoad *>( bc ) ) ) { // body load
                this->computeInternalForcesFromBodyLoadVectorAtPoint( help, bodyLoad, tStep, VM_Total, coords, ds ); // this one is local
                //answer.subtract(help);
            } else if ( ( boundaryLoad = dynamic_cast<BoundaryLoad *>( bc ) ) ) {
                // compute Boundary Edge load vector in GLOBAL CS !!!!!!!
                this->computeInternalForcesFromBoundaryEdgeLoadVectorAtPoint( help, boundaryLoad, ( *it ).boundaryId,
                    ExternalForcesVector, VM_Total, tStep, coords, ds, false );
            }
            answer.add( help );
        }
    }


    if ( subsoilMat ) {
        // @todo: linear subsoil assumed here; more general approach should integrate internal forces
        FloatMatrix k;
        FloatArray u, F;
        this->computeSubSoilStiffnessMatrix( k, TangentStiffness, tStep );
        this->computeVectorOf( VM_Total, tStep, u );
        F.beProductOf( k, u );
        answer.add( F );
    }


    answer.times( -1 );
}


void RBSBeam3d::computeInternalForcesFromBoundaryEdgeLoadVectorAtPoint( FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode, TimeStep *tStep, FloatArray &pointCoords, double ds, bool global )
{
    answer.clear();

    if ( edge != 1 ) {
        OOFEM_ERROR( "Beam3D only has 1 edge (the midline) that supports loads. Attempted to apply load to edge %d", edge );
    }

    if ( type != ExternalForcesVector ) {
        return;
    }

    FloatArray coords, t;
    FloatMatrix T;


    for ( GaussPoint *gp : *this->giveDefaultIntegrationRulePtr() ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        this->computeGlobalCoordinates( coords, lcoords, pointCoords );
        if ( load->giveFormulationType() == Load ::FT_Entity ) {
            load->computeValues( t, tStep, lcoords, { D_u, D_v, D_w, R_u, R_v, R_w }, mode );
        } else {
            load->computeValues( t, tStep, coords, { D_u, D_v, D_w, R_u, R_v, R_w }, mode );
        }

        if ( load->giveCoordSystMode() == Load ::CST_Global ) {
            if ( this->computeLoadGToLRotationMtrx( T ) ) {
                t.rotatedWith( T, 'n' );
            }
        }


        double dl = gp->giveWeight() * 0.5 * ds;
        FloatArray f;
        f = t;
        f.at( 5 ) += f.at( 3 ) * ( lcoords.at( 1 ) + 1 ) * ds / 2;
        f.at( 6 ) -= f.at( 2 ) * ( lcoords.at( 1 ) + 1 ) * ds / 2;
        answer.add( dl, f );
    }

    if ( global ) {
        // Loads from sets expects global c.s.
        this->computeGtoLRotationMatrix( T );
        answer.rotatedWith( T, 't' );
    }
}

void RBSBeam3d::computeInternalForcesFromBodyLoadVectorAtPoint( FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode, FloatArray &pointCoords, double ds )
// Computes numerically the load vector of the receiver due to the body
// loads, at tStep.
// load is assumed to be in global cs.
// load vector is then transformed to coordinate system in each node.
// (should be global coordinate system, but there may be defined
//  different coordinate system in each node)
{
    double dens, dV;
    FloatArray force, ntf;
    FloatMatrix n, T;
    FloatArray lc( 1 );

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
        OOFEM_ERROR( "unknown load type" );
    }

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt( force, tStep, mode );
    force.times( this->giveCrossSection()->give( CS_Area, lc, this ) );
    // transform from global to element local c.s
    if ( this->computeLoadGToLRotationMtrx( T ) ) {
        force.rotatedWith( T, 'n' );
    }

    answer.clear();

    if ( force.giveSize() ) {
        for ( GaussPoint *gp : *this->giveDefaultIntegrationRulePtr() ) {
            const FloatArray &lcoords = gp->giveNaturalCoordinates();
            this->computeNmatrixAt( gp->giveSubPatchCoordinates(), n );
            dV   = gp->giveWeight() * 0.5 * ds;
            dens = this->giveCrossSection()->give( 'd', gp );
            FloatArray iF;
            iF = force;
            iF.at( 5 ) += force.at( 3 ) * ( lcoords.at( 1 ) + 1 ) * ds / 2;
            iF.at( 6 ) -= force.at( 2 ) * ( lcoords.at( 1 ) + 1 ) * ds / 2;
            answer.add( dV * dens, iF );
        }
    } else {
        return;
    }
}


void RBSBeam3d::setNumberOfGaussPoints( int nip )
{
    numberOfGaussPoints = nip;
}

} // end namespace oofem
// <temporary/>
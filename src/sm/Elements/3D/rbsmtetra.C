/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
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

#include "sm/Elements/3D/rbsmtetra.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "fei3dtetlin.h"
#include "classfactory.h"

#include "rigidarmnode.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
 #include "sm/Materials/rcm2.h"

 #include <Etetrawd.h>
#endif

namespace oofem {
REGISTER_Element(RBSMTetra);

// S: todo: update or confirm
FEI3dTetLin RBSMTetra :: interpolation;

// S: todo: update or confirm
RBSMTetra :: RBSMTetra(int n, Domain *aDomain) :
    Structural3DElement(n, aDomain)
    /*,
        ZZNodalRecoveryModelInterface(this),
        NodalAveragingRecoveryModelInterface(),
        SPRNodalRecoveryModelInterface(),
        SpatialLocalizerInterface(this),
        ZZErrorEstimatorInterface(this),
        HuertaErrorEstimatorInterface()
    */
{
    numberOfDofMans  = 4;
    numberOfGaussPoints = 1;
}


void RBSMTetra::initializeFrom( InputRecord &ir )
{
    numberOfGaussPoints = 1;
    numberOfCornerNodes = 4;

    // A safer option to initialize the RBSM Tetra is to
    // make a dummy input record after cloning corner nodes.
    // currently we first initialize the element then clone
    // the corner nodes.
    Structural3DElement::initializeFrom( ir );

    RBSMTetra::setCornerNodesFromIr( ir );
    RBSMTetra::makeDofmanagers( ir );

}


void RBSMTetra::setCornerNodesFromIr( InputRecord &ir )
// Obtains the coordinates of rigid body cell
{

    RBSMTetra::cornerNodes.resize(numberOfCornerNodes );

    // Obtain the node IDs of the rigid body cell
    IR_GIVE_FIELD(ir, cornerNodes, _IFT_Element_nodes);
}

void RBSMTetra::makeDofmanagers( InputRecord &ir )
// Makes central and cloned corner nodes for rigid body
{
    int nextID;
    std::vector<OOFEMTXTInputRecord> rbsmInputRecords;

    // find next available DOF manager's global number
    nextID = RBSMTetra::nextDofmanagerID();
    // create fake input records to make RBSM DOF managers
    RBSMTetra::rbsmDummyIr(ir, rbsmInputRecords, nextID);

    // create central DOF manager (master)
    this->centerDofmanager =
        RBSMTetra::makeDofmanager( rbsmInputRecords.at( 0 ), nextID );

    // create cloned corner nodes (rigid arm)
    // (make sure that the dofManArray is initialized)
    for ( int i = 1; i <= numberOfCornerNodes; ++i ) {
        this->dofManArray.at( i ) =
            RBSMTetra::makeDofmanager( rbsmInputRecords.at( i ) );
    }
}


int RBSMTetra::makeDofmanager( InputRecord &dummyIr )
// Makes central DOF manager and returns given global number
{
    int id;
    id = RBSMTetra::nextDofmanagerID();

    return RBSMTetra::makeDofmanager( dummyIr, id );
}

int RBSMTetra::makeDofmanager( InputRecord &dummyIr, int id )
// Makes central DOF manager and returns given global number
{
    int num, nDofman = 0;
    Domain *d = this->giveDomain();
    std::vector<OOFEMTXTInputRecord> rbsmInputRecords;

    nDofman = d->dofManagerList.size();
    if ( nDofman <= 0 ) {
        OOFEM_ERROR("Domain returned invalid DOF managers count: %d\n", id);
    }

    // number for central DOF manager of the rigid body
    num = nDofman + 1;

    std::unique_ptr<DofManager> dmanCenter( classFactory.createDofManager( _IFT_Node_Name, num, d ) );
    if ( !dmanCenter ) {
        OOFEM_ERROR("Couldn't create node of type: %s\n", _IFT_Node_Name);
    }

    dmanCenter->initializeFrom(dummyIr);

    // make sure that global number is unique
    auto hasSameNum{
        [&id]( std::unique_ptr<oofem::DofManager> &dman ) { return dman->giveGlobalNumber() == id; }
    };
    auto it = std::find_if( d->dofManagerList.begin(), d->dofManagerList.end(), hasSameNum);
    if (it!=d->dofManagerList.end()){
        // the global ID already exists
        OOFEM_ERROR( "Failed to create DOF manager; the global number '%d' is already taken", id );
    }

    dmanCenter->setGlobalNumber( id );

    d->resizeDofManagers( nDofman + 1 );
    //d->dofManagerList[nDofman] = std::move( dmanCenter );
    d->setDofManager(nDofman + 1, std::move(dmanCenter));

    // I think dofManLabelMap should be a property of domain, in that case
    // we could update mapping from label to local numbers for dofmans.
    // Since this is not the case, we override the Element::updateLocalNumbering
    /*
    if ( d->dofManLabelMap.find(id) == d->dofManLabelMap.end() ) {
        // label does not exist yet
        d->dofManLabelMap[num] = nDofman + 1;
    } else {
        OOFEM_ERROR("iDofmanager entry already exist (label=%d)", num);
    }
    */
    return id;
}

void RBSMTetra::updateLocalNumbering(EntityRenumberingFunctor &f)
{
    // How element class do this:
    /*
    for ( auto &dnum : dofManArray ) {
        //Functor:
        std :: map< int, int > :: const_iterator it;
        if ( ( it = dofmanMap.find(dnum) ) != dofmanMap.end() ) {
            dnum it->second;
        }
        OOFEM_ERROR("component label %d not found", dnum);
    }
    */

    // Temporary fix:
    // Replicate EntityRenumberingFunctor behaviour using
    // domain->mDofManPlaceInArray instead of dofmanMap
    for ( auto &dnum : dofManArray ) {
        std :: map< int, int > :: const_iterator it;
        dnum = this->giveDomain()->giveDofManPlaceInArray(dnum);
    }
    return;
}


int RBSMTetra::nextDofmanagerID()
// finds the next available DOF manager's global number
{
    int num, count, nDofman = 0;
    Domain *d = this->giveDomain();

    nDofman = d->dofManagerList.size();
    if ( nDofman <= 0 ) {
        OOFEM_ERROR("Domain returned invalid DOF managers count: %d\n", num);
    }

    // first try for the next DOF manager's global number
    num = nDofman + 1;

    // make sure that global number is unique
    count = 0;
    auto hasSameNum{
        //[num]( std :: unique_ptr< DofManager > dman ) { return dman->giveGlobalNumber() == num; }
        [&num]( std::unique_ptr<oofem::DofManager> &dman ) { return dman->giveGlobalNumber() == num; }
    };
    auto it = std::find_if( d->dofManagerList.begin(), d->dofManagerList.end(), hasSameNum);
    while (it!=d->dofManagerList.end()){
        // try to resolve duplicated number
        num++;
        count++;
        it = std::find_if(d->dofManagerList.begin(), d->dofManagerList.end(), hasSameNum);
        if (count > nDofman ){
            // prevent infinite loop
            OOFEM_ERROR( "Failed to find a unique node id for 'RBSM Tetra element'" );
        }
    }

    return num;
}


std::vector<FloatArray> RBSMTetra::coordsFromIr( InputRecord &ir )
// Obtains the coordinates of rigid body cell
{
    IntArray cornerNodes;
    std::vector<FloatArray> nodeCoords;

    nodeCoords.resize(numberOfCornerNodes+1);
    nodeCoords.at(0).zero();

    // Obtain the node IDs of the rigid body cell
    IR_GIVE_FIELD(ir, cornerNodes, _IFT_Element_nodes);
    // Obtain node coordinates
    for ( int i = 0; i < numberOfCornerNodes; ++i ) {
        nodeCoords.at(i+1) =
            this->giveDomain()->giveDofManager(cornerNodes(i))->giveCoordinates();
        nodeCoords.at(0).add( nodeCoords.at(i+1));
    }
    // center node
    nodeCoords.at( 0 ).times( 1 / numberOfCornerNodes );

    return nodeCoords;
}


void RBSMTetra::rbsmDummyIr(
    InputRecord &irIn, std::vector<OOFEMTXTInputRecord> &irOut, int master )
// Obtains the coordinates of rigid body cell
{
    char buff[256];
    std::vector<FloatArray> nodeCoords;
    std::string irString;
    std::vector<std::string> dummyInputStrings;

    nodeCoords = RBSMTetra::coordsFromIr( irIn );
    irOut.resize( numberOfCornerNodes + 1 );
    dummyInputStrings.resize( numberOfCornerNodes + 1 );

    // central DOF man
    sprintf( buff, "%s %i   %s 3  %f %f %f",
        //  "Node"  num  "Coords"   x   y   z
        _IFT_Node_Name, 0, _IFT_Node_coords,
        nodeCoords.at( 0 ).at( 1 ),
        nodeCoords.at( 0 ).at( 2 ),
        nodeCoords.at( 0 ).at( 3 ) );
    irString = buff;
    // central DOF manager dummy input record
    irOut.at( 0 ) = OOFEMTXTInputRecord( 0, irString );

    // corner (rigid arm) DOF man
    for ( int i = 1; i < numberOfCornerNodes + 1; ++i ) {
        sprintf( buff, "%s %i   %s 3  %f %f %f   %s %i",
            //  "RigidArmNode"  num  "Coords"   x   y   z "Master" id
            _IFT_RigidArmNode_Name, 0, _IFT_Node_coords,
            nodeCoords.at( i ).at( 1 ),
            nodeCoords.at( i ).at( 2 ),
            nodeCoords.at( i ).at( 3 ),
            _IFT_RigidArmNode_master, master );
        irString = buff;
        // rigid arm dummy input records
        irOut.at( i ) = OOFEMTXTInputRecord( 0, irString );
    }
}


// S: todo: update or confirm
Interface *RBSMTetra::giveInterface( InterfaceType interface )
{
    /*
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >(this);
    } else if ( interface == ZZErrorEstimatorInterfaceType ) {
        return static_cast< ZZErrorEstimatorInterface * >(this);
    } else if ( interface == HuertaErrorEstimatorInterfaceType ) {
        return static_cast< HuertaErrorEstimatorInterface * >(this);
    }
    */

    return NULL;
}

// S: todo: update or confirm
FEInterpolation *RBSMTetra::giveInterpolation() const
{
    return & interpolation;
}


// S: todo: update or confirm
void RBSMTetra ::computeLumpedMassMatrix( FloatMatrix &answer, TimeStep *tStep )
// Returns the lumped mass matrix of the receiver.
{
    GaussPoint *gp;
    double dV, mss1;

    answer.resize(12, 12);
    answer.zero();
    gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    dV = this->computeVolumeAround(gp);
    double density = this->giveStructuralCrossSection()->give('d', gp);
    mss1 = dV * density / 4.;

    for ( int i = 1; i <= 12; i++ ) {
        answer.at(i, i) = mss1;
    }
}


// S: todo: update or confirm
/*
void RBSMTetra ::NodalAveragingRecoveryMI_computeNodalValue(
    FloatArray &answer, int node,
    InternalStateType type, TimeStep *tStep )
{
    GaussPoint *gp;

    if ( numberOfGaussPoints != 1 ) {
        answer.clear(); // for more gp's need to be refined
        return;
    }

    gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    giveIPValue(answer, gp, type, tStep);
}
*/


// S: todo: update or confirm
/*
void RBSMTetra::SPRNodalRecoveryMI_giveSPRAssemblyPoints( IntArray &pap )
{
    pap.resize( 4 );
    pap.at( 1 ) = this->giveNode( 1 )->giveNumber();
    pap.at( 2 ) = this->giveNode( 2 )->giveNumber();
    pap.at( 3 ) = this->giveNode( 3 )->giveNumber();
    pap.at( 4 ) = this->giveNode( 4 )->giveNumber();
}
*/


// S: todo: update or confirm
/*
void RBSMTetra ::SPRNodalRecoveryMI_giveDofMansDeterminedByPatch( IntArray &answer, int pap )
{
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i <= 4; i++ ) {
        if ( this->giveNode(i)->giveNumber() == pap ) {
            found = 1;
        }
    }

    if ( found ) {
        answer.at(1) = pap;
    } else {
        OOFEM_ERROR("node unknown");
    }
}
*/


// S: todo: update or confirm
/*
int RBSMTetra ::SPRNodalRecoveryMI_giveNumberOfIP()
{
    return this->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints();
}
*/


// S: todo: update or confirm
/*
SPRPatchType RBSMTetra::SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_3dBiLin;
}
*/


// S: todo: update or confirm
/*
void RBSMTetra ::HuertaErrorEstimatorI_setupRefinedElementProblem(
    RefinedElement *refinedElement, int level, int nodeId,
    IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
    HuertaErrorEstimatorInterface ::SetupMode sMode, TimeStep *tStep,
    int &localNodeId, int &localElemId, int &localBcId,
    IntArray &controlNode, IntArray &controlDof,
    HuertaErrorEstimator ::AnalysisMode aMode )
{
    double x = 0.0, y = 0.0, z = 0.0;
    int nodes = 4, sides = 6, faces = 4;

    static int sideNode[6][2] = { { 1, 2 }, { 2, 3 }, { 3, 1 }, { 1, 4 }, { 2, 4 }, { 3, 4 } };
    static int faceNode[4][3] = { { 1, 2, 3 }, { 1, 2, 4 }, { 2, 3, 4 }, { 3, 1, 4 } };

    */
/* ordering of hexa nodes must be compatible with refinedElement connectivity ordering;
     * generally the ordering is: corner side side face side face face center;
     * however the concrete ordering is element dependent (see refineMeshGlobally source if in doubts) *//*


    int hexaSideNode[4][3] = { { 1, 3, 4 }, { 2, 1, 5 }, { 3, 2, 6 }, { 4, 6, 5 } };
    int hexaFaceNode[4][3] = { { 1, 2, 4 }, { 1, 3, 2 }, { 1, 4, 3 }, { 4, 2, 3 } };

    FloatArray corner[4], midSide[6], midFace[4], midNode;
    if ( sMode == HuertaErrorEstimatorInterface ::NodeMode || ( sMode == HuertaErrorEstimatorInterface ::BCMode && aMode == HuertaErrorEstimator ::HEE_linear ) ) {
        for ( int inode = 0; inode < nodes; inode++ ) {
            corner[inode] = this->giveNode( inode + 1 )->giveCoordinates();

            x += corner[inode].at( 1 );
            y += corner[inode].at( 2 );
            z += corner[inode].at( 3 );
        }

        for ( int iside = 0; iside < sides; iside++ ) {
            midSide[iside].resize( 3 );

            int nd1 = sideNode[iside][0] - 1;
            int nd2 = sideNode[iside][1] - 1;

            midSide[iside].at( 1 ) = ( corner[nd1].at( 1 ) + corner[nd2].at( 1 ) ) / 2.0;
            midSide[iside].at( 2 ) = ( corner[nd1].at( 2 ) + corner[nd2].at( 2 ) ) / 2.0;
            midSide[iside].at( 3 ) = ( corner[nd1].at( 3 ) + corner[nd2].at( 3 ) ) / 2.0;
        }

        midNode.resize( 3 );

        midNode.at( 1 ) = x / nodes;
        midNode.at( 2 ) = y / nodes;
        midNode.at( 3 ) = z / nodes;

        for ( int iface = 0; iface < faces; iface++ ) {
            x = y = z = 0.0;
            for ( int inode = 0; inode < 3; inode++ ) {
                int nd = faceNode[iface][inode] - 1;
                x += corner[nd].at( 1 );
                y += corner[nd].at( 2 );
                z += corner[nd].at( 3 );
            }

            midFace[iface].resize( 3 );

            midFace[iface].at( 1 ) = x / 3;
            midFace[iface].at( 2 ) = y / 3;
            midFace[iface].at( 3 ) = z / 3;
        }
    }

    this->setupRefinedElementProblem3D( this, refinedElement, level, nodeId, localNodeIdArray, globalNodeIdArray,
        sMode, tStep, nodes, corner, midSide, midFace, midNode,
        localNodeId, localElemId, localBcId, hexaSideNode, hexaFaceNode,
        controlNode, controlDof, aMode, "LSpace" );
}
*/

// S: todo: update or confirm
/*
void RBSMTetra :: HuertaErrorEstimatorI_computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    computeNmatrixAt(gp->giveSubPatchCoordinates(), answer);
}
*/

/*
#ifdef __OOFEG
 #define TR_LENGHT_REDUCT 0.3333

void RBSMTetra :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    WCRec p [ 4 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    EASValsSetFillStyle(FILL_SOLID);
#ifdef __PARALLEL_MODE
    if (this->giveParallelMode() == Element_remote) {
      EASValsSetColor( gc.getRemoteElementColor() );
      EASValsSetEdgeColor( gc.getRemoteElementEdgeColor() );
    }
#endif
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3);
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
    p [ 2 ].z = ( FPNum ) this->giveNode(3)->giveCoordinate(3);
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveCoordinate(1);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveCoordinate(2);
    p [ 3 ].z = ( FPNum ) this->giveNode(4)->giveCoordinate(3);

    go =  CreateTetra(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void RBSMTetra :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    WCRec p [ 4 ];
    GraphicObj *go;
    double defScale = gc.getDefScale();

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    EASValsSetFillStyle(FILL_SOLID);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, defScale);
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, defScale);
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 2 ].z = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(3, tStep, defScale);
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 3 ].z = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(3, tStep, defScale);

    go =  CreateTetra(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


void RBSMTetra :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int i, indx, result = 0;
    WCRec p [ 4 ];
    GraphicObj *tr;
    FloatArray v [ 4 ];
    double s [ 4 ], defScale = 0.0;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( gc.giveIntVarMode() == ISM_recovered ) {
        for ( i = 1; i <= 4; i++ ) {
            result += this->giveInternalStateAtNode(v [ i - 1 ], gc.giveIntVarType(), gc.giveIntVarMode(), i, tStep);
        }

        if ( result != 4 ) {
            return;
        }
    } else if ( gc.giveIntVarMode() == ISM_local ) {
      if ( integrationRulesArray [ 0 ]->giveNumberOfIntegrationPoints() != 1 ) return;
      GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
      if ( giveIPValue(v [ 0 ], gp, gc.giveIntVarType(), tStep) == 0 ) {
        return;
      }
      v[1]=v[0];
      v[2]=v[0];
      v[3]=v[0];
    }

    indx = gc.giveIntVarIndx();

    for ( i = 1; i <= 4; i++ ) {
        s [ i - 1 ] = v [ i - 1 ].at(indx);
    }

    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    if ( gc.getScalarAlgo() == SA_ISO_SURF ) {
        for ( i = 0; i < 4; i++ ) {
            if ( gc.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(3, tStep, defScale);
            } else {
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(3);
            }
        }

        gc.updateFringeTableMinMax(s, 4);
        tr = CreateTetraWD(p, s);
        EGWithMaskChangeAttributes(LAYER_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    }
}

void
RBSMTetra :: drawSpecial(oofegGraphicContext &gc, TimeStep *tStep)
{
    int i, j, k;
    WCRec q [ 4 ];
    GraphicObj *tr;
    double defScale = gc.getDefScale();
    FloatArray crackStatuses, cf;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( gc.giveIntVarType() == IST_CrackState ) {
        int crackStatus;
        double xc, yc, zc, length;
        FloatArray crackDir;

        if ( numberOfGaussPoints != 1 ) {
            return;
        }

        //   for (GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        {
            GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
            if ( this->giveIPValue(cf, gp, IST_CrackedFlag, tStep) == 0 ) {
                return;
            }

            if ( ( int ) cf.at(1) == 0 ) {
                return;
            }

            //
            // obtain gp global coordinates - here only one exists
            // it is in centre of gravity.
            xc = yc = zc = 0.;
            for ( i = 0; i < 4; i++ ) {
                if ( gc.getInternalVarsDefGeoFlag() ) {
                    // use deformed geometry
                    xc += ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                    yc += ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                    zc += ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(3, tStep, defScale);
                } else {
                    xc += ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                    yc += ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                    zc += ( FPNum ) this->giveNode(i + 1)->giveCoordinate(3);
                }
            }

            xc = xc / 4.;
            yc = yc / 4.;
            zc = zc / 4.;
            length = TR_LENGHT_REDUCT * pow(this->computeVolumeAround(gp), 1. / 3.) / 2.0;
            if ( this->giveIPValue(crackDir, gp, IST_CrackDirs, tStep) ) {
                this->giveIPValue(crackStatuses, gp, IST_CrackStatuses, tStep);


                for ( i = 1; i <= 3; i++ ) {
                    crackStatus = ( int ) crackStatuses.at(i);
                    if ( ( crackStatus != pscm_NONE ) && ( crackStatus != pscm_CLOSED ) ) {
                        // draw a crack
                        // this element is 3d element

                        if ( i == 1 ) {
                            j = 2;
                            k = 3;
                        } else if ( i == 2 ) {
                            j = 3;
                            k = 1;
                        } else {
                            j = 1;
                            k = 2;
                        }

                        q [ 0 ].x = ( FPNum ) xc + 0.5 * crackDir.at(0 + j) * length + 0.5 * crackDir.at(0 + k) * length;
                        q [ 0 ].y = ( FPNum ) yc + 0.5 * crackDir.at(3 + j) * length + 0.5 * crackDir.at(3 + k) * length;
                        q [ 0 ].z = ( FPNum ) zc + 0.5 * crackDir.at(6 + j) * length + 0.5 * crackDir.at(6 + k) * length;
                        q [ 1 ].x = ( FPNum ) xc + 0.5 * crackDir.at(0 + j) * length - 0.5 * crackDir.at(0 + k) * length;
                        q [ 1 ].y = ( FPNum ) yc + 0.5 * crackDir.at(3 + j) * length - 0.5 * crackDir.at(3 + k) * length;
                        q [ 1 ].z = ( FPNum ) zc + 0.5 * crackDir.at(6 + j) * length - 0.5 * crackDir.at(6 + k) * length;
                        q [ 2 ].x = ( FPNum ) xc - 0.5 * crackDir.at(0 + j) * length - 0.5 * crackDir.at(0 + k) * length;
                        q [ 2 ].y = ( FPNum ) yc - 0.5 * crackDir.at(3 + j) * length - 0.5 * crackDir.at(3 + k) * length;
                        q [ 2 ].z = ( FPNum ) zc - 0.5 * crackDir.at(6 + j) * length - 0.5 * crackDir.at(6 + k) * length;
                        q [ 3 ].x = ( FPNum ) xc - 0.5 * crackDir.at(0 + j) * length + 0.5 * crackDir.at(0 + k) * length;
                        q [ 3 ].y = ( FPNum ) yc - 0.5 * crackDir.at(3 + j) * length + 0.5 * crackDir.at(3 + k) * length;
                        q [ 3 ].z = ( FPNum ) zc - 0.5 * crackDir.at(6 + j) * length + 0.5 * crackDir.at(6 + k) * length;

                        EASValsSetLayer(OOFEG_CRACK_PATTERN_LAYER);
                        EASValsSetLineWidth(OOFEG_CRACK_PATTERN_WIDTH);
                        if ( ( crackStatus == pscm_SOFTENING ) || ( crackStatus == pscm_OPEN ) ) {
                            EASValsSetColor( gc.getActiveCrackColor() );
                        } else {
                            EASValsSetColor( gc.getCrackPatternColor() );
                        }

                        //      EASValsSetFillStyle (FILL_HOLLOW);
                        tr = CreateQuad3D(q);
                        EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
                        EMAddGraphicsToModel(ESIModel(), tr);
                    }
                }
            }
        }
    }
}
#endif
*/

} // end namespace oofem

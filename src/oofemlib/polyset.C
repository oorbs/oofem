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

#include <Elements/3D/rbsmtetra.h>
#include "inputrecord.h"
#include "dynamicinputrecord.h"
#include "domain.h"
#include "classfactory.h"
#include "polyset.h"


/*
#include "set.h"
#include "error.h"
#include "intarray.h"
#include "element.h"
#include "node.h"
#include "range.h"
#include "feinterpol.h"
#include "feinterpol2d.h"
#include "feinterpol3d.h"
#include "contextioerr.h"
#include <list>
*/

namespace oofem {
REGISTER_Set(Polyset);

void Polyset::initializeFrom(InputRecord &ir)
{
    FEMComponent::initializeFrom(ir);
    Set::initializeFrom(ir);

    IntArray inputNodes;
    if ( ir.hasField(_IFT_Polyset_rtetClonesOfGeoNodes) ) {
        IR_GIVE_FIELD( ir, inputNodes, _IFT_Polyset_rtetClonesOfGeoNodes );
        this->clonesOfGeoNodes( this->nodes, inputNodes );
    }
    if ( ir.hasField(_IFT_Polyset_rtetCellsOfGeoNodes) ) {
        IR_GIVE_FIELD( ir, inputNodes, _IFT_Polyset_rtetCellsOfGeoNodes );
        this->cellNodesOfGeoNodes( this->nodes, inputNodes );
    }

    // print postprocess hint
    OOFEM_LOG_INFO( "\n** USE THE FOLLOWING INFORMATION TO POSTPROCESS %s NUM %d:",
        this->giveClassName(), this->giveNumber() );
    for ( int i = 0; i < this->nodes.giveSize(); ++i ) {
        OOFEM_LOG_INFO( "\n#REACTION number %d dof 3", this->nodes[i] );
    }
}


void Polyset :: giveInputRecord(DynamicInputRecord &input)
{
    input.setRecordKeywordField( _IFT_Polyset_Name, this->giveNumber() );

    if ( this->giveNodeList().giveSize() ) {
        input.setField(this->nodes, _IFT_Polyset_rtetClonesOfGeoNodes);
    }
    if ( this->giveNodeList().giveSize() ) {
        input.setField(this->nodes, _IFT_Polyset_rtetCellsOfGeoNodes);
    }

/*
    // for future development
//    if ( this->giveElementList().giveSize() ) {
//        input.setField(this->elements, _IFT_Set_elements);
//    }
//    if ( this->giveBoundaryList().giveSize() ) {
//        input.setField(this->elementBoundaries, _IFT_Set_elementBoundaries);
//    }
//    if ( this->giveEdgeList().giveSize() ) {
//        input.setField(this->elementEdges, _IFT_Set_elementEdges);
//    }
//    if ( this->giveSurfaceList().giveSize() ) {
//        input.setField(this->elementSurfaces, _IFT_Set_elementSurfaces);
//    }
*/
}

void Polyset::clonesOfGeoNodes( IntArray &answer, const IntArray &nodes )
{
    // Find the cell nodes;
    for ( auto node : nodes ) {
        std::set<int> clones = RBSMTetra::giveClonesOfGeoNode( node );
        for ( int clone : clones ) {
            answer.insertSortedOnce( clone );
        }
    }
}

void Polyset::cellNodesOfGeoNodes( IntArray &answer, const IntArray &nodes )
{
    // Find the cell nodes;
    for ( auto node : nodes ) {
        std::set<int> cellsNodesNumber = this->cellNodesOfGeoNode( node );
        for ( int globalNumber : cellsNodesNumber ) {
            answer.insertSortedOnce( globalNumber );
        }
    }
}

std::set<int> Polyset::cellNodesOfGeoNode( int geoNode )
{
    std::set<int> answer;
    std::set<int> cellElements = RBSMTetra::giveCellElementsOfGeoNode( geoNode );

    for ( int index : cellElements ) {
        RBSMTetra *elem =
            dynamic_cast<RBSMTetra *>( domain->elementList[index - 1].get() );
        if ( elem ) {
            answer.insert( elem->giveCellDofmanagerNumber() );
        }
    }
    return answer;
}


}

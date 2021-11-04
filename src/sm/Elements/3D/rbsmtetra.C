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
#include "sm/Elements/Beams/beam3d.h"
#include "sm/Elements/Beams/rbsbeam3d.h"
#include "sm/Elements/Beams/rbsmbeam3d.h"

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#include "oofegutils.h"
#include "sm/Materials/rcm2.h"

#include <Etetrawd.h>
#endif

namespace oofem {
REGISTER_Element( RBSMTetra );

// use Mindlin beam (RBSBeam3d=>RBSMBeam3D)
#define MINDLIN

// S: todo: update or confirm
FEI3dTetLin RBSMTetra::interpolation;
std::map<int, std::set<int>> RBSMTetra::clonesOfGeoNode;
std::map<int, std::set<int>> RBSMTetra::cellElementsOfGeoNode;
std::map<std::vector<int>, std::set<int>> RBSMTetra::mapFacetElement;

// S: todo: update or confirm
RBSMTetra::RBSMTetra( int n, Domain *aDomain ) :
    Structural3DElement( n, aDomain )
/*,
        ZZNodalRecoveryModelInterface(this),
        NodalAveragingRecoveryModelInterface(),
        SPRNodalRecoveryModelInterface(),
        SpatialLocalizerInterface(this),
        ZZErrorEstimatorInterface(this),
        HuertaErrorEstimatorInterface()
    */
{
    numberOfDofMans     = 4;
    numberOfGaussPoints = 1;
}


void RBSMTetra::initializeFrom( InputRecord &ir )
{
    numberOfGaussPoints = 1;
    numberOfCornerNodes = 4;
    int numberOfFacets  = 4;

    // currently we first initialize the element then clone
    // the corner nodes.
    Structural3DElement::initializeFrom( ir );

    RBSMTetra::setGeoNodesFromIr( ir );
    RBSMTetra::makeDofmanagers( ir );

    RBSMTetra::updateClonesOfGeoNode();
    RBSMTetra::updateCellElementsOfGeoNode();

    std::map<int, IntArray> existingSisters =
        RBSMTetra::updateCellElementsOfFacets();

    // make springs beams
    //@todo: move to post-initialization, remove override RBSMTetra::setCrossSection
    springsBeams.resize( numberOfFacets );
    int number, startPoint, endPoint;
    endPoint = centerDofmanager;
    for ( const auto &facetExistingSisters : existingSisters ) {
        int count = 0;
        for ( int neesan : facetExistingSisters.second ) {
            count++;
            if ( count > 1 ) {
                // currently only 1 sister element is accepted for each facet
                OOFEM_ERROR( "There are more than 1 adjacent element for elem %d",
                    this->giveGlobalNumber() );
            }
            if ( RBSMTetra *e = //> make sure sister element is an RBSM element
                dynamic_cast<RBSMTetra *>( domain->giveElement( neesan ) ) ) {
                std::string name;
                // guess next available element number from current element number
                // and number of elements.
                //@todo: should be replaced by a functor that finds next elem number
                IR_GIVE_RECORD_KEYWORD_FIELD( ir, name, number );
                number     = nextElementGlobalNumber( number );
                startPoint = e->giveCellDofmanagerNumber();
                number     = RBSMTetra::makeSpringsBeam( number, startPoint, endPoint );
                springsBeams[facetExistingSisters.first - 1].insertSortedOnce( number );
            }
        }
    }
}

void RBSMTetra::postInitialize()
{
    Element::postInitialize();

    int numberOfFacets = 4;
    for ( int i = 0; i < numberOfFacets; ++i ) {
        if ( springsBeams[i].isEmpty() ) {
            continue;
        }

        //int cs = makeSpringsBeamCrossSection_3TriaDissect( i + 1 );
        int cs = makeSpringsBeamCrossSection_CircularDist( i + 1, 4 );
        for ( int sb : springsBeams[i] ) {
#if defined MINDLIN
            RBSMBeam3d *springsBeam = dynamic_cast<RBSMBeam3d *>( domain->giveElement( sb ) );
#else
            RBSBeam3d *springsBeam = dynamic_cast<RBSBeam3d *>( domain->giveElement( sb ) );
#endif
            if ( !springsBeam ) {
                OOFEM_ERROR(
                    "face %d of element %d points to an invalid element for its springs",
                    i + 1, number );
            }
            springsBeam->setCrossSection( cs );
        }
    }
}

void RBSMTetra ::updateLocalNumbering( EntityRenumberingFunctor &f )
{
    Element::updateLocalNumbering( f );
    // convert facet dof managers numbers to local
    for ( auto &fNum : facetArray ) {
        for ( int &fdNum : fNum ) {
            fdNum = f( fdNum, ERS_DofManager );
        }
    }
}

void RBSMTetra::giveCharacteristicMatrix( FloatMatrix &answer, CharType type, TimeStep *tStep )
{
    // rigid body does not have a characteristic matrix
    answer.resize( 0, 0 );
}

void RBSMTetra::giveCharacteristicVector( FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep )
{
    if ( type == ExternalForcesVector ) {
        StructuralElement::giveCharacteristicVector( answer, type, mode, tStep );
    } else if ( type == InternalForcesVector ) {
        // rigid body does not have InternalForcesVector
        answer.resize( 0 );
    } else if ( ( type == LastEquilibratedInternalForcesVector ) && ( mode == VM_Total ) ) {
        answer.resize( 0 );
    } else {
        OOFEM_ERROR( "RBSM Tetra does not support requested characteristic" )
        //StructuralElement::giveCharacteristicVector( answer, type, mode, tStep );
    }
}
/*
bool RBSMTetra::isActivated(TimeStep *tStep)
{
    return false;
}
*/


void RBSMTetra::setGeoNodesFromIr( InputRecord &ir )
// Obtains the geometry nodes of rigid body cell from the input record
{
    RBSMTetra::geoNodes.resize( numberOfCornerNodes );

    // Obtain the node IDs of the rigid body cell
    IR_GIVE_FIELD( ir, geoNodes, _IFT_Element_nodes );
}

void RBSMTetra::setCrossSection( int csIndx )
{
    Structural3DElement::setCrossSection( csIndx );
    //@todo: this may not be needed if making beams in post-initialization
    for ( IntArray ia : this->springsBeams ) {
        for ( int i : ia ) {
            if ( i > 0 ) {
                domain->elementList[i - 1]->setCrossSection( csIndx );
            }
        }
    }
}

void RBSMTetra::makeDofmanagers( InputRecord &ir )
// Makes central and cloned corner nodes for the rigid body
{
    int nextID;
    std::vector<OOFEMTXTInputRecord> rbsmInputRecords;

    // find next available DOF manager's global number
    nextID = RBSMTetra::nextDofmanagerGlobalNumber();
    // create fake input records to make RBSM DOF managers
    RBSMTetra::rbsmDummyIr( ir, rbsmInputRecords, nextID );

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

void RBSMTetra::updateClonesOfGeoNode( bool dofManArrayIsGlobal )
// adds the cloned geometry nodes of this element to 'clonesOfGeoNode' map
{
    int size = this->dofManArray.giveSize();
    // dof Manager size and geometry nodes size consistency check
    if ( size != this->geoNodes.giveSize() ) {
        OOFEM_ERROR(
            "Inconsistent size of geometry nodes and DOF manager arrays for element %d",
            this->giveGlobalNumber() );
    }

    for ( int i = 0; i < size; ++i ) {
        int geo   = this->geoNodes( i );
        int clone = this->dofManArray( i );
        // at early stage domain has not converted the global numbers to local
        clone = dofManArrayIsGlobal ? domain->giveDofManPlaceInArray( clone ) : clone;
        RBSMTetra::clonesOfGeoNode[geo].insert( clone );
    }
}

std::map<int, IntArray> RBSMTetra::updateCellElementsOfFacets()
{
    int numberOfFacets        = 4;
    int numberOfFacetVertices = 3;
    std::list<std::vector<int>> facets;
    std::map<int, IntArray> existingSisters;
    this->facetArray.resize( numberOfFacets );

    // make facets of the current element
    facets.push_back(
        { this->geoNodes( 0 ), this->geoNodes( 1 ), this->geoNodes( 2 ) } );
    facets.push_back(
        { this->geoNodes( 0 ), this->geoNodes( 1 ), this->geoNodes( 3 ) } );
    facets.push_back(
        { this->geoNodes( 0 ), this->geoNodes( 2 ), this->geoNodes( 3 ) } );
    facets.push_back(
        { this->geoNodes( 1 ), this->geoNodes( 2 ), this->geoNodes( 3 ) } );

    int index = 0;
    for ( std::vector<int> facet : facets ) {
        // sort facets vertices
        std::sort( facet.begin(), facet.end() );

        // keep a copy inside element
        // @todo: can facetArray be iterators of mapFacetElement to save memory?
        this->facetArray[index].resize( numberOfFacetVertices );
        for ( int i = 0; i < numberOfFacetVertices; ++i ) {
            this->facetArray[index][i] = facet[i];
        }

        // search for facet
        index++;
        std::map<std::vector<int>, std::set<int>>::iterator it = mapFacetElement.find( facet );
        if ( it == mapFacetElement.end() ) {
            // facet is new
            mapFacetElement[facet].insert( number );
        } else {
            // facet already exists
            for ( int neesan : it->second ) {
                existingSisters[index].insertSortedOnce( neesan );
            }
            it->second.insert( number );
        }
    }
    return existingSisters;
}

void RBSMTetra::updateCellElementsOfGeoNode()
// adds the index number of this element to 'cellElementsOfGeoNode' map
{
    int size = this->dofManArray.giveSize();
    // dof Manager size and geometry nodes size consistency check
    if ( size != this->geoNodes.giveSize() ) {
        OOFEM_ERROR(
            "Inconsistent size of geometry nodes and DOF manager arrays for element %d",
            this->giveGlobalNumber() );
    }

    for ( int i = 0; i < size; ++i ) {
        int geo = this->geoNodes( i );
        // add 'number' of current rigid body element
        RBSMTetra::cellElementsOfGeoNode[geo].insert( number );
    }
}

int RBSMTetra::makeDofmanager( InputRecord &dummyIr )
// Makes central DOF manager and returns given number
{
    int globalNumber;
    globalNumber = RBSMTetra::nextDofmanagerGlobalNumber();

    return RBSMTetra::makeDofmanager( dummyIr, globalNumber );
}

int RBSMTetra::makeDofmanager( InputRecord &dummyIr, int globalNumber )
// Makes central DOF manager and returns given number
{
    int number, nDofman = 0;
    std::string nodeType;
    Domain *d = this->giveDomain();
    std::vector<OOFEMTXTInputRecord> rbsmInputRecords;

    nDofman = d->dofManagerList.size();
    if ( nDofman <= 0 ) {
        OOFEM_ERROR( "Domain returned invalid DOF managers count: %d\n", globalNumber );
    }

    IR_GIVE_RECORD_KEYWORD_FIELD( dummyIr, nodeType, number );
    // number for central DOF manager of the rigid body
    number = nDofman + 1;

    std::unique_ptr<DofManager> newDman( classFactory.createDofManager( nodeType.c_str(), number, d ) );
    if ( !newDman ) {
        OOFEM_ERROR( "Couldn't create node of type: %s\n", _IFT_Node_Name );
    }

    newDman->initializeFrom( dummyIr );

    // make sure that global number is unique
    auto hasSameNum{
        [&globalNumber]( std::unique_ptr<oofem::DofManager> &dman ) { return dman->giveGlobalNumber() == globalNumber; }
    };
    auto it = std::find_if( d->dofManagerList.begin(), d->dofManagerList.end(), hasSameNum );
    if ( it != d->dofManagerList.end() ) {
        // the global ID already exists
        OOFEM_ERROR( "Failed to create DOF manager; the global number '%d' is already taken", globalNumber );
    }

    newDman->setGlobalNumber( globalNumber );

    d->resizeDofManagers( nDofman + 1 );

    d->setDofManager( number, std::move( newDman ) );

    /* OPTION 1
     * return the local number so **updateClonesOfGeoNode()** does not need to convert it
     * return number;
     * OPTION 2
     * return global number for consistency with dof mangers
     * return globalNumber;
     */
    return globalNumber;
}

int RBSMTetra::nextDofmanagerGlobalNumber()
// finds the next available DOF manager's global number
{
    int num, count, nDofman = 0;
    Domain *d = this->giveDomain();

    nDofman = d->dofManagerList.size();
    if ( nDofman <= 0 ) {
        OOFEM_ERROR( "Domain returned invalid DOF manager count: %d\n", nDofman );
    }

    // first try for the next DOF manager's global number
    num = nDofman + 1;

    // make sure that global number is unique
    count = 0;
    auto hasSameNum{
        //[num]( std :: unique_ptr< DofManager > dman ) { return dman->giveGlobalNumber() == num; }
        [&num]( std::unique_ptr<oofem::DofManager> &dman ) { return dman->giveGlobalNumber() == num; }
    };
    auto it = std::find_if( d->dofManagerList.begin(), d->dofManagerList.end(), hasSameNum );
    while ( it != d->dofManagerList.end() ) {
        // try to resolve duplicated number
        num++;
        count++;
        it = std::find_if( d->dofManagerList.begin(), d->dofManagerList.end(), hasSameNum );
        if ( count > nDofman ) {
            // prevent infinite loop
            OOFEM_ERROR( "Failed to find a unique node id for 'RBSM Tetra element'" );
        }
    }

    return num;
}

int RBSMTetra::nextElementGlobalNumber( int baseNumber )
// finds the next available element's global number
{
    int globalNumber, count, nElem = 0;
    Domain *d = this->giveDomain();

    nElem = d->giveNumberOfElements();
    if ( nElem <= 0 ) {
        OOFEM_ERROR( "Domain returned invalid element count: %d\n", nElem );
    }

    // first try for the next DOF manager's global number
    // @todo: this is not a good method, reliably find available elem number
    globalNumber = nElem + baseNumber;

    // make sure that global number is unique
    count = 0;
    auto hasSameNum{
        [&globalNumber]( std::unique_ptr<oofem::Element> &elem ) { return elem ? elem->giveGlobalNumber() == globalNumber : false; }
    };
    auto it = std::find_if( d->elementList.begin(), d->elementList.end(), hasSameNum );
    while ( it != d->elementList.end() ) {
        // try to resolve duplicated number
        globalNumber++;
        count++;
        it = std::find_if( d->elementList.begin(), d->elementList.end(), hasSameNum );
        if ( count > nElem + baseNumber ) {
            // prevent infinite loop
            OOFEM_ERROR( "Failed to find a unique id to make 'RBSM Tetra element'" );
        }
    }

    return globalNumber;
}


std::vector<FloatArray> RBSMTetra::coordsFromIr( InputRecord &ir )
// Obtains the coordinates of rigid body cell
{
    IntArray cornerNodes;
    std::vector<FloatArray> nodeCoords;

    nodeCoords.resize( numberOfCornerNodes + 1 );
    nodeCoords.at( 0 ).zero();

    // Obtain the node IDs of the rigid body cell
    IR_GIVE_FIELD( ir, cornerNodes, _IFT_Element_nodes );
    // Obtain node coordinates
    for ( int i = 0; i < numberOfCornerNodes; ++i ) {
        cornerNodes( i )       = domain->giveDofManPlaceInArray( cornerNodes( i ) );
        nodeCoords.at( i + 1 ) = domain->giveDofManager( cornerNodes( i ) )->giveCoordinates();
        nodeCoords.at( 0 ).add( nodeCoords.at( i + 1 ) );
    }
    // center node
    nodeCoords.at( 0 ).times( 1. / (double)numberOfCornerNodes );

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
        sprintf( buff,
            //  "RigidArmNode"  num  "Coords"   x   y   z "Master" id
            "%s %i   %s 3  %f %f %f   %s %i  "
            "dofidmask 6  1 2 3 4 5 6   "
            "mastermask 6  1 1 1 1 1 1   doftype 6  2 2 2 2 2 2",
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

double RBSMTetra::giveAreaOfFacet( int nFacet )
{
    auto facet = facetArray[nFacet - 1];
    int A      = facet.at( 1 );
    int B      = facet.at( 2 );
    int C      = facet.at( 3 );

    // A = 0.5 * |AB x AC|
    FloatArray AB =
        domain->giveDofManager( B )->giveCoordinates()
        - domain->giveDofManager( A )->giveCoordinates();
    FloatArray AC =
        domain->giveDofManager( C )->giveCoordinates()
        - domain->giveDofManager( A )->giveCoordinates();

    FloatArray ABxAC;
    ABxAC.beVectorProductOf( AB, AC );
    return 0.5 * ABxAC.computeNorm();
}

std::vector<FloatArray> RBSMTetra::giveFiberZonesOffsetsOfFacet_3TriaDissect( int nFacet )
{
    auto facet                    = facetArray[nFacet - 1];
    int numberOfRidgesAndVertices = 3;
    std::vector<FloatArray> fiberZones( numberOfRidgesAndVertices );
    FloatArray coordsFacetCentroid( numberOfRidgesAndVertices );

    // calculate facet centroid coordinates
    coordsFacetCentroid.zero();
    for ( int i = 0; i < numberOfRidgesAndVertices; ++i ) {
        coordsFacetCentroid.add( domain->giveDofManager( facet[i] )->giveCoordinates() );
    }
    coordsFacetCentroid.times( 1. / (double)numberOfRidgesAndVertices );

    // calculate fiber coordinates
    for ( int i = 0; i < numberOfRidgesAndVertices; ++i ) {
        fiberZones[i].resize( numberOfRidgesAndVertices );
        fiberZones[i].zero();

        // add coordinates of nodes of the ridge
        fiberZones[i].add( domain->giveDofManager( facet[i % 3] )->giveCoordinates() );
        fiberZones[i].add( domain->giveDofManager( facet[( i + 1 ) % 3] )->giveCoordinates() );

        // add coordinates of the centroid of the ridge
        fiberZones[i].add( coordsFacetCentroid );

        // calculate average to obtain centroid of the fiber zone
        fiberZones[i].times( 1. / 3. );

        // calculate the distance from facet centroid
        fiberZones[i].add( -1., coordsFacetCentroid );
    }

    return fiberZones;
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
    return &interpolation;
}

std::set<int> RBSMTetra::giveClonesOfGeoNode( int geoNode )
{
    std::set<int> answer;
    std::map<int, std::set<int>>::iterator it = clonesOfGeoNode.find( geoNode );
    if ( it != clonesOfGeoNode.end() ) {
        answer = it->second;
    }
    return answer;
}

std::set<int> RBSMTetra::giveCellElementsOfGeoNode( int geoNode )
{
    std::set<int> answer;
    std::map<int, std::set<int>>::iterator it = cellElementsOfGeoNode.find( geoNode );
    if ( it != cellElementsOfGeoNode.end() ) {
        answer = it->second;
    }
    return answer;
}

int RBSMTetra::makeSpringsBeam( int globalNumber, int dmanA, int dmanB )
{
    int number, nElements = 0;
    int numberOfIntPnt = 1;
    Domain *d          = this->giveDomain();
#if defined MINDLIN
    auto ELEM_TYPE     = _IFT_RBSMBeam3d_Name;
#else
    auto ELEM_TYPE     = _IFT_RBSBeam3d_Name;
#endif

    nElements = d->giveNumberOfElements();
    if ( nElements < 0 ) {
        OOFEM_ERROR( "Domain returned invalid Elements count: %d\n", globalNumber );
    }

    number = nElements + 1;

    std ::unique_ptr<Element> element( classFactory.createElement( ELEM_TYPE, number, d ) );
    if ( !element ) {
        OOFEM_ERROR( "Couldn't create spring element of type: %s", ELEM_TYPE );
    }

    // instantiate beam element
    {
        // Element
        // should develop for the case materials are different (ITZ)
        element->setMaterial( material );
        element->setCrossSection( crossSection );
        element->setDofManagers( { dmanA, dmanB } );
        //elem->setBodyLoads({});
        //elem->giveBoundaryLoadArray()->clear();
        //elem->elemLocalCS.clear();
        element->setActivityTimeFunctionNumber( 0 );
        //elem->numberOfGaussPoints = 1;
        //elem->partitions.clear();
        element->setParallelMode(Element_local);

        // Springs element

        if (
#if defined MINDLIN
            auto sprElement = dynamic_cast <RBSMBeam3d *>( element.get() )
#else
            auto sprElement = dynamic_cast <RBSBeam3d *>( element.get() )
#endif
            )
        {
            sprElement->setNumberOfGaussPoints( numberOfIntPnt );
            //elem->referenceAngle = 0; // is already set to 0
            // use condensed DOF for failed moment springs?
            //dofsToCondense = ...;
        } else {
            OOFEM_ERROR("Failed to modify Gauss points of springs beam element")
        }
    }


    // make sure that global number is unique
    auto hasSameNum{
        [&globalNumber]( std::unique_ptr<oofem::Element> &elem ) { return elem ? elem->giveGlobalNumber() == globalNumber : false; }
    };
    auto it = std::find_if( d->elementList.begin(), d->elementList.end(), hasSameNum );
    if ( it != d->elementList.end() ) {
        // the global ID already exists
        OOFEM_ERROR( "Failed to create springs beam element;"
                     "the global number '%d' is already taken",
            globalNumber );
    }

    element->setGlobalNumber( globalNumber );

    d->resizeElements( nElements + 1 );

    d->setElement( number, std::move( element ) );

    return number;
}

int RBSMTetra::makeSpringsBeamCrossSection_3TriaDissect( int nFacet )
{
    std::string csType = this->giveCrossSection()->giveClassName();
    if ( csType == "FiberedCrossSection" ) {
        int number                          = domain->giveNumberOfCrossSectionModels() + 1;
        int numberOfFibers                  = 3;
        double area                         = this->giveAreaOfFacet( nFacet );
        double fibArea                      = area / (double)numberOfFibers;
        double fibDim                       = sqrt( fibArea );
        double csDim                        = sqrt( area );
        FloatArray fiberThicks              = { fibDim, fibDim, fibDim };
        FloatArray fiberWidths              = { fibDim, fibDim, fibDim };
        double thick                        = csDim;
        double width                        = csDim;
        std::vector<FloatArray> fiberCoords = giveFiberZonesOffsetsOfFacet_3TriaDissect( nFacet );
        FloatArray fiberYcoords, fiberZcoords;
        FloatMatrix lcs;
        IntArray fiberMaterials;
        // new empty cross-section
        std::unique_ptr<RBSFiberedCrossSection> crossSection = std::make_unique<RBSFiberedCrossSection>( number, domain );
        // hard copy cross-section
        //FiberedCrossSection *rawPtr = static_cast<FiberedCrossSection *> (this->giveCrossSection());
        //std::unique_ptr<FiberedCrossSection> crossSection4 = std::make_unique<FiberedCrossSection>( *rawPtr );

        if ( springsBeams[nFacet - 1].isEmpty() ) {
            OOFEM_ERROR( "Request to make cross-section for facet %d of element %d"
                         "which does not have any springs beam element",
                nFacet, this->number );
        } else {
            // coordinates system of beams (in rare cases of having more than one sister on the same facet) on the same
            // facet could be different, but the effect is negligible so I will use the first beam.
            Element *b = domain->giveElement( springsBeams[nFacet - 1][0] );
            b->giveLocalCoordinateSystem( lcs );
        }

        for ( FloatArray coords : fiberCoords ) {
            coords.rotatedWith( lcs, 'n' );
            fiberYcoords.append( coords.at( 2 ) );
            fiberZcoords.append( coords.at( 3 ) );
        }

        // if material is provided for the element
        if ( this->material ) {
            // use provided element material for the fibers
            fiberMaterials.resize( numberOfFibers );
            for ( int &m : fiberMaterials ) {
                m = this->material;
            }
        } else {
            #if 0
            // do not modify cross section materials (for hard copy CS)
            fiberMaterials = {};
            #else
            // use material of current element (alternative way)
            fiberMaterials.resize( numberOfFibers );
            auto ip = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
            int mat = this->giveCrossSection()->giveMaterial(ip)->giveNumber();
            for ( int &m : fiberMaterials ) {
                m = mat;
            }
            #endif
        }

        // dummy linear fibers
        if ( false ) {
            numberOfFibers += 4;
            fiberMaterials.resizeWithValues( numberOfFibers );
            fiberThicks.resizeWithValues( numberOfFibers );
            fiberWidths.resizeWithValues( numberOfFibers );
            fiberYcoords.resizeWithValues( numberOfFibers );
            fiberZcoords.resizeWithValues( numberOfFibers );
            for ( int i = 0; i < 4; ++i ) {
                fiberMaterials.at( numberOfFibers - i ) = 2;
                fiberThicks.at( numberOfFibers - i )    = 1.;
                fiberWidths.at( numberOfFibers - i )    = 1.;
                fiberYcoords.at( numberOfFibers - i ) = i < 2 ? -1. : +1.;
                fiberZcoords.at( numberOfFibers - i ) = ( i == 0 || i == 3 ) ? 1. : -1.;
            }
        }

        crossSection->initializeFrom( number, fiberMaterials, fiberThicks, fiberWidths,
            numberOfFibers, thick, width, fiberYcoords, fiberZcoords );
        domain->resizeCrossSectionModels( number );
        domain->setCrossSection( number, std::move(crossSection) );

    } else {
        OOFEM_ERROR( "%s does not support %s, please use fibered section",
            this->giveClassName(), csType.c_str() )
    }

    return number;
}
int RBSMTetra::makeSpringsBeamCrossSection_CircularDist( int nFacet, int numberOfFibers )
{
    std::string csType = this->giveCrossSection()->giveClassName();
    if ( csType == "FiberedCrossSection" ) {
        int number                          = domain->giveNumberOfCrossSectionModels() + 1;
        numberOfFibers = numberOfFibers > 2 ? numberOfFibers : 3;
        double area                         = this->giveAreaOfFacet( nFacet );
        double fibArea                      = area / (double)numberOfFibers;
        double fibDim                       = sqrt( fibArea );
        double csDim                        = sqrt( area );
        double thick                        = csDim;
        double width                        = csDim;
        FloatArray fiberThicks( numberOfFibers );
        FloatArray fiberWidths( numberOfFibers );
        FloatArray fiberYcoords( numberOfFibers );
        FloatArray fiberZcoords( numberOfFibers );
        IntArray fiberMaterials( numberOfFibers );
        double fibZoneAngle                 = 2. * M_PI / (double)numberOfFibers;
        double fibZoneOffsetAngle           = fibZoneAngle / 2.;
        double fibZoneOffsetRadii           = 2. / 3. * sqrt( area * M_1_PI ); // 2/3 * sqrt(area/pi)

        for ( int i = 0; i < numberOfFibers; ++i ) {
            fiberThicks[i]  = fibDim;
            fiberWidths[i]  = fibDim;
            fiberYcoords[i] = sin( (double)i * fibZoneAngle + fibZoneOffsetAngle ) * fibZoneOffsetRadii;
            fiberZcoords[i] = cos( (double)i * fibZoneAngle + fibZoneOffsetAngle ) * fibZoneOffsetRadii;
        }

        // new empty cross-section
        std::unique_ptr<RBSFiberedCrossSection> crossSection = std::make_unique<RBSFiberedCrossSection>( number, domain );
        // hard copy cross-section
        //FiberedCrossSection *rawPtr = static_cast<FiberedCrossSection *> (this->giveCrossSection());
        //std::unique_ptr<FiberedCrossSection> crossSection = std::make_unique<FiberedCrossSection>( *rawPtr );

        // if material is provided for the element
        if ( this->material ) {
            // use provided element material for the fibers
            fiberMaterials.resize( numberOfFibers );
            for ( int &m : fiberMaterials ) {
                m = this->material;
            }
        } else {
            // use material of current element (alternative way)
            fiberMaterials.resize( numberOfFibers );
            auto ip = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint( 0 );
            int mat = this->giveCrossSection()->giveMaterial( ip )->giveNumber();
            for ( int &m : fiberMaterials ) {
                m = mat;
            }
        }

        // dummy linear fibers
        if ( false ) {
            int numberOfDummyFibers = numberOfFibers;
            double dummyScale = 0.1;
            double numberAllFibers  = numberOfFibers + numberOfDummyFibers;
            fiberMaterials.resizeWithValues( numberAllFibers );
            fiberThicks.resizeWithValues( numberAllFibers );
            fiberWidths.resizeWithValues( numberAllFibers );
            fiberYcoords.resizeWithValues( numberAllFibers );
            fiberZcoords.resizeWithValues( numberAllFibers );
            for ( int i = 0; i < numberOfDummyFibers; ++i ) {
                fiberMaterials[numberOfFibers + i] = 2;
                fiberThicks[numberOfFibers + i]    = dummyScale * fiberThicks[i];
                fiberWidths[numberOfFibers + i]    = dummyScale * fiberWidths[i];
                fiberYcoords[numberOfFibers + i]   = fiberYcoords[i];
                fiberZcoords[numberOfFibers + i]   = fiberZcoords[i];
            }
        }

        crossSection->initializeFrom( number, fiberMaterials, fiberThicks, fiberWidths,
            numberOfFibers, thick, width, fiberYcoords, fiberZcoords );
        domain->resizeCrossSectionModels( number );
        domain->setCrossSection( number, std::move(crossSection) );

    } else {
        OOFEM_ERROR( "%s does not support %s, please use fibered section",
            this->giveClassName(), csType.c_str() )
    }

    return number;
}

// S: todo: update or confirm
void RBSMTetra ::computeLumpedMassMatrix( FloatMatrix &answer, TimeStep *tStep )
// Returns the lumped mass matrix of the receiver.
{
    GaussPoint *gp;
    double dV, mss1;

    answer.resize( 12, 12 );
    answer.zero();
    gp = integrationRulesArray[0]->getIntegrationPoint( 0 );

    dV             = this->computeVolumeAround( gp );
    double density = this->giveStructuralCrossSection()->give( 'd', gp );
    mss1           = dV * density / 4.;

    for ( int i = 1; i <= 12; i++ ) {
        answer.at( i, i ) = mss1;
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
     * however the concrete ordering is element dependent (see refineMeshGlobally source if in doubts) */
/*


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

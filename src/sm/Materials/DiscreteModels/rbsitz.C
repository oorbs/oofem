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

#pragma clang diagnostic push                           //1 override clang warning
#pragma clang diagnostic ignored "-Woverloaded-virtual" //2 *
#include "rbsitz.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"
#pragma clang diagnostic pop                            //3 *

namespace oofem {
#define ZERO 1.E-6
REGISTER_Material(RBSItz);

RBSItz::RBSItz( int n, Domain *d ) : DummyMaterial( n, d )
{}


void RBSItz::initializeFrom(InputRecord &ir)
{
    DummyMaterial::initializeFrom( ir );

    IR_GIVE_FIELD(ir, mats1, _IFT_RBSItz_mat1);
    IR_GIVE_FIELD(ir, mats2, _IFT_RBSItz_mat2);
    IR_GIVE_FIELD(ir, matsi, _IFT_RBSItz_mati);

    if ( mats1.giveSize() != mats2.giveSize() ||
        mats1.giveSize() != matsi.giveSize() ) {
        OOFEM_ERROR( "The number of designated materials for"
                     "ITZ mat. %d  is not consistent",
            this->number );
    }
}


void RBSItz::giveInputRecord(DynamicInputRecord &ir)
{
    DummyMaterial::giveInputRecord(ir);

    ir.setField(this->mats1, _IFT_RBSItz_mat1 );
    ir.setField(this->mats2, _IFT_RBSItz_mat2);
    ir.setField(this->matsi, _IFT_RBSItz_mati);
}
int RBSItz::findItzMaterial( int primaryMaterial, int secondaryMaterial )
{
    for ( int i = 0; i < mats1.giveSize(); ++i ) {
        if ( this->mats1( i ) == primaryMaterial
            && this->mats2( i ) == secondaryMaterial ) {
            return this->matsi( i );
        }
    }
    // swap primary and secondary trying to find a match
    for ( int i = 0; i < mats1.giveSize(); ++i ) {
        if ( this->mats1( i ) == secondaryMaterial
            && this->mats2( i ) == primaryMaterial ) {
            return this->matsi( i );
        }
    }
    // could not find a match
    OOFEM_ERROR(
        "Could not find a valid ITZ material number between Mat %i and Mat %i",
        primaryMaterial, secondaryMaterial );
    return 0;
}

} // end namespace oofem

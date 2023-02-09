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

#ifndef rbsitz_h
#define rbsitz_h

#include "oofemlib/dummymaterial.h"
#include "intarray.h"

///@name Input fields for rbsitz
//@{
#define _IFT_RBSItz_Name "itz"
#define _IFT_RBSItz_mat1 "primarymaterials"
#define _IFT_RBSItz_mat2 "secondarymaterials"
#define _IFT_RBSItz_mati "interfacematerials"
//@}

namespace oofem {
class Domain;

/**
 * This class implements ITZ (interfacial transition zone) material (interface)
 * between two other  materials (primary and secondary)
 * @author Saeid Mehrpay
 */
class RBSItz : public DummyMaterial
{
protected:
    /// primary materials series
    IntArray mats1;
    /// secondary materials series
    IntArray mats2;
    /// interface materials series
    IntArray matsi;

public:
    RBSItz(int n, Domain * d);

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &ir) override;
    const char *giveInputRecordName() const override { return _IFT_RBSItz_Name; }
    const char *giveClassName() const override { return "RBSItz"; }

    /**
     * Gives designated material for interface of two materials
     * @param primaryMaterial first material.
     * @param secondaryMaterial second material.
     * @returns material number assigned as itz of 1st and 2nd materials,
     * returns zero for undefined conditions.
     */
    int findItzMaterial( int primaryMaterial, int secondaryMaterial );

};

/**
 * The itz interface required by classes (elements/CS) that support itz
 */
class ItzInterface : public Interface
{
protected:
    /// Interfacial material
    int itzMaterial = 0;
public:
    ItzInterface() { }

    /**
     * Gives designated ITZ material number containing information of interfaces
     * @returns ITZ material number
     */
    virtual int ItzInterface_giveItzMaterialNumber()
    {
        return this->itzMaterial;
    }
};

} // end namespace oofem
#endif // rbsitz_h

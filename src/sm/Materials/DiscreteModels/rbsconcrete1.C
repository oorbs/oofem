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

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Woverloaded-virtual"
#include "rbsconcrete1.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"
#include "sm/Elements/structuralelement.h"
#include "mathfem.h"
#pragma clang diagnostic pop

namespace oofem {
#define ZERO 1.E-6
#define CONFINED_STRESS_CONTROL     // override the stress control reduction
//#define USE_DIAGONAL_STIFFNESS    // don't use with CONFINED_STRESS_CONTROL
//#define MAKE_DIAGONAL_STIFFNESS   // only use with USE_DIAGONAL_STIFFNESS
//#define ALLOW_TMODULUS            //
//#define ALLOW_NEGATIVETMODULUS    //
REGISTER_Material(RBSConcrete1);

RBSConcrete1::RBSConcrete1(int n, Domain *d) : StructuralMaterial(n, d), D(n, d)
{}


void
RBSConcrete1::initializeFrom(InputRecord &ir)
{
    StructuralMaterial::initializeFrom( ir );

    D.initializeFrom( ir );
    IR_GIVE_FIELD( ir, this->fc, _IFT_RBSConcrete1_fc );
    IR_GIVE_FIELD( ir, this->Et, _IFT_RBSConcrete1_tangentmodulus ); ///fixme: Et should be calculated automatically
    IR_GIVE_FIELD(ir, this->nu, _IFT_IsotropicLinearElasticMaterial_n);


    if (this->fc <= 0.) {
        OOFEM_ERROR( "f'c = %d value is not valid (use positive value", this->fc );
    }

    this->linearStressRatio = .50;
    this->shearCoef         = 2.; // shear coef. 1:3
    this->criticalStrain    = 0.002; // 0.002~0.003

    this->ft = fc <= 50. ? 0.3 * pow( fc, .6666667 ) : 2.12 * log( 1. + 0.1 * ( fc + 8. ) );
    this->fs = 0.5 * this->fc / shearCoef;
    this->E  = this->D.giveYoungsModulus();
    this->H  = E * Et / ( E - Et );
    this->G  = D.giveShearModulus();

    fs_k.resize(maxNK);
    G_k.resize( maxNK );
    H_k.resize( maxNK );
    eps_k.resize( maxNK + 1 );
    epsP_k.resize( maxNK + 1 );

    // shear yield stress
    this->fs_k = {
        0.5 * this->fc / this->shearCoef * linearStressRatio,
        0.5 * this->fc / this->shearCoef * 1.0,
        0.5 * this->fc / this->shearCoef * 0.5,
        0.
    };

    ///TODO: this one should be inside the status since it depends on GP
    {
        double maxLinearStress      = this->linearStressRatio * this->fc;
        double maxLinearStrain      = maxLinearStress / this->E;
        //double peakShearStress      = 0.5 * this->fc / shearCoef;
        double maxLinearShearStrain = ( 1. + nu ) / 2. * maxLinearStrain;
        double maxLinearShearStress = this->linearStressRatio * fs;
        //double maxLinearShearStrain = maxLinearShearStress / G;
#if 0   // ignore nonlinear shear strains to calculate critical shear strain
        double criticalShearStrain  = ( 1. + nu ) / 2. * this->criticalStrain;
        // corresponds to CEB-FIP's Ec1 (slope of origin to peak)
        //double Gc1 = fs / criticalShearStrain;
#else   // consider nonlinear shear strains
        double criticalShearStrain;
        {
            double criticalShearStrainE, criticalShearStrainP;
            /// TODO: better to update this for Dr. Nagai's PR effect
            double PR_temp       = ( 2 + this->nu ) / ( 2 + 2 * this->nu ); // temporary fix for effect of PR on E_Mac
            criticalShearStrainE = this->fs / this->G;
            // alternatively: criticalShearStrainE = 0.5 * this->fc / ( this->shearCoef * this->G );
            // or for Gsc=0.5, Gs=2: criticalShearStrainE = ( 1. + nu ) / 2. * ( fc / E );
            criticalShearStrainP = ( this->criticalStrain - this->fc / ( this->E * PR_temp ) );
            // alternatively: criticalShearStrainP = criticalStrain - 2. / ( 1. + this->nu ) * criticalShearStrainE
            criticalShearStrain  = criticalShearStrainE + criticalShearStrainP;
        }
#endif
        double dStrain = criticalShearStrain - maxLinearShearStrain;
        double Gt1;
        if ( dStrain < ZERO ) {
            OOFEM_WARNING(
                "The stiffness of the material is too low to provide hardening before the compressive strain of %f",
                criticalStrain );
            Gt1 = 0.;
        } else {
            Gt1 = ( fs - maxLinearShearStress ) / dStrain ;
        }
        // nonlinear shear (hardening/softening) moduli
        this->G_k = {
            G,
            Gt1,
            -500.,
            -50.
        };
    }
    // calculate other multilinear specifications
    eps_k( 0 )  = fs_k( 0 ) / G_k( 0 );
    epsP_k( 0 ) = 0;
    for ( int i = 1; i < maxNK; ++i ) {
        eps_k( i )  = eps_k( i - 1 ) + ( fs_k( i ) - fs_k( i - 1 ) ) / G_k( i );
        epsP_k( i ) = eps_k( i  ) - fs_k( i ) / G_k( 0 );
        if ( eps_k( i ) != eps_k( i ) || epsP_k( i ) != epsP_k( i ) ){
            OOFEM_ERROR("Unable to calculate concrete yield stress due to divide to zero error.")
        }
        H_k( i - 1 ) = ( fs_k( i ) - fs_k( i - 1 ) ) / ( epsP_k( i ) - epsP_k( i - 1 ) );
    }
    eps_k( maxNK )  = 1.e+16; // INFINITY
    epsP_k( maxNK ) = 1.e+16; // INFINITY
    H_k( maxNK - 1 ) = 0.;
    //
}


void RBSConcrete1::giveInputRecord(DynamicInputRecord &ir)
{
    StructuralMaterial::giveInputRecord(ir);
    D.giveInputRecord(ir);
    ir.setField(this->fc, _IFT_RBSConcrete1_fc );
    ir.setField(this->Et, _IFT_RBSConcrete1_tangentmodulus);
    //ir.setField(this->H, _IFT_RBSConcrete1_hardeningmoduli);
}


MaterialStatus *RBSConcrete1::CreateStatus(GaussPoint *gp) const
{
    return new RBSConcrete1Status(gp);
}


FloatArrayF<6> RBSConcrete1::giveRealStressVector_3d( const FloatArrayF<6> &totalStrain, GaussPoint *gp, TimeStep *tStep ) const
{
    auto status = static_cast<RBSConcrete1Status *>( this->giveStatus( gp ) );

    // subtract stress thermal expansion
    auto thermalStrain = this->computeStressIndependentStrainVector_3d( gp, tStep, VM_Total );
    auto strain = totalStrain - thermalStrain;

    // subtract plastic strains
    auto plasticStrain = status->givePlasticStrain();
    FloatArrayF<6> trialElasticStrain;
    // normal elastic strain
    for ( int i = 1; i <= 3; ++i ) { //todo: update this if focal point is used
        trialElasticStrain.at( i ) = strain.at( i ) - plasticStrain.at( i );
    }
    // shear elastic strain
    for ( int i = 4; i <= 6; ++i ) {
        trialElasticStrain.at( i ) = strain.at( i ) - sgn( strain.at( i ) ) * plasticStrain.at( i );
    }

#if !defined USE_DIAGONAL_STIFFNESS // use elastic isotropic stiffness by ref
    const auto &elasticStiffness = D.giveTangent();
#elif defined MAKE_DIAGONAL_STIFFNESS // make diagonal stiffness matrix
    const FloatMatrixF<6,6> elasticStiffness = {
        this->E, 0., 0., 0., 0., 0.,
        0., this->E, 0., 0., 0., 0.,
        0., 0., this->E, 0., 0., 0.,
        0., 0., 0., this->G, 0., 0.,
        0., 0., 0., 0., this->G, 0.,
        0., 0., 0., 0., 0., this->G
    };
#else // use elastic isotropic [D] but make it diagonal
    auto elasticStiffness = D.giveTangent();
    for ( size_t i = 1; i <= elasticStiffness.rows(); ++i ) {
        for ( size_t j = 1; j <= elasticStiffness.cols(); ++j ) {
            if ( i != j ) elasticStiffness.at( i, j ) = 0.;
        }
    }
#endif

    auto trialStress = dot( elasticStiffness, trialElasticStrain );

    // Trial stresses
    double trialNormalStress = trialStress.at( 1 );
    double trialShearStress1 = trialStress.at( 5 );
    double trialShearStress2 = trialStress.at( 6 );

    // Stress signs
    auto signN0 = sgn( trialElasticStrain.at( 1 ) );
    auto signS1 = sgn( trialElasticStrain.at( 5 ) );
    auto signS2 = sgn( trialElasticStrain.at( 6 ) );
    // Check consistency of stress signs
    if ( signN0 * sgn( trialNormalStress ) < 0 ) {
#ifndef USE_DIAGONAL_STIFFNESS
        signN0 = sgn( trialNormalStress );
#else
        OOFEM_ERROR("Inconsistent calculated stress sign for normal spring N at step %d:%d",
            tStep->giveNumber(), tStep->giveSubStepNumber());
#endif
    } else if ( signS1 * sgn( trialShearStress1 ) < 0 ) {
        OOFEM_ERROR("Inconsistent calculated stress sign for shear spring S1");
    } else if ( signS2 * sgn( trialShearStress2 ) < 0 ) {
        OOFEM_ERROR("Inconsistent calculated stress sign for shear spring S2");
    }

    //auto [devTrialStress, meanTrialStress] = computeDeviatoricVolumetricSplit(); // c++17
    //double J2 = this->computeSecondStressInvariant(devTrialStress);

    FloatArrayF<6> stress;
    stress = trialStress;

    /// Spring nonlinear curve state
    int nKn0        = status->giveNormalState();
    int nKs1        = status->giveShearState1();
    int nKs2        = status->giveShearState2();
    double k = status->giveK();
    double ks1 = status->giveKs1();
    double ks2 = status->giveKs2();

    // 1. Normal stress correction

    // evaluate the yield surface
    double sigma_yt = this->ft + H * k;    ///fixme!
    double tr_f = trialNormalStress - sigma_yt;

    if ( nKn0 ) {
        // damaged normal spring
        stress.at( 1 ) = 0;
#ifndef USE_DIAGONAL_STIFFNESS
        // to prevent problems with inner iterative method of structural material
        stress.at( 2 ) = 0.;
        stress.at( 3 ) = 0.;
#endif
    } else if ( tr_f <= ZERO ) { // elastic
        //status->letTempPlasticStrainBe( status->givePlasticStrain() );
    } else { // plastic loading
        // plastic strain inc.
        double dPlStrain = tr_f / ( E + H );
        // radial return
        auto corNormalStress = trialNormalStress - E * dPlStrain;
        if ( corNormalStress < 0. ) {
            corNormalStress = 0.;
            nKn0            = 1;
#ifndef USE_DIAGONAL_STIFFNESS
            // to prevent problems with inner iterative method of structural material
            stress.at( 2 ) = 0.;
            stress.at( 3 ) = 0.;
#endif
        }
        stress.at( 1 ) = corNormalStress;
        k += dPlStrain;

        //auto plasticStrain = status->givePlasticStrain();
        plasticStrain.at( 1 ) += dPlStrain;
        // FOR NON-DIAGONAL [D]: this is wrong for PR is calculated from stress vect
        // + remove the Poisson effect's strains by making them part of plastic strain
        // plasticStrain.at( 2 ) -= this->nu * dPlStrain;
        // plasticStrain.at( 3 ) -= this->nu * dPlStrain;
        // - another option is to remove these strains from strain vector:
        // totalStrain.at( 2 ) -= this->nu * dPlStrain;
        // totalStrain.at( 3 ) -= this->nu * dPlStrain;
    }


    // 2. Shear stress correction

    if ( status->giveNormalState() ) {
        // damaged normal spring
        stress.at( 5 ) = 0;
        stress.at( 6 ) = 0;
    } else {
        // Shear spring 1:
        {
            double sigma_ys1, tr_fs1;
            sigma_ys1 = fs_k( nKs1 ) + H_k( nKs1 ) * ( ks1 - epsP_k( nKs1 ) );   ///fixme!
            tr_fs1    = fabs( trialShearStress1 ) - sigma_ys1;
            //
            if ( nKn0 ) { // if normal spring has failed, there is no shear stress
                stress.at( 5 ) = 0;
            } else if ( tr_fs1 <= ZERO ) { // elastic
                // status->letTempPlasticStrainBe( status->givePlasticStrain() );
            } else { // plastic loading
                double dPlStrain = tr_fs1 / ( G + H_k( nKs1 ) ); // plastic multiplier
                // radial return
                auto corShearStress = trialShearStress1 - sgn( trialShearStress1 ) * G * dPlStrain;
                // Making sure that the corrected stresses make sense
                if ( sgn( signS1 * corShearStress ) < 0 ) {
                    double nextElStrainMagnitude_s1 = fabs(strain.at( 5 )) - ( plasticStrain.at( 5 ) + dPlStrain );
                    if ( nextElStrainMagnitude_s1 < 0. && H_k( nKs1 ) < 0. ) {
                        // plastic strain is larger than strain due to softening (H_k < 0)
                        dPlStrain      = dPlStrain - fabs( corShearStress ) / G;
                        corShearStress = 0.;
                        nextElStrainMagnitude_s1 = fabs( strain.at( 5 ) ) - ( plasticStrain.at( 5 ) + dPlStrain );
#ifdef VERBOSE
                        if ( fabs( nextElStrainMagnitude_s1 ) > ZERO ) { // zero=~0.
                            OOFEM_WARNING( "Invalid plastic strain for shear S1@El%d.Gp%d",
                                gp->giveElement()->giveGlobalNumber(), gp->giveNumber() );
                        }
#endif
                        //normalState++;
                    } else {
                        OOFEM_ERROR( "Invalid stress for shear S1@El%d",
                            gp->giveElement()->giveGlobalNumber() );
                    }
                }
                ks1 += dPlStrain;
                plasticStrain.at( 5 ) += dPlStrain;

                int c = 0;
                // this can be a loop "WHILE" or deliberately done once "IF"
                while ( ks1 > epsP_k( nKs1 + 1 ) ) {
                    if ( ++c > 100 ) {
                        OOFEM_ERROR("Plastic stress S1 calculation failed after %i corrections", c);
                    }
                    nKs1++;
                    // Evaluate new yield stress from extra plastic strain
                    sigma_ys1 = fs_k( nKs1 ) + H_k( nKs1 ) * ( ks1 - epsP_k( nKs1 ) );
                    tr_fs1    = fabs( corShearStress ) - sigma_ys1;
                    // Calculate complimentary part of S1 plastic strain
                    dPlStrain      = tr_fs1 / ( G + H_k( nKs1 ) );
                    corShearStress = corShearStress - signS1 * G * dPlStrain;
                    // Making sure that the corrected stresses make sense
                    if ( sgn( signS1 * corShearStress ) < 0 ) {
                        double nextElStrainMagnitude_s1 = fabs(strain.at( 5 )) - ( plasticStrain.at( 5 ) + dPlStrain );
                        if ( nextElStrainMagnitude_s1 < 0. && H_k( nKs1 ) < 0. ) {
                            // plastic strain is larger than strain due to softening (H_k < 0)
                            dPlStrain      = dPlStrain - fabs( corShearStress ) / G;
                            nextElStrainMagnitude_s1 = fabs(strain.at( 5 )) - ( plasticStrain.at( 5 ) + dPlStrain );
                            if ( fabs(nextElStrainMagnitude_s1) < ZERO) { // zero=~0.
                                corShearStress = 0.;
                            } else {
                                OOFEM_ERROR( "Invalid elastic strain for shear spring S1" );
                            }
                            // normalState++;
                        } else {
                            OOFEM_ERROR( "Invalid shear stress / strain for spring S1" );
                        }
                    }

                    ks1 += dPlStrain;
                    plasticStrain.at( 5 ) += dPlStrain;
                }
                if ( corShearStress * trialShearStress1 < 0. ) {
                    //corShearStress = 0.;
                    //normalState++;
                    OOFEM_ERROR( "Invalid shear stress for spring S1" );
                }
                stress.at( 5 ) = corShearStress;
            }
        }
        // Shear spring 2:
        {
            double sigma_ys2, tr_fs2;
            sigma_ys2 = fs_k( nKs2 ) + H_k( nKs2 ) * ( ks2 - epsP_k( nKs2 ) );
            tr_fs2    = fabs( trialShearStress2 ) - sigma_ys2;
            //
            if ( nKn0 ) { // if normal spring has failed, there is no shear stress
                stress.at( 6 ) = 0;
            } else if ( tr_fs2 <= ZERO ) { // elastic
                // status->letTempPlasticStrainBe( status->givePlasticStrain() );
            } else { // plastic loading
                double dPlStrain = tr_fs2 / ( G + H_k( nKs2 ) ); // plastic multiplier
                // radial return
                auto corShearStress = trialShearStress2 - sgn( trialShearStress2 ) * G * dPlStrain;
                // Making sure that the corrected stresses make sense
                if ( sgn( signS2 * corShearStress ) < 0 ) {
                    double nextElStrainMagnitude_s2 = fabs(strain.at( 6 )) - ( plasticStrain.at( 6 ) + dPlStrain );
                    if ( nextElStrainMagnitude_s2 < 0. && H_k( nKs2 ) < 0. ) {
                        // plastic strain is larger than strain due to softening (H_k < 0)
                        dPlStrain      = dPlStrain - fabs( corShearStress ) / G;
                        corShearStress = 0.;
                        nextElStrainMagnitude_s2 = fabs(strain.at( 6 )) - ( plasticStrain.at( 6 ) + dPlStrain );
#ifdef VERBOSE
                        if ( fabs(nextElStrainMagnitude_s2) > ZERO) { // zero=~0.
                            OOFEM_WARNING( "Invalid plastic strain for shear S2@El%d.Gp%d",
                                gp->giveElement()->giveGlobalNumber(), gp->giveNumber() );
                        }
#endif
                        // normalState++;
                    } else {
                        OOFEM_ERROR( "Invalid stress for shear S2@El%d",
                                gp->giveElement()->giveGlobalNumber() );
                    }
                }
                ks2 += dPlStrain;
                plasticStrain.at( 6 ) += dPlStrain;

                int c = 0;
                // this can be a loop "WHILE" or deliberately done once "IF"
                while ( ks2 > epsP_k( nKs2 + 1 ) ) {
                    if ( ++c > 100 ) {
                        OOFEM_ERROR("Plastic stress S2 calculation failed after %i corrections", c);
                    }
                    nKs2++;
                    // Evaluate new yield stress from extra plastic strain
                    sigma_ys2 = fs_k( nKs2 ) + H_k( nKs2 ) * ( ks2 - epsP_k( nKs2 ) );
                    tr_fs2    = fabs( corShearStress ) - sigma_ys2;
                    // Calculate complimentary part of S2 plastic strain
                    dPlStrain      = tr_fs2 / ( G + H_k( nKs2 ) );
                    corShearStress = corShearStress - signS2 * G * dPlStrain;
                    // Making sure that the corrected stresses make sense
                    if ( sgn( signS2 * corShearStress ) < 0 ) {
                        double nextElStrainMagnitude_s2 = fabs(strain.at( 6 )) - ( plasticStrain.at( 6 ) + dPlStrain );
                        if ( nextElStrainMagnitude_s2 < 0. && H_k( nKs2 ) < 0. ) {
                            // plastic strain is larger than strain due to softening (H_k < 0)
                            dPlStrain      = dPlStrain - fabs( corShearStress ) / G;
                            nextElStrainMagnitude_s2 = fabs(strain.at( 6 )) - ( plasticStrain.at( 6 ) + dPlStrain );
                            if ( fabs(nextElStrainMagnitude_s2) < ZERO) { // zero=~0.
                                corShearStress = 0.;
                            } else {
                                OOFEM_ERROR( "Invalid elastic strain for shear spring S2" );
                            }
                            // normalState++;
                        } else {
                            OOFEM_ERROR( "Invalid shear stress / strain for spring S2" );
                        }
                    }

                    ks2 += dPlStrain;
                    plasticStrain.at( 6 ) += dPlStrain;
                }
                if ( corShearStress * trialShearStress2 < 0. ) {
                    //corShearStress = 0.;
                    //normalState++;
                    OOFEM_ERROR( "Invalid shear stress for spring S2" );
                }
                stress.at( 6 ) = corShearStress;
            }
        }
    }

    // 3. Store the temporary values for the given iteration
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(stress);
    status->letTempPlasticStrainBe( plasticStrain );
    status->letTempKBe(k);
    status->letTempKs1Be(ks1);
    status->letTempShearState1(nKs1);
    status->letTempKs2Be(ks2);
    status->letTempShearState2(nKs2);
    //status->letTempDevTrialStressBe(devTrialStress);
    status->letTempNormalStateBe( nKn0 );
    return stress;
}

// (v-2.5) possibly not needed, delete
FloatArrayF<6> RBSConcrete1::giveRealStressVector_3dBeamSubSoil( const FloatArrayF<6> &reducedStrain, GaussPoint *gp, TimeStep *tStep ) const
{
    auto status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    FloatArray vE, vS;
    FloatMatrix tangent, reducedTangent;

    // Initial guess;
    vE = status->giveStrainVector();
    vS = this->giveRealStressVector_3d(vE, gp, tStep);

    return vS;
}

// (v-2.5) keep confined stresses for redistribution in RBSM element
FloatArray RBSConcrete1::giveRealStressVector_StressControl(
    const FloatArray &reducedStrain, const IntArray &strainControl, GaussPoint *gp, TimeStep *tStep ) const
{
#ifdef CONFINED_STRESS_CONTROL
    auto status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    IntArray stressControl;
    FloatArray vE, increment_vE, vS;
    FloatMatrix tangent, reducedTangent;
    // Stores reduced strain vector
    FloatArray answer;

    stressControl.resize(6 - strainControl.giveSize() );
    for ( int i = 1, j = 1; i <= 6; i++ ) {
        if ( !strainControl.contains(i) ) {
            stressControl.at(j++) = i;
        }
    }

    // Initialize strain vector
    vE = status->giveStrainVector();
    for ( int i = 1; i <= strainControl.giveSize(); ++i ) {
        vE.at(strainControl.at(i) ) = reducedStrain.at(i);
    }

    // Calculate reduced stress vector
    vS = this->giveRealStressVector_3d(vE, gp, tStep);

    // Pick out the (response) stresses for the controlled strains
    answer.beSubArrayOf(vS, strainControl);

    // the default implementation is to eliminate stress control,
    // reducedvS.beSubArrayOf(vS, stressControl);
    // do while ( reducedvS.computeNorm() <= 1e-6 * vS.computeNorm() )
    // {
    //    tangent = this->give3dMaterialStiffnessMatrix(TangentStiffness, gp, tStep);
    //    reducedTangent.beSubMatrixOf(tangent, stressControl, stressControl);
    //    reducedTangent.solveForRhs(reducedvS, increment_vE);
    //    increment_vE.negated();
    //    vE.assemble(increment_vE, stressControl);
    //    etc.
    // }
    // but for RBSM we keep these stresses (confined condition):
    return answer;
#else
    return StructuralMaterial::giveRealStressVector_StressControl(reducedStrain, strainControl, gp, tStep);
#endif
}


FloatMatrixF<6,6>
RBSConcrete1::give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{

#if !defined USE_DIAGONAL_STIFFNESS // use elastic isotropic stiffness
    auto elasticStiffness = D.giveTangent();
#elif defined MAKE_DIAGONAL_STIFFNESS // make diagonal stiffness matrix
    FloatMatrixF<6,6> elasticStiffness = {
        this->E, 0., 0., 0., 0., 0.,
        0., this->E, 0., 0., 0., 0.,
        0., 0., this->E, 0., 0., 0.,
        0., 0., 0., this->G, 0., 0.,
        0., 0., 0., 0., this->G, 0.,
        0., 0., 0., 0., 0., this->G
    };
#else // use elastic isotropic [D] but make it diagonal
    auto elasticStiffness = D.giveTangent();
    for ( size_t i = 1; i <= elasticStiffness.rows(); ++i ) {
        for ( size_t j = 1; j <= elasticStiffness.cols(); ++j ) {
            if ( i != j ) elasticStiffness.at( i, j ) = 0.;
        }
    }
#endif

#ifdef ALLOW_TMODULUS
    auto status = static_cast< RBSConcrete1Status * >( this->giveStatus(gp) );

    double trialNormalStress    = status->giveStressVector().at( 1 );
    double trialShearStress1    = status->giveStressVector().at( 5 );
    double trialShearStress2    = status->giveStressVector().at( 6 );

    double k = status->giveK();
    double sigma_y = this->ft + H * k;

    // 1. Normal spring

    // evaluate the yield surface
    double tr_f = trialNormalStress - sigma_y;

    if(status->giveNormalState()) {
        // prevent instability by not using 0.
        elasticStiffness.at( 1, 1 ) = 1.; // almost zero
    } else if ( tr_f < -ZERO ) {
        // elastic loading
    } else {
        // plastic loading: Et
        // 1. use Et > 0 or small positive value
        /*elasticStiffness.at(1, 1 ) = max(this->Et,
            min( 1., elasticStiffness.at( 1, 1 ) / 10000. ));*/
        // 2. use 1.
        /*elasticStiffness.at(1, 1 ) = 1.;*/
#ifdef ALLOW_NEGATIVETMODULUS
        // 3. use Et
        elasticStiffness.at(1, 1 ) = this->Et;
#else
        // 4. use Et > 0 or Ee
        elasticStiffness.at(1, 1 ) = this->Et > 0. ? this->Et : this->E;
#endif
    }

    // 2. Shear springs

    if ( status->giveNormalState() ) {
        // damaged normal spring
        elasticStiffness.at( 5, 5 ) = 1.;   // almost zero
        elasticStiffness.at( 6, 6 ) = 1.;   // almost zero
    } else {
        int nKs1        = status->giveShearState1();
        int nKs2        = status->giveShearState2();
        double ks1 = status->giveKs1();
        double ks2 = status->giveKs2();
        //double sigma_ys1 = this->fc + Hs1 * ks1;
        double sigma_ys1 = fs_k( nKs1 ) + H_k( nKs1 ) * ( ks1 - epsP_k( nKs1 ) );       ///FIXME!
        double tr_fs1 = fabs( trialShearStress1 ) - sigma_ys1;
        double sigma_ys2 = fs_k( nKs2 ) + H_k( nKs2 ) * ( ks2 - epsP_k( nKs2 ) );       ///FIXME!
        double tr_fs2 = fabs( trialShearStress2 ) - sigma_ys2;

        // Shear 1
        if ( tr_fs1 <= -ZERO ) { // elastic loading
            // continue
        } else { // plastic loading;
#ifdef ALLOW_NEGATIVETMODULUS
            // 1. use Gt
            elasticStiffness.at( 5, 5 ) = this->G_k( nKs1 );
#else
            // 2. use Gt for hardening (faster convergence) and elastic G for softening (avoid divergence)
            elasticStiffness.at( 5, 5 ) = this->G_k( nKs1 ) > 0. ? this->G_k( nKs1 ) : G;
#endif
        }

        // Shear 2
        if ( tr_fs2 <= -ZERO
            ) { // elastic loading
            // continue
        } else { // plastic loading
#ifdef ALLOW_NEGATIVETMODULUS
            // 1. use Gt
            elasticStiffness.at( 6, 6 ) = this->G_k( nKs2 );
#else
            // 2. use Gt for hardening (faster convergence) and elastic G for softening (avoid divergence)
            elasticStiffness.at( 6, 6 ) = this->G_k( nKs2 ) > 0. ? this->G_k( nKs2 ) : G;
#endif
        }
    }
#endif

    return elasticStiffness;
}


FloatArrayF<6>
RBSConcrete1::giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const
{
    return D.giveAlpha();
}


int
RBSConcrete1::giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    RBSConcrete1Status *status = static_cast< RBSConcrete1Status * >( this->giveStatus(gp) );
    if ( type == IST_PlasticStrainTensor ) {
        answer = status->givePlasticStrain();
        return 1;
    } else {
        return StructuralMaterial::giveIPValue(answer, gp, type, tStep);
    }
}


//=============================================================================


RBSConcrete1Status::RBSConcrete1Status(GaussPoint * g) :
    StructuralMaterialStatus(g)
{
    strainVector.resize(6);
    stressVector.resize(6);
    tempStressVector = stressVector;
    tempStrainVector = strainVector;
}

void
RBSConcrete1Status::initTempStatus()
{
    //StructuralMaterialStatus::initTempStatus();

    // reset temp vars.
    tempStressVector = stressVector;
    tempStrainVector = strainVector;
    tempPVector      = PVector;
    tempFVector      = FVector;


    tempPlasticStrain = plasticStrain;
    tempK = k;
    tempNormalState = normalState;
    tempDevTrialStress = zeros<6>();
}


void 
RBSConcrete1Status::updateYourself(TimeStep *tStep)
{
    // Copy the temporary values to the convered ones. 
    // This method is called after equilibrium has been reached and should 
    // set all...
    StructuralMaterialStatus::updateYourself(tStep);

    plasticStrain = tempPlasticStrain;
    k   = tempK;
    ks1 = tempKs1;
    ks2 = tempKs2;
    normalState   = tempNormalState;
    shearState1   = tempShearState1;
    shearState2   = tempShearState2;
    // deviatoric trial stress is not really a state variable and was used not to repeat some code...
}

} // end namespace oofem

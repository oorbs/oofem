//                      ____ ____  ____  ____ _____
//                     / __   __ \/ __ \/ __ ) ___/
//                    / / /  / / / /_/ / __  \__ \
//                   / /_/  /_/ / _, _/ /_/ /__/ /
//                   \____ ____/_/ |_/_____/____/
//
// 3D prism mesoscopic prism specimen generator
// File author: Saeid Mehrpay
// Script Language: GMSH geo script, 
//     (for more information see https://gmsh.info/)
// Project: OORBS 3D-RBSM - Shenzhen University 2022
//
/// -----------------------------------------
// .: Options and Workspace :.
/// -----------------------------------------
//General.NumThreads = 2;
// adjust value here for correct merge result
Geometry.Tolerance = 1.;
Geometry.MatchMeshTolerance = Geometry.Tolerance;
Coherence Mesh;
General.AbortOnError = 1;
SetFactory("OpenCASCADE");
///
/// -----------------------------------------
// .: MESHING ALGORITHM FOR Surface :.
/// -----------------------------------------
h = 1.; // local meshing dimension
Mesh.MeshSizeMin = 1.;
Mesh.MeshSizeMax = 30.;
//MeshAlgorithm Surface {17} = 2;
//Mesh.Voronoi = 1;
//Mesh.MeshSizeFactor = 1.;
userMeshSizeofAllGeometry = 1.; // new (set to 0 to inactivate)
If(userMeshSizeofAllGeometry > 0)
    MeshSize{ PointsOf{ Volume{:}; } } = userMeshSizeofAllGeometry;
EndIf
/// not effective:
//set global mesh dimension = 0.2; /// not valid anymore
//Mesh.ScalingFactor = 0.2;
///
// devisions of circle element
Mesh.MinimumCirclePoints = 3;
//Mesh.CharacteristicLengthExtendFromBoundary = 0;
//Mesh.CharacteristicLengthFromPoints = 0;
//Mesh.CharacteristicLengthFromCurvature = 0;
///
Mesh.RefineSteps = 20;
Mesh.Smoothing = 10; //1; smoothing steps to final mesh
/// -----------------------------------------
// .: Specimen Dimensions :.
/// -----------------------------------------
height = 200.;
width  = 100.;
depth  = 100.;
// Min Aggregates - Mold gap
wallEffect = .5; // 1.5
// wall effect between aggregates
aggwalleff = 5.; // 2.
// gap error tolerance
Tol = Geometry.Tolerance * 0.1;
///
moveMax       = 8;      // allowed hits for existing aggregate before calling macro
move_tryMax   = 24;     // move macro's max tries tolerance
shiftCountMax = 24;     // for maximum shift tries
debug         = 1;      // print debug information
verbose       = 0;      // print full detailed information
conservative  = 0;      // check for modeling errors occasionally
///
/// -----------------------------------------
// .: Number of sieves divisions :.
/// -----------------------------------------
// For aggregate's radii use R(sieveindex)
// For aggregates' count use numAgg(sieveNum)
sieveNum = 1; // use 6 of 10
D      = {20.80, 17.50, 14.75, 11.95, 10.36,  8.76,  7.37,  6.54,  6.01,  5.21};//
numAgg = {   16,    26,    40,    52,   100,   165,   261,   120,   299,   655};//
///
// convert Diameter to Radii
R(sieveNum) = 0;
For i In  {0:sieveNum-1}
    R[i] = D[i] * 0.5;
EndFor
//
Sleep 0.02; //pause
Draw;
/// -----------------------------------------
// .: View Options :.
/// -----------------------------------------
General.AxesMikado= 1;
General.AxesTicsX = 6;
General.AxesTicsY = 11;
General.AxesTicsZ = 6;
/// =========================================
/// =========================================
/// 
/// > Making Specimen <
iPnt = 1;
iLin = 1;
Point(iPnt++) = {0.,    0.,     0.,    h};
Point(iPnt++) = {width, 0.,     0.,    h};
Line(iLin++) = {iPnt-2, iPnt-1};
Point(iPnt++) = {width, height, 0.,    h};
Line(iLin++) = {iPnt-2, iPnt-1};
Point(iPnt++) = {0.,    height, 0.,    h};
Line(iLin++) = {iPnt-2, iPnt-1};
Line(iLin++) = {iPnt-1, iPnt-4};
///
Point(iPnt++) = {0.,    0.,     depth, h};
Point(iPnt++) = {width, 0.,     depth, h};
Line(iLin++) = {iPnt-2, iPnt-1};
Point(iPnt++) = {width, height, depth, h};
Line(iLin++) = {iPnt-2, iPnt-1};
Point(iPnt++) = {0.,    height, depth, h};
Line(iLin++) = {iPnt-2, iPnt-1};
Line(iLin++) = {iPnt-1, iPnt-4};
///
Line(iLin++) = {iPnt-8, iPnt-4};
Line(iLin++) = {iPnt-7, iPnt-3};
Line(iLin++) = {iPnt-6, iPnt-2};
Line(iLin++) = {iPnt-5, iPnt-1};
///
/// ----------------------------------------------
keep_flag = 0; // if aggregate was moved outside specimen is true
keep_target = 0; // checking target
keep_radii = 0;  // radious of target
keep_error = 0;
Macro Keep_Macro
        //Printf("    .: Keep Macro :. " );
        keep_flag = 0; //okay
        If ( x[keep_target] < keep_radii + wallEffect || x[keep_target] > width  - keep_radii - wallEffect || 
             y[keep_target] < keep_radii + wallEffect || y[keep_target] > height - keep_radii - wallEffect || 
             z[keep_target] < keep_radii + wallEffect || z[keep_target] > depth  - keep_radii - wallEffect )
            keep_flag = 1;
            // alternative:
            //Printf("Aggregated shifted out of specimen, replacing aggregate...");
            //confLen = 0.; // means new random coordinates to be calculated
            // or:
            If(keep_error)
                Error ( "Aggregate %g was placed out of specimen", keep_target+1);
                Abort;
            Else
                Printf("  @ Aggregate %g shifted out of specimen, will be placed back on boundary...", keep_target+1);
            EndIf
            If (x[keep_target] < keep_radii + wallEffect) x[keep_target] = keep_radii + wallEffect; EndIf
            If (y[keep_target] < keep_radii + wallEffect) y[keep_target] = keep_radii + wallEffect; EndIf
            If (z[keep_target] < keep_radii + wallEffect) z[keep_target] = keep_radii + wallEffect; EndIf
            If (x[keep_target] > width  - keep_radii - wallEffect) x[keep_target] = width  - keep_radii - wallEffect; EndIf
            If (y[keep_target] > height - keep_radii - wallEffect) y[keep_target] = height - keep_radii - wallEffect; EndIf
            If (z[keep_target] > depth  - keep_radii - wallEffect) z[keep_target] = depth  - keep_radii - wallEffect; EndIf
        EndIf
        //Printf("                            ");
Return
/// -----------------------------------------
///
///
/// ----------------------------------------------
check_Target = 0;
Macro Check_Macro
        If(verbose) Printf("    .: Check Macro :. " ); EndIf
        If(verbose) Printf("    Check target aggregate %g : [%g %g %g]", check_Target+1, x[check_Target], y[check_Target], z[check_Target] ); EndIf
        Check_C2 = 0;
        /// find radii of target
        Check_TargetSieve = 0;
        temp = -1; // it's 0-based
        Check_Flag = 0;
        For M1 In {0:N} //sieves
            temp += numAgg[M1]; // last aggregate in sieve M1
            If(verbose) Printf("    N=%g M1=%g Aggregates in Check Sieve=%g temp=%g, check_Target=%g", N, M1, numAgg[M1], temp, check_Target); EndIf
            If(temp >= check_Target)
                Check_TargetSieve = M1;
                If(verbose) Printf("    Target Sieve #: %g, R: %g mm Current-Last Sieve (N)=%g Check Sieve (M1)=%g", M1+1, R[M1], N, M1 ); EndIf
                // exit loop:
                M1 = N; 
                Check_Flag = 1;
            EndIf
        EndFor
        If (!Check_Flag)
            Error("    Failed to find sieve of aggregate %g : Current Sieve (N)=%g M1=%g", check_Target+1, N, M1 );
            Abort;
        EndIf
        // loop previous sieves
        For M1 In {0:N-1}
            For j1 In {0:numAgg[M1]-1}
                d1 = ((x[check_Target] - x[Check_C2])^2 + (y[check_Target] - y[Check_C2])^2 + (z[check_Target] - z[Check_C2])^2)^0.5;
                confLen1 = R[Check_TargetSieve] + R[M1] + aggwalleff - d1;
                If(Check_C2 == check_Target)
                    If(verbose) Printf("        OK: skip self R1: %g", R[Check_TargetSieve]); EndIf
                ElseIf(confLen1 > Tol)
                    Error("    CHECK FAILED ::: Aggregate %g [%g %g %g] (%g in the sieve) conflicts (%g mm) with target aggregate %g [%g %g %g] ", 
                            Check_C2+1, x[Check_C2], y[Check_C2], z[Check_C2], j1+1, confLen1, 
                            check_Target+1, x[check_Target], y[check_Target], z[check_Target]);
                        Printf("    R1: %g    R2: %g   R1+R2: %g", R[Check_TargetSieve], R[M1], R[Check_TargetSieve] + R[M1]);
                        Printf("    distance: %g    d-R1+R2: %g", d1, d - R[Check_TargetSieve] + R[M1]);
                    Abort;
                Else
                    If(verbose) Printf("        OK: Aggregate %g [%g %g %g] conflict %g mm target aggregate %g [%g %g %g] R1: %g R2: %g", 
                            Check_C2+1, x[Check_C2], y[Check_C2], z[Check_C2], confLen1, 
                            check_Target+1, x[check_Target], y[check_Target], z[check_Target], R[Check_TargetSieve], R[M1]); EndIf
                EndIf
                Check_C2++;
            EndFor
        EndFor
        //M1++;
        If(verbose) Printf("    M1++ : Current Sieve (N)=%g M1=%g", N, M1 ); EndIf
        For j1 In {0:i-1}
            d1 = ((x[check_Target] - x[Check_C2])^2 + (y[check_Target] - y[Check_C2])^2 + (z[check_Target] - z[Check_C2])^2)^0.5;
            confLen1 = R[Check_TargetSieve] + R[M1] + aggwalleff - d1; // both are from current sieve N
            If(Check_C2 == check_Target)
                    If(verbose) Printf("        OK: skip self"); EndIf
            ElseIf(confLen1 > Tol)
                    Error("    REPORT FAILED ::: Aggregate %g [%g %g %g] (%g in the sieve) conflicts (%g mm) with target aggregate %g [%g %g %g] ", 
                            Check_C2+1, x[Check_C2], y[Check_C2], z[Check_C2], j1+1, confLen1, 
                            check_Target+1, x[check_Target], y[check_Target], z[check_Target]);
                    Printf("    R1: %g    R2: %g   R +R : %g", R[Check_TargetSieve], R[M1], R[Check_TargetSieve] + R[M1]);
                    Printf("    distance: %g     d-R +R : %g", d1,     d1 - R[Check_TargetSieve] + R[M1]);
                    Abort;
            Else
                If(verbose) Printf("        OK: Aggregate %g [%g %g %g] conflict %g mm target aggregate %g [%g %g %g] R1: %g R2: %g", 
                            Check_C2+1, x[Check_C2], y[Check_C2], z[Check_C2], confLen1, 
                            check_Target+1, x[check_Target], y[check_Target], z[check_Target], R[Check_TargetSieve], R[M1]); EndIf
            EndIf
            Check_C2++;
        EndFor
        //Printf("    >>> End of Macro <<<    ");
Return
/// -----------------------------------------
///
/// ----------------------------------------------
Macro Move_Macro
        Printf("    .: Move Macro :. " );
        Printf("    move target aggregate %g ", move_Target+1 );

        hits[move_Target] = 0;
        move_flagConf = 0;
        confLen0 = 0;
        C02 = 0;
        // also see move_tryMax = 25;
        move_try = 0;
        
        X0 = x[move_Target];
        Y0 = y[move_Target];
        Z0 = z[move_Target];

        x[move_Target] += moveDirection[0]; // shift x
        y[move_Target] += moveDirection[1]; // shift y
        z[move_Target] += moveDirection[2]; // shift z
   
        /// find radii of target
        Move_TargetSieve = 0;
        temp = -1; // it's 0-based
        Move_Flag = 0;
        For M0 In {0:N} //sieves
            temp += numAgg[M0]; // last aggregate in sieve M
            If(verbose) Printf("    N=%g M0=%g Aggregates in Check Sieve=%g temp=%g, move_Target=%g", N, M0, numAgg[M0], temp, move_Target); EndIf
            If(temp >= move_Target)
                Move_TargetSieve = M0;
                If(verbose) Printf("    Move Target Sieve #: %g, R: %g mm Current-Last Sieve (N)=%g Check Sieve (M0)=%g", M0+1, R[M0], N, M0 ); EndIf
                // exit loop:
                M0 = N; 
                Move_Flag = 1;
            EndIf
        EndFor
        If (!Move_Flag)
            Error("    Failed to find sieve of aggregate %g : Current Sieve (N)=%g M1=%g", move_Target+1, N, M1 );
            Abort;
        EndIf

        // Keep Macro paramters:
        keep_radii  = R[Move_TargetSieve];
        keep_target = move_Target;
        // check specimen boundaries:
        If (verbose) Printf ("@@@"); EndIf
        Call Keep_Macro;

        // loop previous sieves
        For M0 In {0:N-1} // sieve number
            If (verbose) Printf("    +++ "); EndIf
            If (M0 < 0) Error ("M0 is smaller than zero"); Abort; EndIf
            //If (verbose) Printf("   Start Loop M0 = %g", M0); EndIf
            If (verbose) Printf("    Checking sieve %g of %g with %g aggregates, Check target aggregate %g ", M0+1, N+1, numAgg[M0], move_Target+1); EndIf

            For j0 In {0:numAgg[M0]-1} // aggregate in sieve
                If(C02 == move_Target)
                    If (verbose) Printf("    Aggregate %g (Self) (%g/%g sieve %g)", C02+1, j0+1, numAgg[M0], M0+1); EndIf
                    // skip self
                    C02++;
                    confLen0 = 0;
                Else
                    If (verbose) Printf("    Aggregate %g (%g/%g sieve %g)", C02+1, j0+1, numAgg[M0], M0+1); EndIf
                    d0 = ((x[move_Target] - x[C02])^2 + (y[move_Target] - y[C02])^2 + (z[move_Target] - z[C02])^2)^0.5;
                    confLen0 = R[Move_TargetSieve] + R[M0] + aggwalleff - d0;
                    If(confLen0 > Tol)
                        If (debug) Printf("    +++ "); EndIf
                        If (debug) Printf("    Aggregate %g (%g sieve %g) conflicted with move target %g", C02+1, j0+1, M0+1, move_Target+1); EndIf
                        move_try++;
                        If(move_try > move_tryMax)
                            Printf("    Move try exceeds");
                            move_flagConf = 1;
                            // exit both loops
                            j0 = numAgg[M0] - 1;
                            M0 = N - 1;
                        Else
                            If (verbose) Printf("    shift coordinate and reset C02=0"); EndIf
                            If (verbose) Printf("    "); EndIf
        
                            x[move_Target] += ( x[move_Target] - x[C02] )/ d0 * confLen0; // shift x
                            y[move_Target] += ( y[move_Target] - y[C02] )/ d0 * confLen0; // shift y
                            z[move_Target] += ( z[move_Target] - z[C02] )/ d0 * confLen0; // shift z

                            // check specimen boundaries:
                            If (verbose) Printf ("    @@@"); EndIf
                            Call Keep_Macro;

                            // reset check aggregate ID
                            C02 = 0;
                            // exit current loop
                            j0 = numAgg[M0] - 1;
                            // reset outer loop
                            M0 = -1;
                        EndIf
                    Else
                        If (verbose) Printf("    Aggregate %g no conflict, M0=%g increase C02 ++ %g", C02+1, M0, C02); EndIf
                        C02++;
                        confLen0 = 0;
                    EndIf
                EndIf
            EndFor
            //  loop in current sieve (N) 
            If( move_flagConf == 0 && M0 == N-1)
                M0++; // make M0 = N which is current sieve
                If (verbose) Printf("    +++ "); EndIf
                If (verbose) Printf("    Checking sieve %g (current sieve)", N+1); EndIf
                If (verbose) Printf("  last aggregate of current sieve  i = %g  ", i); EndIf
                For j0 In {0:i-1}
                    If (verbose) Printf("    Aggregate %g (%g/%g of current sieve)  ", C02+1, j0+1, i); EndIf
                    d0 = ((x[move_Target] - x[C02])^2 + (y[move_Target] - y[C02])^2 + (z[move_Target] - z[C02])^2)^0.5;
                    confLen0 = R[Move_TargetSieve] + R[N] + aggwalleff - d0;
                    If(confLen0 > Tol && C02 != move_Target)
                        If (debug) Printf("    +++ "); EndIf
                        If (debug) Printf("    Aggregate %g (%g current sieve) conflicted with move target %g", C02+1, j0+1, move_Target+1); EndIf
                        move_try++;
                        If(move_try > move_tryMax)
                                            Printf("    move try exceed");

                            move_flagConf = 1;
                            // exit both loop
                            j0 = i-1;
                            M0 = N - 1;
                        Else
                            If (verbose) Printf("    reset C02 and loop"); EndIf
                            x[move_Target] += ( x[move_Target] - x[C02] )/ d0 * confLen0; // shift x
                            y[move_Target] += ( y[move_Target] - y[C02] )/ d0 * confLen0; // shift y
                            z[move_Target] += ( z[move_Target] - z[C02] )/ d0 * confLen0; // shift z

                            // check specimen boundaries:
                            If (verbose) Printf("    @@@ "); EndIf
                            Call Keep_Macro;

                            // reset check aggregate ID
                            C02 = 0;
                            // reset outer loop
                            M0 = -1;
                            // exit current loop
                            j0 = i-1;
                        EndIf
                    Else
                        If (verbose) Printf("    no conflict, increase C02=%g+1", C02); EndIf
                        C02++;
                        confLen0 = 0;
                    EndIf
                EndFor
            Else
                If (move_flagConf) 
                    //Printf("    Resetting conflict flag");
                    //move_flagConf = 0;
                    Printf("    Conflict flag is true");
                EndIf
            EndIf
        EndFor
        // M0 should be N+1:
        If (verbose) Printf("End Loop M0 = %g, N = %g", M0, N); EndIf
        // check move success:
        If ( move_flagConf )
            // move failed, revert
            x[move_Target] = X0; // cancel shift x
            y[move_Target] = Y0; // cancel shift y
            z[move_Target] = Z0; // cancel shift z
            Printf("    Moving conflicted with %g and conf.Len of %g", C02+1, confLen0);
            Printf("    Original position: [%g %g %g] Target %g", x[move_Target], y[move_Target], z[move_Target], move_Target+1 );
            Printf("    Moving  Direction: [%g %g %g] ", moveDirection[0], moveDirection[1], moveDirection[2]);
            Printf("    Conflict position: [%g %g %g] of conflict aggregate %g ", x[C02], y[C02], z[C02], C02 + 1 );

            hits[C02]++;
        Else
            Printf("    Moving Succeed!!!   ");

            If (conservative) 
                Printf(" ::: check after move");
                Printf("    *** Start...");
                //check_Target = C1;
                check_Target = move_Target;
                Call Check_Macro;
                // check boundary
                keep_target = move_Target;
                keep_radii = R[Move_TargetSieve];
                keep_error = 1;
                Call Keep_Macro;
                keep_error = 0;
                Printf("    *** End; Succeed!");
            EndIf



        EndIf
            Printf("    >>> End of Move Macro <<<    ");
Return
/// -----------------------------------------
///
numAggTotal = 0;
For N In {0:sieveNum-1}
    numAggTotal += numAgg[N];
EndFor
Printf("Starting random distribution of %g total aggregates...", numAggTotal);
///
// Generate Aggregates
// initiate coordinates
x(numAggTotal) = 0;
y(numAggTotal) = 0;
z(numAggTotal) = 0;
hits(numAggTotal) = 0;
///
k = 0;  // conflict flag 1
shiftCount = 0;         // for maximum shift tries
move_Target = 0;         // move target
moveDirection(3) = 0;   // moving vector
C1 = 0; // current aggregate count (base 0)
C2 = 0; // inner loop aggregate count
confLen = 0.; // conflict length
For N In {0:sieveNum-1}         // sieve loop: M & N
    flagCalculateAgg = 1;
    Printf("> Next sieve %g mm", D[N]);
    For i In {0:numAgg[N]-1}    // sieve aggregate loop: i & j
        //Printf(">> Aggregate %g", C1);
        k = 0;
        If ( confLen > 0 )
            keep_target = C1;
            keep_radii = R[N];
            Call Keep_Macro;
        EndIf
        If (confLen <= 0.)
            Printf("Assigning random coordinates to aggregate %g ...", C1+1);
            shiftCount = 0; 
            W = width  - R[N] * 2. - wallEffect * 2.; // left and right => -R x 2
            H = height - R[N] * 2. - wallEffect * 2.; // top and bottom => -R x 2
            D = depth  - R[N] * 2. - wallEffect * 2.; // front and back => -R x 2
            x[C1] = Rand(W) + R[N] + wallEffect; //  left shift => +R
            y[C1] = Rand(H) + R[N] + wallEffect; // bottom shift=> +R
            z[C1] = Rand(D) + R[N] + wallEffect; // front shift=> +R
        EndIf
        // check if location is available
        C2 = 0;
        //  loop in past sieves
        For M In {0:N-1}
            For j In {0:numAgg[M]-1}
                d = ((x[C1] - x[C2])^2 + (y[C1] - y[C2])^2 + (z[C1] - z[C2])^2)^0.5;
                confLen = R[N] + R[M] + aggwalleff - d;
                If(confLen > Tol)
                    factor = 1.;
                    hits[C2]++;
                    If(hits[C2] > moveMax)
                        factor = 0.5;
                        Printf("Aggregate %g made too many conflicts, calling move macro...", C2+1);
                        move_Target = C2;
                        moveDirection = {(x[C2] - x[C1]) / d * confLen,
                                         (y[C2] - y[C1]) / d * confLen, 
                                         (z[C2] - z[C1]) / d * confLen };
                        Call Move_Macro;
                        Printf("    " );
                    EndIf
                    // update xyz
                    x[C1] += (x[C1] - x[C2]) / d * confLen * factor; // shift x
                    y[C1] += (y[C1] - y[C2]) / d * confLen * factor; // shift y
                    z[C1] += (z[C1] - z[C2]) / d * confLen * factor; // shift z
                    // update flag and report
                    k++;
                    shiftCount++;
                    If(shiftCount >= shiftCountMax) 
                        Printf("Too many failures, replacing aggregate %g ...", C1+1);
                        confLen = -1.;
                    Else
                        If (debug) Printf("Fail: aggregate %g of %g conflicted with aggregate %g with conf.Len. %g:%.3g", C1+1, numAggTotal, C2+1, confLen, Tol); EndIf
                    EndIf
                    // exist loop j, M
                    j = numAgg[M] - 1;
                    M = N - 1;
                    // reset C2
                    C2 = 0;
                    // repeat loop i
                    i--;
                Else
                    confLen = 0;
                    C2++;
                EndIf
            EndFor
        EndFor
        //  loop in current sieve
        If(k==0)
            For j In {0:i-1}
                d = ((x[C1] - x[C2])^2 + (y[C1] - y[C2])^2 + (z[C1] - z[C2])^2)^0.5;
                If ( d ==0)
                    Printf("d = 0, C1 = %g, C2 = %g", C1, C2);
                EndIf
                confLen = R[N] + R[M] + aggwalleff - d;
                If(confLen > Tol)
                    If(d == 0)
                        Warning ("line 479~481 d length is zero, changed to 1e-6");
                        d = 0.000001;
                    EndIf
                    // update xyz
                    x[C1] += (x[C1] - x[C2]) / d * confLen; // shift x
                    y[C1] += (y[C1] - y[C2]) / d * confLen; // shift y
                    z[C1] += (z[C1] - z[C2]) / d * confLen; // shift z
                    // exiting loop
                    k++;
                    shiftCount++;
                    If(shiftCount >= shiftCountMax) 
                        Printf("Too many failures, replacing aggregate %g ...", C1+1);
                        confLen = -1.;
                    Else
                        If (debug) Printf("Fail: aggregate %g of %g conflicts with aggregate %g with conf.Len. %g:%.3g", C1+1, numAggTotal, C2+1, confLen, Tol); EndIf
                    EndIf
                    j = i - 1;
                    M = N; // no need?
                    i--;
                    C2 = 0;
                Else
                    confLen = 0;
                    C2++;
                EndIf
            EndFor
        EndIf
        If(k==0)
        //if no conflict found in past and current loops then create circle
            //Sphere(C1+1) = {x[C1], y[C1], z[C1], R[N], -Pi/2, Pi/2, 2*Pi};
            Printf("    -----------------   ");
            Printf("SUCCESS: aggregate %g (Agg(%g) in Sieve(%g)) was created", C1+1, i+1, N+1);
            Printf("         X %g , Y %g , Z %g ", x[C1], y[C1], z[C1]);
            Printf("    -----------------   ");
            Printf("       ");
            If (conservative) 
                Printf(" ::: check after success");
                Printf("    *** Start...");
                //check_Target = C1;
                For check_Target In {0:C1}
                    Call Check_Macro; 
                EndFor
                Printf("    *** End; Succeed!");
            EndIf
            Printf("       ");
            Printf("       ");
            Printf("       ");
            Printf("       ");
            Printf("       ");
            Printf("       ");
            Printf("       ");
            C1++;
            //Draw;
        EndIf
    EndFor
EndFor
///
Printf("    ::: FINAL CHECK ::: ");
For j In {0:C1-1}
    check_Target = j;
    N = sieveNum-1;
    i = numAgg[sieveNum-1]-1;
    If (conservative) Call Check_Macro; EndIf

    /// make aggregate if check passes:
    // find radii of target
    Target = j;
    TargetSieve = 0;
    temp = -1; // it's 0-based
    Flag = 0;
    For M In {0:N} //sieves
        temp += numAgg[M]; // last aggregate in sieve M
        If(verbose) Printf("    N=%g M0=%g Aggregates in Check Sieve=%g temp=%g, Target=%g", N, M, numAgg[M], temp, Target); EndIf
        If(temp >= Target)
            TargetSieve = M;
            If(debug) Printf("    Aggregate %g/%g of Sieve %g/%g, R: %g mm ", j+1, numAggTotal, M+1, N+1, R[M] ); EndIf
            // exit loop:
            M = N; 
            Flag = 1;
        EndIf
    EndFor
    If (!Flag)
        Error("    Failed to find sieve of aggregate %g : Sieve %g/%g", Target+1, N+1, M+1 );
        Abort;
    EndIf
    Sphere(j+1) = {x[j], y[j], z[j], R[TargetSieve], -Pi/2, Pi/2, 2*Pi};

    Sleep 0.02;
    Draw;
EndFor
///
C2 = C1++; // last aggregate global number
Printf("Last aggregate global number =%g", C2);
///
// .: Make volume :.
// make specimen surfaces
// - back
Curve Loop(C1) = {3, 4, 1, 2};
Plane Surface(C1) = {C1};
Physical Surface("specimen_surface_bak", 11) = {C1++};
// - top
Curve Loop(C1) = {11, 7, -12, -3};
Plane Surface(C1) = {C1};
Physical Surface("specimen_surface_top", 12) = {C1++};
// - right
Curve Loop(C1) = {10, 6, -11, -2};
Plane Surface(C1) = {C1};
Physical Surface("specimen_surface_rit", 13) = {C1++};
// - bottom
Curve Loop(C1) = {10, -5, -9, 1};
Plane Surface(C1) = {C1};
Physical Surface("specimen_surface_bot", 14) = {C1++};
// - front
Curve Loop(C1) = {5, 6, 7, 8};
Plane Surface(C1) = {C1};
Physical Surface("specimen_surface_frn", 15) = {C1++};
// - left
Curve Loop(C1) = {9, -8, -12, 4};
Plane Surface(C1) = {C1};
Physical Surface("specimen_surface_lft", 16) = {C1++};
//
///
C3 = C1; // Start of surface loop IDs
///
// Specimen's loop
Printf("C1=%g", C1);
Printf("C2=%g", C2);
Printf("C3=%g", C3);
Surface Loop(C1++) = {C2+1:C3-1};
///
// Aggregates' Surface loop
For i In {1:C2}
    Surface Loop(C1++) = {i};
    //Physical Volume("Aggregates", 2) += {i};
    //Physical Volume("Aggregates", 2) += {C1-1};
EndFor
// Volume = specimen surface loop : last aggregate surface loop
Volume(C2+1) = {C3:C1-1};
// make physical groups
Physical Volume("Matrix", 1) = {C2+1};
Physical Volume("Aggregates", 2) = {1:C2};
// Set mesh size of all geometry
If(userMeshSizeofAllGeometry>0)
    MeshSize{ PointsOf{ Volume{:}; } } = userMeshSizeofAllGeometry;
EndIf
For i In {1:iPnt-1}
    MeshSize {i} = h;
EndFor
//
Printf("");
Printf("NUMBER OF AGGREGATES: %g", numAggTotal);
Printf("");
Printf("SCRIPT FINISHED SUCCESSFULLY");
//
// ===================================================
// MAKING ATTACHED PLATES I & II
// ===================================================
Physical Volume("SteelCaps", 3) = {};
// Bottom Plate I Points
PL1_Bot_Points(4) = 0;
PL1_Bot_Points = Translate {0, -10, 0} { Duplicata { Point{1}; Point{2}; Point{5}; Point{6}; } };
//For i In {1:4}
//    Printf("%g", PL1_Bot_Points[i-1]);
//EndFor
//45:46
Translate {40, -10, 0} {
  Duplicata { Point{1}; Point{5}; }
}
Translate {-40, -10, 0} {
  Duplicata { Point{2}; Point{6}; }
}
// Bottom Plate II
/*
Translate {40, -20, 0} {
  Duplicata { Point{1}; Point{5}; }
}
Translate {-40, -20, 0} {
  Duplicata { Point{2}; Point{6}; }
}
*/
// Bottom Plate I Lines
//
PL1_Bot_corner_lines(4) = 0;
// temp: 61
temp = newl;
Line(newl) = {1, PL1_Bot_Points[0]};
Line(newl) = {2, PL1_Bot_Points[1]};
Line(newl) = {5, PL1_Bot_Points[2]};
Line(newl) = {6, PL1_Bot_Points[3]};
//
//65
Line(newl) = {PL1_Bot_Points[0], PL1_Bot_Points[3]+1};
//66
Line(newl) = {PL1_Bot_Points[3]+1, PL1_Bot_Points[3]+3};
//67
Line(newl) = {PL1_Bot_Points[3]+3, PL1_Bot_Points[1]};
//68
Line(newl) = {PL1_Bot_Points[1], PL1_Bot_Points[3]};
//69
Line(newl) = {PL1_Bot_Points[3], PL1_Bot_Points[3]+4};
//70
Line(newl) = {PL1_Bot_Points[3]+4, PL1_Bot_Points[3]+2};
//71
Line(newl) = {PL1_Bot_Points[3]+2, PL1_Bot_Points[2]};
//72
Line(newl) = {PL1_Bot_Points[2], PL1_Bot_Points[0]};
// ***
Line(newl) = {PL1_Bot_Points[3]+1, PL1_Bot_Points[3]+2};
Line(newl) = {PL1_Bot_Points[3]+4, PL1_Bot_Points[3]+3};
// ***
// temp_cl/s: 23
temp_cl = newcl;
temp1 = news;
Curve Loop(temp_cl) = {1, temp+1, -6-temp, -temp-5, -temp-4, -temp};
Plane Surface(temp1) = {temp_cl};
temp_cl=newcl;
temp2 = news;
Curve Loop(temp_cl) = {temp+1, temp+7, -temp-3, -10};
Plane Surface(temp2) = {temp_cl};
temp_cl=newcl;
temp3 = news;
Curve Loop(temp_cl) = {5, temp+3, temp+8, temp+9, temp+10, -temp-2};
Plane Surface(temp3) = {temp_cl};
temp_cl=newcl;
temp4 = news;
Curve Loop(temp_cl) = {9, temp+2, temp+11, -temp};
Plane Surface(temp4) = {temp_cl};
temp_cl=newcl;
temp5 = news;
Curve Loop(temp_cl) = {temp+4, temp+12, temp+10, temp+11};
Plane Surface(temp5) = {temp_cl};
// ***
temp_cl=newcl;
temp6 = news;
Curve Loop(temp_cl) = {temp+5, -temp-13, temp+9, -temp-12};
Plane Surface(temp6) = {temp_cl};
temp_cl=newcl;
temp7 = news;
Curve Loop(temp_cl) = {temp+6, temp+7, temp+8, temp+13};
Plane Surface(temp7) = {temp_cl};
// ***
//
PL2_Extrude_Bot = temp6;
temp = newsl;
volNum = newv;
Surface Loop(temp) = {temp1, temp2, temp7, temp3, temp6, temp5, temp4, numAggTotal+4};
Volume(volNum) = {temp};
Physical Surface("plateI_surface_bot", 17) = {temp5, temp6, temp7};
Physical Volume("SteelCaps", 3) += {volNum};
//
// ------------------------------------------
// Top Plate I Points
PL1_Top_Points(4) = 0;
PL1_Top_Points = Translate {0, 10, 0} { Duplicata { Point{3}; Point{4}; Point{7}; Point{8}; }};
Translate {-40, 10, 0} {
  Duplicata { Point{3}; Point{7}; }
}
Translate {40, 10, 0} {
  Duplicata { Point{4}; Point{8}; }
}
// Top Plate II
/*
Translate {-40, 20, 0} {
  Duplicata { Point{3}; Point{7}; }
}
Translate {40, 20, 0} {
  Duplicata { Point{4}; Point{8}; }
}
*/
// Top Plate I Lines
//
PL1_Top_corner_lines(4) = 0;
temp = newl;
Line(temp) = {3, PL1_Top_Points[0]};
Line(newl) = {4, PL1_Top_Points[1]};
Line(newl) = {7, PL1_Top_Points[2]};
Line(newl) = {8, PL1_Top_Points[3]};
//
Line(newl) = {PL1_Top_Points[0], PL1_Top_Points[3]+1};
Line(newl) = {PL1_Top_Points[3]+1, PL1_Top_Points[3]+3};
Line(newl) = {PL1_Top_Points[3]+3, PL1_Top_Points[1]};
Line(newl) = {PL1_Top_Points[1], PL1_Top_Points[3]};
Line(newl) = {PL1_Top_Points[3], PL1_Top_Points[3]+4};
Line(newl) = {PL1_Top_Points[3]+4, PL1_Top_Points[3]+2};
Line(newl) = {PL1_Top_Points[3]+2, PL1_Top_Points[2]};
Line(newl) = {PL1_Top_Points[2], PL1_Top_Points[0]};
// ***
Line(newl) = {PL1_Top_Points[3]+1, PL1_Top_Points[3]+2};
Line(newl) = {PL1_Top_Points[3]+4, PL1_Top_Points[3]+3};
// ***
//
temp_cl = newcl;
temp1 = news;
Curve Loop(temp_cl) = {3, temp+1, -6-temp, -temp-5, -temp-4, -temp};
Plane Surface(temp1) = {temp_cl};
temp_cl=newcl;
temp2 = news;
Curve Loop(temp_cl) = {temp+1, temp+7, -temp-3, -12};
Plane Surface(temp2) = {temp_cl};
temp_cl=newcl;
temp3 = news;
Curve Loop(temp_cl) = {7, temp+3, temp+8, temp+9, temp+10, -temp-2};
Plane Surface(temp3) = {temp_cl};
temp_cl=newcl;
temp4 = news;
Curve Loop(temp_cl) = {11, temp+2, temp+11, -temp};
Plane Surface(temp4) = {temp_cl};
temp_cl=newcl;
temp5 = news;
Curve Loop(temp_cl) = {temp+4, temp+12, temp+10, temp+11};
Plane Surface(temp5) = {temp_cl};
// ***
temp_cl=newcl;
temp6 = news;
Curve Loop(temp_cl) = {temp+5, -temp-13, temp+9, -temp-12};
Plane Surface(temp6) = {temp_cl};
temp_cl=newcl;
temp7 = news;
Curve Loop(temp_cl) = {temp+6, temp+7, temp+8, temp+13};
Plane Surface(temp7) = {temp_cl};
// ***
//
PL2_Extrude_Top = temp6;
temp = newsl;
Surface Loop(temp) = {temp1, temp2, temp7, temp3, temp6, temp5, temp4, numAggTotal+2};
volNum = newv;
Volume(volNum) = {temp};
Physical Surface("plateI_surface_top", 18) = {temp5, temp6, temp7};
Physical Volume("SteelCaps", 3) += {volNum};
//
// ------------------------------------------------------------------------------------
// Plate II
out[] = Extrude {0, -10, 0} {
  Surface{PL2_Extrude_Bot}; 
};
// out[0] is the new surface
Physical Surface("plateII_surface_bot", 19) = {out[0]};
Physical Volume("SteelCaps", 3) += {out[1]};
//  - - - - - - - - - - - - - - - - - - - - - 
out[] = Extrude {0, 10, 0} {
  Surface{PL2_Extrude_Top}; 
};
// out[0] is the new surface
Physical Surface("plateII_surface_top", 20) = {out[0]};
Physical Volume("SteelCaps", 3) += {out[1]};
//
// =============================================

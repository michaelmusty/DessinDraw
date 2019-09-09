174,0
V,BelyiNewton,1
S,GetAssignedAttributes,return assigned attributes of Gamma in a sequence,0,1,0,0,0,0,0,0,0,GrpPSL2,,82,-38,-38,-38,-38,-38
S,NewtonGenusOne,,1,0,1,82,0,GrpPSL2,1,0,0,0,0,0,0,0,82,,GrpPSL2,-38,-38,-38,-38,-38
S,NewtonGenusOne,Less Naive wrapper..,0,1,0,0,0,0,0,0,0,GrpPSL2,,GrpPSL2,-38,-38,-38,-38,-38
S,NewtonGetNumericalData,"Computes numerical data necessary for Newton, writes it to Gamma and returns Gamma",0,1,0,0,0,0,0,0,0,GrpPSL2,,GrpPSL2,-38,-38,-38,-38,-38
S,NewtonGetRamificationPoints,"Assigns TriangleNewtonRamificationPoints0,1,oo to Gamma, a list of pairs [x_p,y_p] (for each of 0,1,oo) on the curve over CC",0,1,0,0,0,0,0,0,0,GrpPSL2,,GrpPSL2,-38,-38,-38,-38,-38
S,TwoTorsionTest,"Given a point w in the hyperbolic disc associated to Gamma, test if w is a 2-torsion point of the associated elliptic curve",0,3,0,0,0,0,0,0,0,82,,0,0,GrpPSL2,,0,0,SpcHydElt,,-1,-38,-38,-38,-38,-38
S,NewtonVanishingEquations,,0,1,0,0,0,0,0,0,0,GrpPSL2,,-1,-38,-38,-38,-38,-38
S,NewtonGetBasicEquations,"Computes basic Newton equations (ramification, order of vanishing) and assigns them to Gamma",0,1,0,0,0,0,0,0,0,GrpPSL2,,GrpPSL2,-38,-38,-38,-38,-38
S,NewtonGetBasicInitializationValues,"Assigns start_vector [c4, c6, points0, points1, pointsoo, extra_points, lc, num_coeffs, den_coeffs] to Gamma",0,1,0,0,0,0,0,0,0,GrpPSL2,,GrpPSL2,-38,-38,-38,-38,-38
S,NewtonGetRescalingEquation,assign (polynomial equation for rescaling) to Gamma,0,1,0,0,0,0,0,0,0,GrpPSL2,,GrpPSL2,-38,-38,-38,-38,-38
S,NewtonIterate,Newton iterate starting solution to equations (polynomials) to get a solution to precision precNewton,2,0,1,82,0,63,1,1,82,0,172,3,0,0,0,0,0,0,0,148,,0,0,82,,0,0,82,,82,-38,-38,-38,-38,-38
S,NewtonIterate,uses equations and initial values assigned to Gamma,0,2,0,0,0,0,0,0,0,148,,0,0,GrpPSL2,,GrpPSL2,-38,-38,-38,-38,-38
S,NewtonRecognize,Recognize elements of solution (complex numbers) with power relations up to bound,0,1,0,0,0,0,0,0,0,GrpPSL2,,GrpPSL2,-38,-38,-38,-38,-38
S,NewtonMakeBelyiMaps,Assigns Belyi curve and Belyi map to Gamma after some sanity checks,0,1,0,0,0,0,0,0,0,GrpPSL2,,GrpPSL2,-38,-38,-38,-38,-38

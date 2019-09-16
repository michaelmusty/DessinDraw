174,0
S,GetAssignedAttributes,return assigned attributes of Gamma in a sequence,0,1,0,0,0,0,0,0,0,519,,82,-38,-38,-38,-38,-38
S,NewtonHyperelliptic,wrapper..,0,1,0,0,0,0,0,0,0,519,,519,-38,-38,-38,-38,-38
S,NewtonHyperellipticGetNumericalData,"Computes numerical data necessary for Newton, writes it to Gamma and returns Gamma",0,1,0,0,0,0,0,0,0,519,,519,-38,-38,-38,-38,-38
S,NewtonHyperellipticGetRamificationPoints,"Assigns TriangleNewtonRamificationPoints0,1,oo to Gamma, a list of pairs [x_p,y_p] (for each of 0,1,oo) on the curve over CC",0,1,0,0,0,0,0,0,0,519,,519,-38,-38,-38,-38,-38
S,HyperellipticTwoTorsionTest,,0,3,0,0,0,0,0,0,0,82,,0,0,519,,0,0,521,,-1,-38,-38,-38,-38,-38
S,PolarPart,,0,1,0,0,0,0,0,0,0,285,,285,-38,-38,-38,-38,-38
S,RiemannRochBasisHyperellipticSimple,,0,2,0,0,0,0,0,0,0,519,,0,0,148,,-1,-38,-38,-38,-38,-38
S,RiemannRochBasisHyperellipticSimpleAnalytic,,0,2,0,0,0,0,0,0,0,519,,0,0,148,,-1,-38,-38,-38,-38,-38
S,RiemannRochBasisHyperellipticFormal,Basis for L(m*(infinity_1+infinity_2)) as function field elements,0,3,0,0,0,0,0,0,0,82,,0,0,519,,0,0,148,,-1,-38,-38,-38,-38,-38
S,NewtonHyperellipticGenericBelyiMap,Make generic hyperelliptic Belyi map,0,1,0,0,0,0,0,0,0,519,,-1,-38,-38,-38,-38,-38
S,NewtonHyperellipticEchelonizationEquations,,0,1,0,0,0,0,0,0,0,519,,-1,-38,-38,-38,-38,-38
S,NewtonHyperellipticVanishingEquations,,0,1,0,0,0,0,0,0,0,519,,-1,-38,-38,-38,-38,-38
S,NewtonHyperellipticSpecialEquations,,0,1,0,0,0,0,0,0,0,519,,-1,-38,-38,-38,-38,-38
S,NewtonHyperellipticGetBasicEquations,"Computes basic Newton equations (ramification, order of vanishing) and assigns them to Gamma",0,1,0,0,0,0,0,0,0,519,,519,-38,-38,-38,-38,-38
S,NewtonHyperellipticNumericalBelyiMap,,0,1,0,0,0,0,0,0,0,519,,-1,-38,-38,-38,-38,-38
S,NewtonHyperellipticCommonZeroes,,0,1,0,0,0,0,0,0,0,519,,-1,-38,-38,-38,-38,-38
S,NewtonHyperellipticGetBasicInitializationValues,"Assigns start_vector [curve_coeffs, points0, points1, pointsoo, extra_points, lc, num_coeffs, den_coeffs, special_points] to Gamma",0,1,0,0,0,0,0,0,0,519,,519,-38,-38,-38,-38,-38
S,NewtonHyperellipticGetRescalingEquation,assign (polynomial equation for rescaling) to Gamma,0,1,0,0,0,0,0,0,0,519,,519,-38,-38,-38,-38,-38
S,NewtonIterate,Newton iterate starting solution to equations (polynomials) to get a solution to precision precNewton,2,0,1,82,0,63,1,1,82,0,172,3,0,0,0,0,0,0,0,148,,0,0,82,,0,0,82,,82,-38,-38,-38,-38,-38
S,NewtonIterate,uses equations and initial values assigned to Gamma,0,2,0,0,0,0,0,0,0,148,,0,0,519,,519,-38,-38,-38,-38,-38
S,NewtonHyperellipticRecognize,Recognize elements of solution (complex numbers) with power relations up to bound,0,1,0,0,0,0,0,0,0,519,,519,-38,-38,-38,-38,-38
S,NewtonHyperellipticMakeBelyiMaps,Assigns Belyi curve and Belyi map to Gamma after some sanity checks,0,1,0,0,0,0,0,0,0,519,,519,-38,-38,-38,-38,-38

174,0
S,TriangleRescaleCoefficients,"Given a GrpPSL2 Gamma, and SeqEnums coeffs and vals, return the rescaling factor lamba, as well as the rescaled coefficients. coeffs should be of the form [curve_coeffs, num_coeffs, denom_coeffs], and vals [curve_vals, num_vals, denom_vals]",0,3,0,0,0,0,0,0,0,82,,0,0,82,,0,0,GrpPSL2,,-1,-38,-38,-38,-38,-38
S,TriangleRecordCoefficients,"Given a sequence of triangle subgroups (GrpPSL2s), write FldComElt invariants and coefficients for the curve and BelyiMap to the objects",1,0,1,82,0,GrpPSL2,1,0,0,0,0,0,0,0,82,,-1,-38,-38,-38,-38,-38
S,TriangleRecordCoefficients,"Given a sequence of permutation triples, write FldComElt invariants and coefficients for the curve and BelyiMap to the objects",1,0,1,82,1,82,0,222,1,0,0,0,0,0,0,0,82,,-1,-38,-38,-38,-38,-38
S,TriangleRootMatcher5000,"Given a complex number z, thought to be a root of one of polynomials in fact_list (in the format output by Factorization()), identify which polynomial it satisfies. Returns the polynomial f, the error err, an optimized representation of the number field defined by f, and the embedding data given by a place v and boolean conj",0,2,0,0,0,0,0,0,0,82,,0,0,172,,-1,-38,-38,-38,-38,-38
S,GaloisMinPoly,"Given a list of Galois conjugate complex numbers returns their (possibly) reducible ""minpoly""",1,0,1,82,0,172,2,0,0,0,0,0,0,0,402,,0,0,82,,-1,-38,-38,-38,-38,-38
S,TriangleLengthSort,,1,0,1,82,0,GrpPSL2,1,0,0,0,0,0,0,0,82,,-1,-38,-38,-38,-38,-38
S,TriangleRecognizeCoefficients,"Given a sequence of triangle subgroups with all the numerical stuff computed, to each one assign a number field, coefficients of the curve and Belyi map as elements of this number field. ExactAl specifies algorithm by which coefficients are recognized. Options are GaloisOrbits (requires the sequence of Gammas to be closed under the action of Galois) and AlgebraicNumbers which attempts to recognize the algebraic coefficients directly over a number field with degree bounded by DegreeBound, if given (otherwise taken to be the size of the pointed passport for Gamma)",1,0,1,82,0,GrpPSL2,1,0,0,0,0,0,0,0,82,,-1,-38,-38,-38,-38,-38
S,TriangleRecognizeAlgebraicCoefficients,"Given a triangle subgroup with all the numerical stuff computed, assign it a number field, coefficients of the curve and Belyi map as elements of this number field",0,1,0,0,0,0,0,0,0,GrpPSL2,,-1,-38,-38,-38,-38,-38
S,TriangleMakeBelyiMap,"Given a GrpPSL2 with the number field K recognized and exact data computed construct the Belyi curve and Belyi map, then return and assign them",0,1,0,0,0,0,0,0,0,GrpPSL2,,-1,-38,-38,-38,-38,-38
S,TriangleMakeBelyiMaps,"Given a sequence of GrpPSL2 with the number field K recognized and exact data computed for each Gamma, construct the Belyi curves and Belyi maps",1,0,1,82,0,GrpPSL2,1,0,0,0,0,0,0,0,82,,-1,-38,-38,-38,-38,-38

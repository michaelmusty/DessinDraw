intrinsic HyperellipticDerivativeY(f::RngUPolElt, N::RngIntElt) -> Any
  {Given a polynomial f, returns the N^th derivative of y with respect to x as a function on the hyperelliptic curve y^2 = f(x).}
  R<x,y> := RationalFunctionField(BaseRing(Parent(f)), 2);
  h := hom<Parent(f)->R|[R.1]>;
  f := h(f);
  derivs := [y];
  for n := 1 to N do
    derivs_new := Derivative(f, 1, n);
    for k := 1 to n-1 do
      //derivs_new -:= &+[Binomial(n, k)*derivs[k+1]*derivs[n-k+1]];
      derivs_new -:= Binomial(n, k)*derivs[k+1]*derivs[n-k+1];
    end for;
    derivs_new *:= 1/(2*y);
    Append(~derivs, derivs_new);
  end for;
  return derivs[#derivs];
end intrinsic;

intrinsic OneVarTwoVar(f::FldFunRatUElt, R::FldFunRat) -> Any
  {Convert element of QQ(x)(y) to QQ(x,y)}
  f_num := Numerator(f);
  f_den := Denominator(f);
  num_cs := Coefficients(f_num);
  den_cs := Coefficients(f_den);
  assert Rank(R) eq 2;
  R<x,y> := R;
  num_cs := [Evaluate(el,x) : el in num_cs];
  den_cs:= [Evaluate(el,x) : el in den_cs];
  f_num_new := R!0;
  for i := 1 to #num_cs do
    f_num_new +:= num_cs[i]*y^(i-1);
  end for;
  f_den_new := R!0;
  for i := 1 to #den_cs do
    f_den_new +:= den_cs[i]*y^(i-1);
  end for;
  return f_num_new/f_den_new;
end intrinsic;

intrinsic OneVarTwoVarPoly(f::RngUPolElt, R::Rng) -> Any
  {Convert element of BFR(x)(y) to BFR[x,y]}
  f_num := Numerator(f);
  assert Denominator(f) eq 1;
  num_cs := Coefficients(f_num);
  R<x,y> := R;
  num_cs := [Evaluate(el,x) : el in num_cs];
  f_num_new := R!0;
  for i := 1 to #num_cs do
    f_num_new +:= num_cs[i]*y^(i-1);
  end for;
  return Numerator(f_num_new);
end intrinsic;

intrinsic HyperellipticDerivativePhi(Gamma::GrpPSL2, N::RngIntElt, k::RngIntElt) -> Any
  {Take the Nth derivative of the generic Belyi map of Gamma}
  assert k in {1,2,3};
  g := Genus(Gamma);
  curve_vars := Gamma`TriangleNewtonVariablesHyperellipticCurveCoefficients;
  vprint Shimura : "Constructing generic Belyi map";
  phi_num, phi_den, lc_var := NewtonHyperellipticGenericBelyiMap(Gamma);
  phi := phi_num/phi_den;
  R<x,y> := Parent(phi);
  S0<X> := RationalFunctionField(BaseRing(R));
  S<Y> := RationalFunctionField(S0);
  h := hom< R -> S | [S0.1, S.1] >;
  phi_num := Numerator(h(Numerator(phi)));
  phi_den := Numerator(h(Denominator(phi)));
  if k eq 1 then
    phi_root := phi_num;
  elif k eq 2 then
    phi_root := phi_num - phi_den;
  elif k eq 3 then
    phi_root := phi_den;
  else
    error "bad index";
  end if;
  // set up various coefficients
  vprint Shimura : "Setting up coefficients";
  a0 := Coefficient(phi_root,0);
  a1 := Coefficient(phi_root,1);
  // make curve poly
  // FIXME: include odd case
  v := S0!0;
  for i in [0..2*g] do
    v +:= curve_vars[i+1]*X^i;
  end for;
  v +:= X^(2*g+2);
  // write phi = a(x) + b(x)*y
  vprint Shimura : "Computing derivatives";
  phi_root_ds := [phi_root];
  for i := 1 to N do
    deriv_old := phi_root_ds[i];
    b0 := Coefficient(deriv_old,0);
    b1 := Coefficient(deriv_old,1);
    deriv_new := S!(Derivative(b0)) + S!(Derivative(b1))*h(y) + S!(b1)*h(HyperellipticDerivativeY(Numerator(v),1));
    Append(~phi_root_ds, Numerator(deriv_new)); // can we just take numerator here?
  end for;
  return phi_root_ds;
end intrinsic;

intrinsic EvaluateCoefficients(f::RngMPolElt, vals::SeqEnum) -> Any
  {Given a multivariable polynomial f whose base ring is also a polynomial ring, evaluate the coefficients of f at vals.}
  R := Parent(f);
  n := Rank(R);
  R0 := BaseRing(R);
  cs, mons := CoefficientsAndMonomials(f);
  cs_eval := [Evaluate(el, vals) : el in cs];
  //F := Parent(vals[1]);
  F := Parent(cs_eval[1]);
  S := PolynomialRing(F,n);
  if Rank(S) eq 2 then // hacky
    S<X,Y> := S;
  end if;
  mons_S := [S!mon : mon in mons];
  assert #mons_S eq #cs_eval;
  f_eval := S!0;
  for i := 1 to #cs_eval do
    f_eval +:= cs_eval[i]*mons_S[i];
  end for;
  return f_eval;
end intrinsic;

/*
intrinsic HyperellipticDerivativePhi(Gamma::GrpPSL2, N::RngIntElt) -> Any
  {Take the Nth derivative of the generic Belyi map of Gamma}
  g := Genus(Gamma);
  curve_vars := Gamma`TriangleNewtonVariablesHyperellipticCurveCoefficients;
  vprint Shimura : "Constructing generic Belyi map";
  phi_num, phi_den, lc_var := NewtonHyperellipticGenericBelyiMap(Gamma);
  phi := phi_num/phi_den;
  R<x,y> := Parent(phi);
  S0<X> := RationalFunctionField(BaseRing(R));
  S<Y> := RationalFunctionField(S0);
  h := hom< R -> S | [S0.1, S.1] >;
  phi1_num := Numerator(h(Numerator(phi)));
  phi1_den := Numerator(h(Denominator(phi)));
  // set up various coefficients
  vprint Shimura : "Setting up coefficients";
  a0 := Coefficient(phi1_num,0);
  a1 := Coefficient(phi1_num,1);
  b0 := Coefficient(phi1_den,0);
  b1 := Coefficient(phi1_den,1);
  // make curve poly
  // FIXME: include odd case
  v := S0!0;
  for i in [0..2*g] do
    v +:= curve_vars[i+1]*X^i;
  end for;
  v +:= X^(2*g+2);
  // write phi = a(x) + b(x)*y
  vprint Shimura : "Writing phi = a(x) + b(x)*y";
  a := (a0*b0 - a1*b1*v)/(b0^2 - b1^2*v);
  b := (-a0*b1 + a1*b0)/(b0^2 - b1^2*v);
  // now take the derivative
  vprint Shimura : "Calculating derivative";
  a_deriv := Derivative(a,N); // this is REALLY slow
  b_deriv := Derivative(b,N)*Y + b*h(HyperellipticDerivativeY(Numerator(v),N)); // product rule
  phi1_deriv := a_deriv + b_deriv;
  phi_deriv := OneVarTwoVar(phi1_deriv, R);
  phi_deriv *:= lc_var; // don't forget the lc_var!!!
  return phi_deriv;
end intrinsic;
*/

/*
intrinsic HyperellipticDerivativePhi(Gamma::GrpPSL2, N::RngIntElt, D_ind::RngIntElt) -> Any
  {}
  g := Genus(Gamma);
  t0 := Explode(Gamma`TriangleRiemannRochParameters);
  curve_vars := Gamma`TriangleNewtonVariablesHyperellipticCurveCoefficients;
  R0 := Parent(curve_vars[1]);
  R0frac := FieldOfFractions(R0);
  names := Names(R0);
  AssignNames(~R0frac,names);
  S<t> := PolynomialRing(R0);
  //h := S!0;
  v := S!0;
  for i in [0..2*g] do
    v +:= curve_vars[i+1]*t^i;
  end for;
  v +:= t^(2*g+2);
  vprint Shimura: "Forming basis for Riemann-Roch space";
  lc_var := Gamma`TriangleNewtonVariablesLeadingCoefficient;
  num_vars := Gamma`TriangleNewtonVariablesNumeratorCoefficients;
  den_vars := Gamma`TriangleNewtonVariablesDenominatorCoefficients;

  RR_basis := RiemannRochBasisHyperellipticSimple(t0, Gamma);
  assert #RR_basis ge #num_vars;
  assert #RR_basis ge #den_vars;
  R<x,y> := Parent(RR_basis[1]);
  R1<X> := PolynomialRing(R0);
  R2<Y> := PolynomialRing(R1);
  poly_map := hom< R -> R2 | 
end intrinsic;
*/

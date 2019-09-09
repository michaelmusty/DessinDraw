// ================================================================================
// Hyperelliptic
// ================================================================================

RemoveLeadingZeros := function(f,eps);
  fes := AbsEltseq(f);
  sf := Degree(LeadingTerm(f));
  while Abs(Coefficient(f,sf)) lt eps do
    fes[sf+1] := 0;
    sf +:= 1;
  end while;
  return Parent(f)!fes, sf;
end function;

RemoveLeadingZerosLaurent := function(f,eps);
  // returns error for machine zero series: leading term reaches precision threshold
  lc := LeadingCoefficient(f);
  while Abs(lc) lt eps do
    f := f - LeadingTerm(f);
    lc := LeadingCoefficient(f);
  end while;
  return f;
end function;

function ZeroifyCoeffs(coeffs, eps)
  assert Type(coeffs) eq SeqEnum;
  //assert #coeffs gt 0;
  CC := Parent(coeffs[1]);
  prec := Precision(CC);
  coeffs0 := [];
  for i in [1..#coeffs] do
    if Abs(coeffs[i]) lt eps then
      Append(~coeffs0, CC!0);
    else
      Append(~coeffs0, CC!coeffs[i]);
    end if;
  end for;
  return coeffs0;
end function;

function PolarPart(f)
  // Input: A Laurent series f
  // Output: The polar part of f, i.e., the Laurent tail of f
  new := f + BigO(Parent(f).1^0);
  new := ChangePrecision(new, Infinity());
  return new;
end function;

// should this really be an intrinsic?  Or should we just keep it as a function?
//function RiemannRochBasisHyperellipticCurveAnalytic(m, x_CC, y_CC, X_CC, Gamma)
intrinsic RiemannRochBasisHyperellipticAnalytic(m::RngIntElt, x_CC::RngSerLaurElt, y_CC::RngSerLaurElt, X_CC::Crv, Gamma::GrpPSL2) -> Any
  {input: m integer, X_CC complex curve, x_CC, y_CC power series in w. output: basis for L(m*(infty_1+infty_2)}
  g := Genus(Gamma);
  // require g ge 2: "Only for genus >= 2";
  assert g ge 2;
  assert g eq Genus(X_CC);
  //sigma := Gamma`TriangleSigma;
  //a := Order(sigma[1]);
  //e := a div #CycleDecomposition(sigma[1])[1];
  CCw<w> := Parent(x_CC);
  CC<I> := BaseRing(CCw);
  basis := [];
  if m le g then
    basis := basis cat [x_CC^i : i in [0..m]];
  else
    basis := basis cat [x_CC^i : i in [0..g]];
    for i := g+1 to m do
      Append(~basis, x_CC^i);
      Append(~basis, x_CC^(i-(g+1))*y_CC);
    end for;
  end if;
  // seems like even vs. odd no longer matters...
  // distinguish even and odd hyperelliptic cases
  // TODO: Now with stabilizers: check if works
  /*
  if Valuation(y_CC) eq -(2*g+1)*e then // this is the odd case TODO test this case
    if m le g then
      basis := basis cat [x_CC^i : i in [0..m]];
    else
      basis := basis cat [x_CC^i : i in [0..g]];
      for i := g+1 to m do
        Append(~basis, x^(i-(g+1))*y);
        Append(~basis, x^i);
      end for;
    end if;
  elif Valuation(y_CC) eq -(g+1)*e then // this is the even case
    if m le g then
      basis := basis cat [x_CC^i : i in [0..m]];
    else
      basis := basis cat [x_CC^i : i in [0..g]];
      for i := g+1 to m do
        Append(~basis, x^i);
        Append(~basis, x^(i-(g+1))*y);
      end for;
  else
    error "what the heck? not even or odd?!?!";
  end if;
  */
  return basis;
end intrinsic;
//end function;

intrinsic RiemannRochBasisHyperellipticExact(m::RngIntElt, X::CrvHyp) -> Any
  {}
  // input: m integer, X curve with algebraic coefficients
  // output: basis for L(m*(infty_1+infty_2))
  g := Genus(X);
  //K := BaseRing(X);
  KX<x,y> := FunctionField(X);
  basis := [];
  if m le g then
    basis := basis cat [x^i : i in [0..m]];
  else
    basis := basis cat [x^i : i in [0..g]];
    for i := g+1 to m do
      Append(~basis, x^(i-(g+1))*y);
      Append(~basis, x^i);
    end for;
  end if;
  return basis;
end intrinsic;

intrinsic TriangleHyperellipticNormalizeSeries(Gamma::GrpPSL2) -> Any
  {Normalize the series for x_CC and y_CC by completing the square and canceling unnecessary terms in Laurent expansion.}
  g := Genus(Gamma);
  x_CC, y_CC := Explode(Gamma`TriangleNewtonCoordinateSeries);
  curve_coeffs := Gamma`TriangleNumericalCurveCoefficients;
  CC<I> := Parent(curve_coeffs[1]);
  // make polys and complete the square
  CCt<t> := PolynomialRing(CC);
  u := CCt!0;
  v := CCt!0;
  for i in [0..2*g+2] do
    v -:= curve_coeffs[i+1]*t^i;
  end for;
  for i in [0..g+1] do
    u := curve_coeffs[2*g+4+i]*t^i; // should negative be positive?
  end for;
  v := v+u^2/4; // completing the square
  printf "After completing the square, v = %o\n", v;
  v_deg := Degree(v);
  c_tr := Coefficient(v, v_deg-1);
  v := Evaluate(v,t - c_tr/(v_deg-1));
  printf "After translating x, v = %o\n", v;
  X_CC := HyperellipticCurve(v);
  // remake curve coefficients
  curve_coeffs := Coefficients(-v);
  if #curve_coeffs eq 2*g+2 then
    vprintf Shimura : "Now odd model with y^2 = %o\n", v;
    Append(~curve_coeffs, 0);
  elif #curve_coeffs eq 2*g+3 then
    vprintf Shimura : "Now even model with y^2 = %o\n", v;
  else
    error "wtf, wrong degree";
  end if;
  assert #curve_coeffs eq 2*g+3;
  curve_coeffs cat:= [ 0 : i in [0..g+1]] cat [1];
  Gamma`TriangleNumericalCurveCoefficients := curve_coeffs;
  y_CC := y_CC + Evaluate(u,x_CC)/2; // change y_CC to completed square version
  x_CC := x_CC - c_tr/(v_deg-1); // change x_CC to translated version
  printf "After completing square, y_CC = %o\n", y_CC + O(Parent(y_CC).1^6);
  printf "After translating, x_CC = %o\n", x_CC + O(Parent(x_CC).1^6);
  /*
  // cancel terms in Laurent series
  y_val := Valuation(y_CC);
  y_CC := y_CC - PolarPart(y_CC);
  y_CC := y_CC - LeadingTerm(y_CC);
  y_CC := y_CC + Parent(y_CC).1^y_val;
  printf "After canceling terms in Laurent series, y_CC = %o\n", y_CC + O(Parent(y_CC).1^6);
  */
  Gamma`TriangleNewtonCoordinateSeries := [x_CC, y_CC];
  series := [x_CC^i : i in [0..2*g+2]] cat [x_CC^i*y_CC : i in [0..g+1]] cat [y_CC^2];
  Remove(~series, #series);
  Append(~series, y_CC^2);
  curve_vals := [Valuation(h) : h in series]; // reassign last val since y_CC has changed
  Gamma`TriangleCurveValuations := curve_vals;
  // TODO: probably also need to translate x to kill next to leading term in hyperelliptic poly
  return "Series for x and y normalized";
end intrinsic;

intrinsic TriangleHyperellipticTest(Sk::SeqEnum, Gamma::GrpPSL2) -> Any
  {
  Test if a given Belyi map is hyperelliptic.
  Input: Sk, an echelonized basis for the space of wt k modular forms, given as the output of TrianglePoweSeriesBasis;
         Gamma, triangle subgroup
  Output: hyp_bool, a boolean, true if exactly one hyperelliptic relation is found;
          curve_coeffs, a SeqEnum containing the coefficients of the hyperelliptic relation found;
          curve_vals, a SeqEnum containing the valuations of the functions in the basis for the curve coefficients.
  }
  prec := Precision(BaseRing(Parent(Sk[1][1])));
  eps := 10^(-prec/2);
  g := Genus(Gamma);
  //require g ge 2: "Only for genus >= 2";
  assert g ge 2;
  Delta := ContainingTriangleGroup(Gamma);
  phi, kappa := TrianglePhi(Delta);

  vprint Shimura: "Creating series for coordinate functions x and y...";
  CCw<w> := Parent(Sk[1][1]);
  CC<I> := BaseRing(CCw);
  f1 := Sk[1][1];
  f2 := Sk[2][1];
  assert g eq #Sk;
  f_g := Sk[g][1];
  //make f1, f2, f_g series in w/kappa, just like phi
  f1 := Evaluate(f1, kappa*w);
  // should make f1 monic?
  f2 := Evaluate(f2, kappa*w);
  f2 := RemoveLeadingZeros(f2, eps);
  f2 := f2/LeadingCoefficient(f2);
  f_g := Evaluate(f_g, kappa*w);
  f_g := RemoveLeadingZeros(f_g,eps);
  f_g := f_g/LeadingCoefficient(f_g);
  x_CC := f1/f2;
  y_CC := Derivative(x_CC)/f_g;
  y_CC := y_CC/LeadingCoefficient(y_CC);
  vprint Shimura: "Computing curve...";
  vprintf Shimura: "\tForming matrix of series coefficients...\n";
  // change to take odd and even into consideration: if odd only go up to x^2g+1
  series := [x_CC^i : i in [0..2*g+2]] cat [x_CC^i*y_CC : i in [0..g+1]] cat [y_CC^2];
  curve_vals := [Valuation(h) : h in series];
  Gamma`TriangleCurveValuations := curve_vals;
  // 2*g + 3 + g + 2 + 1 = 3*g + 6 things in series
  //M := Matrix([[Coefficient(f, i) : i in [2*Valuation(y_CC)..(3*g+5 + 2*Valuation(y_CC))] ] : f in series]);
  // TODO: How many columns?  The + 10 is a hack.
  // M := Matrix([[Coefficient(f, i) : i in [2*Valuation(y_CC)..(3*g+5 + 2*Valuation(y_CC)+10)] ] : f in series]);
  M := Matrix([[Coefficient(f, i) : i in [2*Valuation(y_CC)..(3*g+5 + 2*Valuation(y_CC)+15)] ] : f in series]); // just for one example
  //do I need to worry about e, the order of the stabilizer?
  assert Ncols(M) ge Nrows(M);
  vprintf Shimura: "\tComputing numerical kernel...\n";
  kerX := NumericalKernel(M : Epsilon := eps);
  kerX := KSpaceWithBasis(kerX);
  vprintf Shimura : "Basis of kernel = %o\n", Basis(kerX);
  if Dimension(kerX) eq 0 then
    print "No hyperelliptic relation found!";
    hyp_bool := false;
    Gamma`TriangleIsHyperelliptic := hyp_bool;
    return hyp_bool, _, _;
  elif Dimension(kerX) eq 1 then
    print "Exactly one hyperelliptic relation found!";
    hyp_bool := true;
    kerX := Basis(kerX);
    kerX := Eltseq(kerX[1]);
    curve_coeffs := [kerX[i]/kerX[#kerX] : i in [1..#kerX]];
    curve_coeffs := ZeroifyCoeffs(curve_coeffs,eps);
    Gamma`TriangleNumericalCurveCoefficients := curve_coeffs;
    // normalize series for x and y
    Gamma`TriangleNewtonCoordinateSeries := [x_CC, y_CC];
    TriangleHyperellipticNormalizeSeries(Gamma);
    // assign coordinate series to Gamma
    Gamma`TriangleIsHyperelliptic := hyp_bool;
    return hyp_bool, Gamma`TriangleNumericalCurveCoefficients, Gamma`TriangleNumericalCurveCoefficientsCurveValuations;
  else
    // TODO maybe this should print kernel and dimension?
    error "Multiple hyperelliptic relations found! Try a higher precision.";
  end if;
end intrinsic;

// function TriangleHyperellipticNumericalCoefficients(Sk, Gamma)
intrinsic TriangleHyperellipticNumericalCoefficients(Sk::SeqEnum, Gamma::GrpPSL2 : curve_coeffs := [], curve_vals := []) -> Any
  {Input: Sk, an echelonized basis for the space of wt k modular forms, given as the output of TrianglePoweSeriesBasis; Gamma, triangle subgroup
  Output: Coefficients of the curve, numerator, and denominator of the Belyi map, sequences of valuations of elements of the basis used in computing
  curve, numerator, and denominator. All coefficients are assigned to Gamma.}

  prec := Precision(BaseRing(Parent(Sk[1][1])));
  eps := 10^(-prec/2);
  g := Genus(Gamma);
  d := Gamma`TriangleD;
  sigma := Gamma`TriangleSigma;
  //require g ge 2: "Only for genus >= 2";
  assert g ge 2;
  CCw<w> := Parent(Sk[1][1]);
  CC<I> := BaseRing(CCw);
  Delta := ContainingTriangleGroup(Gamma);
  phi, kappa := TrianglePhi(Delta);
  vprint Shimura: "Creating series for coordinate functions x and y...";
  f1 := Sk[1][1];
  f2 := Sk[2][1];
  assert g eq #Sk;
  f_g := Sk[g][1];
  //make f1, f2, f_g series in w/kappa, just like phi
  f1 := Evaluate(f1, kappa*w);
  f2 := Evaluate(f2, kappa*w);
  f2 := RemoveLeadingZeros(f2, eps);
  f2 := f2/LeadingCoefficient(f2);
  f_g := Evaluate(f_g, kappa*w);
  f_g := RemoveLeadingZeros(f_g,eps);
  f_g := f_g/LeadingCoefficient(f_g);
  x_CC := f1/f2;
  y_CC := Derivative(x_CC)/f_g;
  y_CC := y_CC/LeadingCoefficient(y_CC);
  printf "x_CC = %o\n", x_CC + O(Parent(x_CC).1^3);
  printf "y_CC = %o\n", y_CC + O(Parent(y_CC).1^3);
  //printf "x_CC = %o\n", x_CC; // printing
  //printf "y_CC = %o\n\n", y_CC; // printing
  if #curve_coeffs eq 0 then // no pre-assigned curve_coeffs
    vprint Shimura: "Computing curve...";
    vprintf Shimura: "\tForming matrix of series coefficients...\n";
    series := [x_CC^i : i in [0..2*g+2]] cat [x_CC^i*y_CC : i in [0..g+1]] cat [y_CC^2];
    curve_vals := [Valuation(h) : h in series];
    Gamma`TriangleCurveValuations := curve_vals;
    // 2*g + 3 + g + 2 + 1 = 3*g + 6 things in series
    //M := Matrix([[Coefficient(f, i) : i in [2*Valuation(y_CC)..(3*g+5 + 2*Valuation(y_CC))] ] : f in series]);
    // TODO: How many columns?  The + 10 is a hack.
    M := Matrix([[Coefficient(f, i) : i in [2*Valuation(y_CC)..(3*g+5 + 2*Valuation(y_CC)+10)] ] : f in series]);
    //do I need to worry about e, the order of the stabilizer?
    assert Ncols(M) ge Nrows(M);
    //print Nrows(M), Ncols(M); // printing
    vprintf Shimura: "\tComputing numerical kernel...\n";
    kerX := NumericalKernel(M : Epsilon := eps);
    kerX := KSpaceWithBasis(kerX);
    assert Dimension(kerX) eq 1;
    curve_coeffs := Basis(kerX);
    //printf "curve_coeffs = %o\n", curve_coeffs; //printing
    curve_coeffs := Eltseq(curve_coeffs[1]);
    curve_coeffs := [curve_coeffs[i]/curve_coeffs[#curve_coeffs] : i in [1..#curve_coeffs]];
    curve_coeffs := ZeroifyCoeffs(curve_coeffs,eps);
    Gamma`TriangleNewtonCoordinateSeries := [x_CC, y_CC];
    Gamma`TriangleNumericalCurveCoefficients := curve_coeffs;
    printf "curve_coeffs = %o\n", curve_coeffs;
    vprintf Shimura : "...curve found!\n";
    vprint Shimura : "...normalizing series for x and y";
    TriangleHyperellipticNormalizeSeries(Gamma);
  end if;

  CCt<t> := PolynomialRing(CC);
  u := CCt!0;
  v := CCt!0;
  for i in [0..2*g+2] do
    v -:= curve_coeffs[i+1]*t^i;
  end for;
  for i in [0..g+1] do
    u := curve_coeffs[2*g+4+i]*t^i; // +/-?
  end for;
  assert u eq 0;
  X_CC := HyperellipticCurve(v,u);
  //printf "u = %o\n", u; // printing
  //printf "v = %o\n", v; // printing
  //printf "X_CC is %o\n\n", X_CC; // printing

  vprint Shimura: "Computing Belyi map...";
  t := Ceiling((d+g)/2);
  vprintf Shimura: "\tComputing Riemann-Roch spaces...\n";
  numbasis := RiemannRochBasisHyperellipticAnalytic(t, x_CC, y_CC, X_CC, Gamma);
  num_vals := [Valuation(h) : h in numbasis];
  denombasis := numbasis;
  denom_vals := num_vals;
  printf "\tt = %o\n", t; // printing
  a := DefiningABC(Gamma)[1];
  e := a div #CycleDecomposition(sigma[1])[1];
  printf "\te = %o\n", e; // printing

  vprintf Shimura: "\tForming matrix of series coefficients...\n";
  series := [f*phi : f in denombasis] cat numbasis;
  minval := Min([Valuation(f) : f in series]);
  printf "\tminval = %o\n", minval; // printing
  M := Matrix([[Coefficient(f,e*n) : n in [minval..minval + #numbasis + #denombasis + 20]] : f in series]);
  printf "\tnrows = %o, ncols = %o\n", Nrows(M), Ncols(M); // printing
  vprint Shimura : "\tComputing numerical kernel...";
  kerphi := NumericalKernel(M : Epsilon := eps);
  kerphi := KSpaceWithBasis(kerphi);
  printf "\tdimension of kernel = %o\n", Dimension(kerphi);
  if Dimension(kerphi) eq 0 then
    error "Matrix had trivial kernel. Try a higher precision.";
  end if;
  assert Dimension(kerphi) ge 1;
  Gamma`TriangleRiemannRochParameters := [t];
  kerphi := Basis(kerphi);
  kerphi := [Eltseq(el) : el in kerphi];
  printf "kerphi has length %o\n", #kerphi; //printing
  kerphi := [ZeroifyCoeffs(el,eps) : el in kerphi];

  kerphi_perm := [];
  for v in kerphi do
    v_new := [];
    for i := 1 to #numbasis do // reorder to put basis vectors with same order pole consecutively
      Append(~v_new, v[i]);
      Append(~v_new, v[#numbasis+i]); //denombasis?
    end for;
    Append(~kerphi_perm, v_new);
  end for;
  K_perm := Matrix(kerphi_perm);
  col_num := Ncols(K_perm);
  while (Abs(Norm(Transpose(K_perm)[col_num])) lt eps) and (col_num ge 1) do
    col_num := col_num - 1;
  end while;
  K := Submatrix(K_perm,1,1,Nrows(K_perm),col_num);
  Q, L := QLDecomposition(K);
  printf "lower triangle matrix L = %o\n", L;
  v := Rows(L)[1]; // top row of L is distinguished
  v := Eltseq(v);
  v := v cat [0 : i in [1..(Ncols(K_perm)-col_num)]];
  v := ZeroifyCoeffs(v,eps);
  lc_ind := #v;
  // find last nonzero coeff and divide through
  while (lc_ind ge 1) and (v[lc_ind] eq 0) do
    lc_ind := lc_ind - 1;
  end while;
  //kerphi := [kerphi[i]/kerphi[#kerphi] : i in [1..#kerphi]];
  v := [v[i]/v[lc_ind] : i in [1..#v]];

  // put basis functions back in original order
  v_new := [];
  for i := 1 to Ceiling(#v div 2) do
    Append(~v_new, v[2*i-1]);
  end for;
  for i := 1 to Floor(#v div 2) do
    Append(~v_new, v[2*i]);
  end for;
  // split into num and denom
  denom_coeffs := v_new[1..#denombasis];
  num_coeffs := v_new[#denombasis+1..#denombasis+#numbasis];
  //printf "#denom_coeffs before removing leading zeroes = %o\n", #denom_coeffs;
  //printf "#num_coeffs before removing leading zeroes = %o\n", #num_coeffs;
  while denom_coeffs[#denom_coeffs] eq 0 do
    Remove(~denom_coeffs,#denom_coeffs);
  end while;
  while num_coeffs[#num_coeffs] eq 0 do
    Remove(~num_coeffs,#num_coeffs);
  end while;
  denom_zero_inds := [];
  for i := 1 to #denom_coeffs do
    if denom_coeffs[i] eq 0 then
      Append(~denom_zero_inds, i);
    end if;
  end for;
  num_zero_inds := [];
  for i := 1 to #num_coeffs do
    if num_coeffs[i] eq 0 then
      Append(~num_zero_inds, i);
    end if;
  end for;
  Gamma`TriangleNewtonDenominatorZeroIndices := denom_zero_inds;
  Gamma`TriangleNewtonNumeratorZeroIndices := num_zero_inds;
  //printf "#denom_coeffs after removing leading zeroes = %o\n", #denom_coeffs;
  //printf "#num_coeffs after removing leading zeroes = %o\n", #num_coeffs;
  printf "unscaled curve coeffs = %o\n", curve_coeffs;
  printf "unscaled numerator coeffs = %o\n", num_coeffs;
  printf "unscaled denominator coeffs = %o\n", denom_coeffs;
  Gamma`TriangleUnscaledNumericalCurveCoefficients := curve_coeffs;
  Gamma`TriangleUnscaledNumericalBelyiMapNumeratorCoefficients := num_coeffs;
  Gamma`TriangleUnscaledNumericalBelyiMapDenominatorCoefficients := denom_coeffs;

  num_vals := num_vals[1..#num_coeffs];
  denom_vals := denom_vals[1..#denom_coeffs];
  //assign the vals
  Gamma`TriangleBelyiMapNumeratorValuations := num_vals;
  Gamma`TriangleBelyiMapDenominatorValuations := denom_vals;

  // rescale
  lambda, curve_coeffs, lc, num_coeffs, denom_coeffs := TriangleRescaleCoefficients(Gamma, [curve_coeffs, num_coeffs, denom_coeffs], [curve_vals, num_vals, denom_vals]);
  // write attributes for NewtonHyperelliptic

  // FIXME: I don't think this is right: think I need to do x_CC := Evaluate(x_CC, lambda*Parent(x_CC).1)
  // and then renormalize to make monic...
  // Or maybe need to take kappa out again, as in lines 455-457 of genusone.m
  // TODO: need to rescale series...check if right
  /*
  x_CC := Evaluate(x_CC,(1/kappa)*Parent(x_CC).1);
  y_CC := Evaluate(y_CC,(1/kappa)*Parent(x_CC).1);
  x_CC := Evaluate(x_CC,(1/lambda)*Parent(x_CC).1);
  y_CC := Evaluate(y_CC,(1/lambda)*Parent(x_CC).1);
  x_CC := x_CC/LeadingCoefficient(x_CC);
  y_CC := y_CC/LeadingCoefficient(y_CC);
  */
  x_CC := x_CC*lambda^(-Valuation(x_CC));
  y_CC := y_CC*lambda^(-Valuation(y_CC));
  printf "x_CC = %o\n", x_CC + O(Parent(x_CC).1^6);
  printf "y_CC = %o\n", y_CC + O(Parent(x_CC).1^6);
  Gamma`TriangleNewtonCoordinateSeries := [x_CC, y_CC];
  // write numerical attributes
  Gamma`TriangleNumericalCurveCoefficients := curve_coeffs;
  Gamma`TriangleNumericalBelyiMapLeadingCoefficient := lc;
  Gamma`TriangleNumericalBelyiMapNumeratorCoefficients := num_coeffs;
  Gamma`TriangleNumericalBelyiMapDenominatorCoefficients := denom_coeffs;
  Gamma`TriangleRescalingFactor := lambda;
  return Gamma, curve_coeffs, lc, num_coeffs, denom_coeffs, curve_vals, num_vals, denom_vals;
end intrinsic;
// end function;

/* local sanity check */

intrinsic VarText(var::MonStgElt, lower::RngIntElt, upper::RngIntElt) -> MonStgElt
  {returns text "varlower, varlower+1, ..., varupper-1, varupper".}
  assert upper ge lower;
  var_text := "";
  if upper eq lower then
    var_text *:= Sprintf("%o%o", var, lower);
  else
    for i := lower to upper-1 by 1 do
      var_text *:= Sprintf("%o%o, ", var, i);
    end for;
    var_text *:= Sprintf("%o%o", var, upper);
  end if;
  return var_text;
end intrinsic;

intrinsic VarSeq(var::MonStgElt, lower::RngIntElt, upper::RngIntElt) -> SeqEnum[MonStgElt]
  {returns SeqEnum ["varlower", "varlower+1", ..., "varupper-1", "varupper"].}
  assert upper ge lower;
  var_seq := [];
  for i := lower to upper do
    Append(~var_seq, Sprintf("%o%o", var, i));
  end for;
  return var_seq;
end intrinsic;

intrinsic HomText(var::MonStgElt, lower::RngIntElt, upper::RngIntElt) -> MonStgElt
  {returns text "var.lower, var.lower+1, ..., var.upper"}
  assert upper ge lower;
  var_text := "";
  if upper eq lower then
    var_text *:= Sprintf("%o.%o", var, lower);
  else
    for i := lower to upper-1 by 1 do
      var_text *:= Sprintf("%o.%o, ", var, i);
    end for;
    var_text *:= Sprintf("%o.%o", var, upper);
  end if;
  return var_text;
end intrinsic;

intrinsic BelyiPolynomialReduction(poly::RngMPolElt, mp::Map, P::RngMPol) -> RngMPolElt
  {Given a poly and a mp on coefficients, return poly with new coeffs.}
  assert Codomain(mp) eq BaseRing(P);
  rank := Rank(P);
  assert rank eq Rank(Parent(poly));
  h := eval Sprintf("h := hom<Parent(poly)->P|[%o]>; return h;", VarText("P.", 1, rank));
  coeffs, mons := CoefficientsAndMonomials(poly);
  poly_pp := P!0;
  for i := 1 to #coeffs do
    poly_pp +:= P!(mp(coeffs[i])*h(mons[i]));
  end for;
  return poly_pp;
end intrinsic;

intrinsic BelyiReduceCurve(X::Crv, f::FldFunFracSchElt, p::RngIntElt) -> Crv, FldFunFracSchElt
  {Reduce X and f in QQ(X) mod p and return X mod p, f mod p.}
  // setup
    K := BaseField(X);
    if not IsProjective(X) then
      X := ProjectiveClosure(X);
    end if;
    ZK := Integers(K); // works for any K
    pp := Factorization(p*ZK)[1][1];
    FFq, mpZKtoFFq := ResidueClassField(pp);
    I := Ideal(X);
  // reduce I mod pp
    equations := Basis(I);
    equations_pp := []; // equations for Ipp
    grading := Grading(I);
    P := PolynomialRing(FFq, grading); // grading for CrvHyp
    for eqn in equations do
      eqn_pp := BelyiPolynomialReduction(eqn, mpZKtoFFq, P);
      Append(~equations_pp, eqn_pp);
    end for;
    Ipp := ideal<P|equations_pp>;
  // make new curve and coordinate rings
    PP := ProjectiveSpace(Generic(Ipp));
    Xpp := Curve(PP, Ipp);
    KXpp := FunctionField(Xpp);
    num := Numerator(f);
    den := Denominator(f);
    Aff := Parent(num);
    Affpp := Parent(Numerator(KXpp.1));
    h := eval Sprintf("hompp := hom<Aff->Affpp|[%o]>; return hompp;", VarText("Affpp.", 1, Rank(Affpp)));
  // make num in Affpp
    num_coeffs, num_mons := CoefficientsAndMonomials(num);
    numpp := Affpp!0;
    for i := 1 to #num_coeffs do
      numpp +:= Affpp!(mpZKtoFFq(num_coeffs[i])*h(num_mons[i]));
    end for;
  // make den in Affpp
    den_coeffs, den_mons := CoefficientsAndMonomials(den);
    denpp := Affpp!0;
    for i := 1 to #den_coeffs do
      denpp +:= Affpp!(mpZKtoFFq(den_coeffs[i])*h(den_mons[i]));
    end for;
  // coerce f into KXpp and return
    fpp := KXpp!(numpp)/KXpp!(denpp);
    return Xpp, fpp;
end intrinsic;

/*
intrinsic BelyiMapSanityCheck(sigma::SeqEnum[GrpPermElt], X::Crv, phi::FldFunFracSchElt) -> Any
  {}
end intrinsic;
*/

/* for a single prime p */

intrinsic BelyiLocalSanityCheck(sigma::SeqEnum[GrpPermElt], X::Crv, phi::FldFunFracSchElt, p::RngIntElt) -> BoolElt
  {BelyiMapSanityCheck...Localified...no lax!}
  C, mapp := BelyiReduceCurve(X, phi, p);
  if not BelyiMapSanityCheck(sigma, C, mapp) then
    vprintf Shimura : "Local Sanity Check Failed for p=%o:\n", p;
    vprintf Shimura : "sigma = \n%o.\n", sigma;
    supp, mult := Support(Divisor(mapp));
    vprintf Shimura : "supp(phi) = \n%o\n%o.\n", supp, mult;
    supp1, mult1 := Support(Divisor(mapp-1));
    vprintf Shimura : "supp(phi-1) = \n%o\n%o.\n", supp1, mult1;
    return false;
  end if;
  // if we make it here then we passed!
  return true;
end intrinsic;

intrinsic BelyiLocalSanityCheck(Gamma::GrpPSL2, p::RngIntElt) -> BoolElt
  {BelyiMapSanityCheck...Localified...no lax!}
  if assigned Gamma`TriangleBelyiCurve then
    X := Gamma`TriangleBelyiCurve;
    phi := Gamma`TriangleBelyiMap;
    sigma := Gamma`TriangleSigma;
    return BelyiLocalSanityCheck(sigma, X, phi, p);
  else
    return false;
  end if;
end intrinsic;

/* for several primes */

intrinsic BelyiLocalSanityCheck(sigma::SeqEnum[GrpPermElt], X::Crv, phi::FldFunFracSchElt) -> BoolElt
  {}
  try
    bool := BelyiLocalSanityCheck(sigma, X, phi, 9721);
    if bool then
      return true;
    else
      return false;
    end if;
  catch e
    print "p=9721 bad";
    try
      bool := BelyiLocalSanityCheck(sigma, X, phi, 101);
      if bool then
        return true;
      else
        return false;
      end if;
    catch e
      print "p=101 bad";
      try
        bool := BelyiLocalSanityCheck(sigma, X, phi, 17);
        if bool then
          return true;
        else
          return false;
        end if;
      catch e
        print "p=17 bad";
        error "curve can't be reduced at all 3 primes!";
      end try;
    end try;
  end try;
end intrinsic;

intrinsic BelyiLocalSanityCheck(Gamma::GrpPSL2, p::RngIntElt) -> BoolElt
  {}
  if assigned Gamma`TriangleBelyiCurve then
    X := Gamma`TriangleBelyiCurve;
    phi := Gamma`TriangleBelyiMap;
    sigma := Gamma`TriangleSigma;
    return BelyiLocalSanityCheck(sigma, X, phi);
  else
    return false;
  end if;
end intrinsic;

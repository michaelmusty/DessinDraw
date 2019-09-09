intrinsic TriangleTheta(Sk::SeqEnum[SeqEnum[RngSerPowElt]], Gamma::GrpPSL2 :
                        Al := "ByValues", Recompute := false) -> Any
{}
  if not Recompute and assigned Gamma`TriangleTheta and
      Precision(Parent(Gamma`TriangleTheta)) ge Precision(BaseRing(Parent(Sk[1][1]))) then
    return Gamma`TriangleTheta;
  end if;

  prec := Precision(BaseRing(Parent(Sk[1][1])));
  eps := 10^(-prec/2);

  aut := #AutomorphismGroup(Gamma);
  if aut gt 1 then
    vprintf Shimura : "Gulp, automorphism group size > 1!\n";
  end if;

  if Al eq "Ratio" then
    // First take the one with the smallest stabilizer
    // TODO: Use finer information from automorphism group here (PointedReps?)
    sigma := DefiningPermutation(Gamma);
    _, gap1 := Max([#c : c in CycleDecomposition(sigma[1])]);
    Sk1 := [f[gap1] : f in Sk];

    // Find the element in the basis with the smallest "gap" between nonzero coefficients
    gaps := [];
    ss := [];
    for f in Sk1 do
      s := Degree(LeadingTerm(f));
      while Abs(Coefficient(f,s)) lt eps do
        s +:= 1;
      end while;
      Append(~ss, s);
      sj := s+1;
      while Abs(Coefficient(f,sj)) lt eps do
        sj +:= 1;
      end while;
      Append(~gaps, sj-s);
    end for;
    mingap, gapind := Min(gaps);
    vprintf Shimura : "Taking gap = %o with index %o\n", mingap, gapind;

    f := Sk1[gapind];
    s := ss[gapind];
    n := s+mingap;

    if mingap gt 1 and Abs(Coefficient(f,n+1)) gt eps then
      s := n;
      n := s+1;
      mingap := 1;
    end if;
    // assert mingap le 2;  // otherwise have to be clever about a root of unity?

    // check if "next" coefficient is nonvanishing, since simpler
    // TODO: Otherwise use weighted scaling. Also, can deal with larger
    // mingap in same way.
    Theta := (Factorial(n-s)*Coefficient(f,n)/Coefficient(f,s))^(1/mingap);
  else // Al eq "ByValues"
    // Normalize so that the sum of the values is 1.
    vals := TriangleRamificationValues(Sk, Gamma : NormalizeByTheta := false);
    Thetainvs := &cat[ [ &+[vals[s][i][j+1]/vals[s][i][j] : i in [1..#vals[1]] | Abs(vals[s][i][j]) gt eps] : j in [1..#vals[1][1]-1]] : s in [1..3]];
    Thetainvs := [Thetainv : Thetainv in Thetainvs | Abs(Thetainv) gt eps and Abs(Thetainv) lt 1/eps];
    if #Thetainvs eq 0 then
      vprintf Shimura : "Tried Al := 'ByValues' but failed, so now using 'Ratio'\n";
      Theta := TriangleTheta(Sk, Gamma : Al := "Ratio", Recompute := true);
    else
      _, thind := Min([Abs(Abs(Thetainv)-1) : Thetainv in Thetainvs]);  // choose close to 1 in absval
      Theta := Thetainvs[thind]^-1;
    end if;
  end if;

  Gamma`TriangleTheta := Theta;

  vprintf Shimura : "Theta = %o\n", Theta;

  return Theta;
end intrinsic;

/*
intrinsic TriangleMakeNumberField(uCC::FldComElt, deg::RngIntElt : K := "") -> Any
  {Given a complex number uCC thought to be algebraic, return the optimized representation of the number field Kop generated by uCC, along with embedding data consisting of a place of Kop and a boolean for complex conjugation.}

  
  if m eq 1 then
    vprint Shimura : "  ...m eq 1, so taking QQ";
    K := RationalsAsNumberField();
    v := RealPlaces(K)[1];
    return true, K, v, false, Parent(uCC)!1;
  end if;

  vprintf Shimura : "  ...Trying to MakeK with m = %o\n", m;

  CC := Parent(uCC);
  eps := 10^(-Precision(CC)/2);

  u_pol := PowerRelation(uCC, m : Al := "LLL");
  if Degree(u_pol) ne m or not IsIrreducible(u_pol) then
    if Degree(u_pol) gt 1 then
      vprint Shimura : "  ...coefficient not degree m, looked like degree", [<Degree(c[1]),c[2]> : c in Factorization(u_pol)];
    end if;
    return false, _, _, _, _;
  end if;
  lc := LeadingCoefficient(u_pol);
  K0 := NumberField(u_pol);
  u0 := K0.1;
  K := NumberField(MinimalPolynomial(lc*u0));
  u := K.1/lc;
  ps, foobar := TrialDivision(Discriminant(K));
  if #foobar gt 0 and not &and[IsSquare(foobarf) : foobarf in foobar] then
    return false, _, _, _, _;
  end if;

  if foobar eq [] then
    vprintf Shimura : "  ...coefficient found, with ps = %o and foobar = []\n", ps;
  else
    vprintf Shimura : "  ...coefficient found, with ps = %o and sqrt(foobar) = %o\n", ps, Round(Sqrt(CC!foobar[1])); 
  end if;
  vprintf Shimura : "  ...%o\n", K;
  vprint Shimura : "  ...Trying to optimize"; 
  
  Kop, iotaop := OptimizedRepresentation(K : Ramification := [p[1] : p in ps]);
  uop := iotaop(u);
  
  vprintf Shimura : "  ...successfully optimized\n  Kop = %o\n  now finding complex embedding\n", Kop;
  minv, vind := Min([Abs(uCC-Evaluate(uop,v : Precision := Precision(CC))) : v in InfinitePlaces(Kop)]);
  if minv gt eps then
    // Need complex conjugate
    minv, vind := Min([Abs(Conjugate(uCC)-Evaluate(uop,v : Precision := Precision(CC))) : v in InfinitePlaces(Kop)]);
    conj := true;
  else
    conj := false;
  end if;
  if minv gt eps then
    vprint Shimura : "  ...failed to match; got a polynomial but it apparently sux";
    return false, _, _, _, _;
  end if;
  vprint Shimura : "  ...quality of match", RealField(4)!minv;
  v := InfinitePlaces(Kop)[vind];
  // return true, Kop, v, conj, uop;
  return true, Kop, v, conj, uCC;
end intrinsic;
*/

intrinsic MakeK(uCC::Any, m::Any : UsePolredabs := true) -> Any, Any, Any, Any, Any
  {MakeK!  What more to say?}

  if m eq 1 then
    vprint Shimura : "  ...m eq 1, so taking QQ";
    K := RationalsAsNumberField();
    v := RealPlaces(K)[1];
    return true, K, v, false, Parent(uCC)!1;
  end if;

  vprintf Shimura : "  ...Trying to MakeK with m = %o\n", m;

  CC := Parent(uCC);
  eps := 10^(-Precision(CC)/2);

  u_pol := PowerRelation(uCC, m : Al := "LLL");
  if Degree(u_pol) ne m or not IsIrreducible(u_pol) then
    if Degree(u_pol) gt 1 then
      vprint Shimura : "  ...coefficient not degree m, looked like degree", [<Degree(c[1]),c[2]> : c in Factorization(u_pol)];
    end if;
    return false, _, _, _, _;
  end if;
  lc := LeadingCoefficient(u_pol);
  K0 := NumberField(u_pol);
  u0 := K0.1;
  K := NumberField(MinimalPolynomial(lc*u0));
  u := K.1/lc;
  ps, foobar := TrialDivision(Discriminant(K));
  if #foobar gt 0 and not &and[IsSquare(foobarf) : foobarf in foobar] then
    return false, _, _, _, _;
  end if;

  if foobar eq [] then
    vprintf Shimura : "  ...coefficient found, with ps = %o and foobar = []\n", ps;
  else
    vprintf Shimura : "  ...coefficient found, with ps = %o and sqrt(foobar) = %o\n", ps, Round(Sqrt(CC!foobar[1])); 
  end if;
  vprintf Shimura : "  ...%o\n", K;
  vprint Shimura : "  ...Trying to optimize"; 

  if UsePolredabs then
    Kop, iotaop, bl := Polredabs(K);
    vprint Shimura : "Polredabs done!";
  end if;
  if not UsePolredabs or (UsePolredabs and not bl) then
    vprint Shimura : "Polredabs failed, using optimized representation";
    Kop, iotaop := OptimizedRepresentation(K : Ramification := [p[1] : p in ps]);
  end if;
  uop := iotaop(u);
  
  vprintf Shimura : "  ...successfully optimized\n  Kop = %o\n  now finding complex embedding\n", Kop;

  /*  
  CCdefault := Evaluate(uop,InfinitePlaces(Kop)[1]);
  precval := Min(Precision(CCdefault),Precision(CC));
  
  minv, vind := Min([Abs(uCC-Evaluate(uop,v : Precision := precval)) : v in InfinitePlaces(Kop)]);
  if minv gt 10^(-precval/2) then
    // Need complex conjugate
    minv, vind := Min([Abs(Conjugate(uCC)-Evaluate(uop,v : Precision := precval)) : v in InfinitePlaces(Kop)]);
    conj := true;
  else
    conj := false;
  end if;
  if minv gt 10^(-precval/2) then
    vprint Shimura : "  ...failed to match; got a polynomial but it apparently sux";
    return false, _, _, _, _;
  end if;
  */
  minv, vind := Min([Abs(uCC-Evaluate(uop,v : Precision := Precision(CC))) : v in InfinitePlaces(Kop)]);
  if minv gt eps then
    // Need complex conjugate
    minv, vind := Min([Abs(Conjugate(uCC)-Evaluate(uop,v : Precision := Precision(CC))) : v in InfinitePlaces(Kop)]);
    conj := true;
  else
    conj := false;
  end if;
  if minv gt eps then
    vprint Shimura : "  ...failed to match; got a polynomial but it apparently sux";
    return false, _, _, _, _;
  end if;
  vprint Shimura : "  ...quality of match", RealField(4)!minv;
  v := InfinitePlaces(Kop)[vind];
  // return true, Kop, v, conj, uop;
  return true, Kop, v, conj, uCC;
end intrinsic;

intrinsic TriangleK(Gamma::GrpPSL2) -> Any, Any, Any
{}
  return Gamma`TriangleK, Gamma`TriangleKv, Gamma`TriangleKIsConjugated;
end intrinsic;

intrinsic TriangleK(Sk::SeqEnum[SeqEnum[RngSerPowElt]], Gamma::GrpPSL2 : m := 0) -> Any, Any, Any
{}
  CCw := Parent(Sk[1][1]);
  CC := BaseRing(CCw);
  eps := 10^(-Precision(CC)/2);

  Theta := TriangleTheta(Sk, Gamma);

  if m eq 0 then
    sigma := DefiningPermutation(Gamma);
    m := #PassportRepresentatives(sigma : Pointed := true);
  end if;
  if m eq 1 then
    Gamma`TriangleK := RationalsAsNumberField();
    Gamma`TriangleKv := RealPlaces(Gamma`TriangleK)[1];
    Gamma`TriangleKIsConjugated := false;
    Gamma`TriangleKNumericalGenerator := CC!1;
    return Gamma`TriangleK, Gamma`TriangleKv, Gamma`TriangleKIsConjugated, Gamma`TriangleKNumericalGenerator;
  end if;

  Skc := [];
  for f in Sk do
    s := Degree(LeadingTerm(f[1]));
    while Abs(Coefficient(f[1],s)) lt eps do
      s +:= 1;
    end while;
    fc := [Coefficient(f[1],n)/Theta^(n-s) : n in [0..Precision(Parent(f[1]))-1]];
    Append(~Skc, fc);
  end for;

  // Compute number field
  fc := Skc[1];
  n := 0;
  repeat
    n +:= 1;
    cn := fc[n];
    if m eq 1 then
      assert Abs(Im(cn)) lt eps;  // if rational, should be real!
      cn := Re(cn);
    end if;

    bl, K, v, conj := MakeK(cn, m);
  until bl;

  Gamma`TriangleK := K;
  Gamma`TriangleKv := v;
  Gamma`TriangleKIsConjugated := conj;
  Gamma`TriangleKNumericalGenerator := fc[n];

  return K, v, conj, fc[n];
end intrinsic;

intrinsic RecognizeSeries(Sk::SeqEnum[SeqEnum[RngSerPowElt]], Gamma::GrpPSL2 : KeepTheta := true, m := 0) -> Any
{}
  CCw := Parent(Sk[1][1]);
  CC := BaseRing(CCw);
  eps := 10^(-Precision(CC)/2);
  NN := Precision(CCw);

  Theta := TriangleTheta(Sk, Gamma);
  K, v, conj, uCC := TriangleK(Sk, Gamma : m := m);
  if m eq 0 then
    m := Degree(K);
  end if;

  ZK := Integers(K);
  ZKbCC := [Evaluate(b,v : Precision := Precision(CC)) : b in Basis(ZK)];
  if conj then
    ZKbCC := [Conjugate(zkb) : zkb in ZKbCC];
  end if;
  if not KeepTheta then
    // Adjust Theta
    Skc := [];
    for f in Sk do
      s := Degree(LeadingTerm(f[1]));
      while Abs(Coefficient(f[1],s)) lt eps do
        s +:= 1;
      end while;
      fc := [Coefficient(f[1],n)/Theta^(n-s) : n in [0..Precision(Parent(f[1]))-1]];
      Append(~Skc, fc);
    end for;

    kappapass := 0;
    repeat
      fb := [];
      for n := 0 to #fc-1 do
        if m eq 1 then
          lin := LinearRelation(ZKbCC cat [-fc[n+1]] cat [CC.1] : Al := "LLL");
          // ignore imaginary part   TODO: Check it in range
        else
          lin := LinearRelation(ZKbCC cat [-fc[n+1]] : Al := "LLL");
        end if;
        if lin[m+1] eq 0 then
          break;
        end if;
        Append(~fb, (ZK!lin[1..m])/lin[m+1]);  // no factorial here!
      end for;

      // only take first two recognized coefficients at first
      cntrec := [];
      n := 1;
      while fb[n] eq 0 do
        n +:= 1;
      end while;
      s := n-1;
      n +:= 1;
      while #cntrec lt 2 do
        if fb[n] ne 0 then
          Append(~cntrec, <n-1, fb[n]>);
        end if;
        n +:= 1;
      end while;
      a := DefiningABC(Gamma)[1];
      t := cntrec[2][1]-cntrec[1][1];
      toa := Gcd(t,a);
      cntrec[1][1] := (cntrec[1][1]-s) div toa;
      cntrec[2][1] := (cntrec[2][1]-s) div toa;

      if m eq 1 then
        bb1 := &*([1] cat [ pp[1]^(Floor(pp[2]/cntrec[1][1])) : pp in Factorization(Numerator(cntrec[1][2]))])*
               &*([1] cat [ pp[1]^(Floor(-pp[2]/cntrec[1][1])) : pp in Factorization(Denominator(cntrec[1][2]))]);
        bb2 := &*([1] cat [ pp[1]^(Floor(pp[2]/cntrec[2][1])) : pp in Factorization(Numerator(cntrec[2][2]))])*
               &*([1] cat [ pp[1]^(Floor(-pp[2]/cntrec[2][1])) : pp in Factorization(Denominator(cntrec[2][2]))]);
        beta := Gcd(Numerator(bb1),Numerator(bb2))/Lcm(Denominator(bb1),Denominator(bb2));
      else
        bb1 := &*([1*ZK] cat [ pp[1]^(Floor(pp[2]/cntrec[1][1])) : pp in Factorization(cntrec[1][2]*ZK)]);
        bb2 := &*([1*ZK] cat [ pp[1]^(Floor(pp[2]/cntrec[2][1])) : pp in Factorization(cntrec[2][2]*ZK)]);
        bb := bb1+bb2;
        ClK, mCl := ClassGroup(ZK);
        bb0 := mCl((bb^-1)@@mCl);
        bb *:= bb0;
        bl, beta := IsPrincipal(bb);
      end if;
      betaCC := Evaluate(beta, v : Precision := Precision(CC));
      Theta *:= betaCC^(1/toa);

      vprintf Shimura : "With beta = %o, now Theta = %o\n", beta, Theta;

      if t/Gcd(a,t) gt 1 and kappapass eq 0 then
        at := Lcm(a,t);
        _, kappa := TrianglePhi(Gamma);
        alphat := 1/(Theta*kappa)^at;
        if m eq 1 then
          lin := LinearRelation(ZKbCC cat [-alphat] cat [CC.1] : Al := "LLL");
        else
          lin := LinearRelation(ZKbCC cat [-alphat] : Al := "LLL");
        end if;
        alphat := (ZK!lin[1..m])/lin[m+1];
        if not IsPower(alphat, at div t) then
          Theta /:= (CC!Evaluate(alphat, v : Precision := Precision(CC)))^(1/(at div a));
        end if;

        vprintf Shimura : "Kappa adjustment, now Theta = %o\n", Theta;
        kappapass := 1;

        f := Sk[1];
        s := Degree(LeadingTerm(f[1]));
        while Abs(Coefficient(f[1],s)) lt eps do
          s +:= 1;
        end while;
        fc := [Factorial(n)*Coefficient(f[1],n)/Theta^(n-s) : n in [0..Precision(Parent(f[1]))-1]];
      else
        kappapass := 2;
      end if;
    until kappapass eq 2;

    Gamma`TriangleTheta := Theta;
  end if;

  Skc := [];
  for f in Sk do
    s := Degree(LeadingTerm(f[1]));
    while Abs(Coefficient(f[1],s)) lt eps do
      s +:= 1;
    end while;
    fc := [Coefficient(f[1],n)/Theta^(n-s) : n in [0..Precision(Parent(f[1]))-1]];
    Append(~Skc, fc);
  end for;

  Skb := RecognizeOverK(Skc, K, v, conj : escapeOK := true);

  KTw<Tw> := PowerSeriesRing(K, NN+1);
  Skb := [ KTw!fb + BigO(Tw^(#fb)) : fb in Skb];

  return Skb;
end intrinsic;

intrinsic RecognizeOverK(Skc::SeqEnum[SeqEnum[FldComElt]], K::FldAlg, v::PlcNumElt, conj::BoolElt : escapeOK := false) -> SeqEnum
{Recognizes a sequence of sequences as elements of K with respect to the embedding v, conjugated if conj.}

  CC := Parent(Skc[1][1]);
  eps := 10^(-Precision(CC)/3);

  ZK := Integers(K);
  ZKbCC := [CC!Evaluate(b,v : Precision := Precision(CC)) : b in Basis(ZK)];
  if conj then
    ZKbCC := [Conjugate(zkb) : zkb in ZKbCC];
  end if;

  m := Degree(K);

  Skb := [];
  for fc in Skc do
    fb := [K | ];
    for n := 0 to #fc-1 do
      if m eq 1 then
        lin := LinearRelation(ZKbCC cat [-Re(fc[n+1])] : Al := "LLL");
      else
        lin := LinearRelation(ZKbCC cat [-fc[n+1]] : Al := "LLL");
      end if;
      if lin[m+1] eq 0 then
        error "Failed to find linear relation in RecognizeOverK";
      end if;
      vprint Shimura : "Denominator", lin[m+1], "should be not too large!";
      Append(~fb, ((ZK!lin[1..m])/lin[m+1]));
      fbv := CC!Evaluate(fb[#fb],v : Precision := Precision(CC));
      if conj then fbv := Conjugate(fbv); end if;
      err := Abs(fbv - fc[n+1]);
      vprintf Shimura : "Coefficient %o of %o in array %o recognized, error = %o \n", n+1, #fc, Index(Skc, fc),
        RealField(4)!err;
      if Abs(err) gt eps then
        if escapeOK then fb := fb[1..#fb-1]; break n;
        else error "Insufficient precision in RecognizeOverK, and no way to abort!  Increase precision"; end if;
      end if;
    end for;
    Append(~Skb, fb);
  end for;

  return Skb;
end intrinsic;

//why FldAlg instead of FldNum?
intrinsic RecognizeOverK(z::FldComElt, K::FldAlg, v::PlcNumElt, conj::BoolElt) -> Any
{Recognizes a complex number as an element of K with respect to the embedding v, conjugated if conj.}
  CC := Parent(z);
  ZK := Integers(K);
  ZKbCC := [CC!Evaluate(b,v : Precision := Precision(CC)) : b in Basis(ZK)];
  if conj then
    ZKbCC := [Conjugate(zkb) : zkb in ZKbCC];
  end if;
  m := Degree(K);

  if m eq 1 then
    lin := LinearRelation(ZKbCC cat [-Re(z)] : Al := "LLL");
  else
    lin := LinearRelation(ZKbCC cat [-z] : Al := "LLL");
  end if;
  z_K := (ZK!lin[1..m])/lin[m+1];
  z_K := K!z_K;

  return z_K;
end intrinsic;

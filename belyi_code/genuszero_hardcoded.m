import "genuszero.m" : TriangleGenusZeroMakeCurveAndMap;

// We hard-coded a few nice families of genus zero Belyi maps; these are
// often corner cases for our algorithms, so it saves us some trouble
IsGenusZeroHardCoded := function(sigma);
  abc := [Order(s) : s in sigma];
  d := Degree(Universe(sigma));
  abc_sorted := Sort(abc);
  as, bs, cs := Explode(abc_sorted);
  
  cycs := [Sort(CycleStructure(sigmai)) : sigmai in sigma];
  
  // Anderson--Bouw--Ejder--Girgin--Karemaker--Manes, Dynamical Belyi maps, 2017, Proposition 3.4.
  // cycle structure (d-k,2k+1,d-k)
  if &and[cyc[1][1] eq 1 : cyc in cycs] and &and[#cyc eq 2 : cyc in cycs] 
     and &and[cyc[2][2] eq 1 : cyc in cycs] 
     and ((as eq bs and 2*(d-as)+1 eq cs) or (bs eq cs and 2*(d-bs)+1 eq as)) then
    if as eq bs then k := d-as; else k := d-bs; end if;
    _<x> := FunctionField(Rationals());
    acoeffs := [Factorial(k)*Binomial(d,i)*Binomial(d-k-i-1,k-i) : i in [0..k]];
    phi := x^(d-k)*&+[ (-1)^i*acoeffs[i+1]*x^(k-i) : i in [0..k]]/&+[ (-1)^(k-i)*acoeffs[k-i+1]*x^(k-i) : i in [0..k]];
    assert phi eq Evaluate(1/phi,1/x);
    if abc[1] eq abc[2] then
      phi := 1/(1-phi);
    elif abc[2] eq abc[3] then
      phi := 1-phi;
    end if;
    return true, phi;

  // Cameron--Kemp--Maslak--Melamed--Moy--Pham--Wei, Shabat polynomials..., 2018
  elif cs eq d then // try Shabat polynomials
    lambdas := [CycleStructure(s) : s in sigma];
    found := [];
    
    // need partitions [<s,r-1>,<t,1>] and [<1, (r-1)*(s-1)+(t-1)>, <r,1>] (and [<d,1>])
    // may assume r > 1, else cyclic
    for ijk in [[1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1]] do
      if lambdas[ijk[3]] ne [<d,1>] then continue; end if;
      if #lambdas[ijk[2]] eq 2 then continue; end if;
      tryr1 := Sort(lambdas[ijk[2]]); 
      if tryr1[1][1] eq 1 and tryr1[2][2] eq 1 then
        r := tryr1[2][1];
      else
        continue;
      end if;
      if #lambdas[ijk[1]] gt 2 then continue; end if;  // allow t = 0
      if {lam[2] : lam in lambdas[ijk[1]]} ne {1,r-1} then continue; end if;
      tryrst := lambdas[ijk[1]];
      if #tryrst eq 1 then
        s := tryrst[1][1];
        t := 0;
      else
        if tryrst[1][2] eq 1 then
          t := tryrst[1][1];
          s := tryrst[2][1];
        else
          t := tryrst[2][1];
          s := tryrst[1][1];
        end if;
      end if;
      assert d eq s*(r-1)+t;
      if t eq 0 then 
        assert lambdas[ijk[1]] eq [<s,r-1>];
      else
        assert SequenceToSet(lambdas[ijk[1]]) eq {<s,r-1>,<t,1>};
      end if;
      assert SequenceToSet(lambdas[ijk[2]]) eq {<1,(r-1)*(s-1)+(t-1)>,<r,1>};
      
      _<x> := PolynomialRing(Rationals());
      poch := function(a,k);
        return &*([1] cat [a+i : i in [0..k-1]]);
      end function;
      phi := (1-x)^t*(&+[poch(t/s,k)*x^k/Factorial(k) : k in [0..r-1]])^s;
      found := ijk;
    end for;
    
    if found eq [] then
    // need partitions [<1,r+t-2>,<r,1>,<t,1>] and [<2, r+t-1>] (and [<d,1>])
    for ijk in [[1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1]] do
      if lambdas[ijk[3]] ne [<d,1>] then continue; end if;
      if #lambdas[ijk[2]] ne 1 or lambdas[ijk[2]][1][1] ne 2 then continue; end if;
      if #lambdas[ijk[1]] ne 3 then continue; end if;
      tryr1 := Sort(lambdas[ijk[1]]); 
      if tryr1[1][1] eq 1 and tryr1[2][2] eq 1 and tryr1[3][2] eq 1 then
        r := tryr1[2][1];
        t := tryr1[3][1];
      else
        continue;
      end if;
      assert d eq 2*(r+t-1);
      assert SequenceToSet(lambdas[ijk[1]]) eq {<1,r+t-2>,<r,1>,<t,1>};
      assert SequenceToSet(lambdas[ijk[2]]) eq {<2,r+t-1>};
      
      _<x> := PolynomialRing(Rationals());
      phi := 4*x^r*(1-x)^t*(&+[Binomial(t-1+j,t-1)*x^j : j in [0..r-1]])*
                (&+[Binomial(r-1+j,r-1)*Binomial(r+t-1,r+j)*(-x)^j : j in [0..t-1]]);
      found := ijk;
    end for;
    end if;

    if found eq [] then
    // need partitions [<1,4*r-3>,<r,2>] and [<3, 2*r-1>] (and [<d,1>])
    for ijk in [[1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1]] do
      if lambdas[ijk[3]] ne [<d,1>] then continue; end if;
      if #lambdas[ijk[2]] ne 1 then continue; end if;
      tryr2 := lambdas[ijk[2]][1];
      if tryr2[1] ne 3 or tryr2[2] mod 2 eq 0 then continue; end if;
      r := (tryr2[2]+1) div 2;
      if #lambdas[ijk[1]] ne 2 then continue; end if;
      tryr1 := Sort(lambdas[ijk[1]]); 
      if tryr1 ne [<1,4*r-3>,<r,2>] then continue; end if;
      assert d eq 3*(2*r-1);
      
      K<rho> := NumberField(Polynomial([1,1,1]));
      s3 := 2*rho+1;
      _<xK> := PolynomialRing(K);
      Sr := (1-xK)^r*&+[Binomial(r-1+j,r-1)*xK^j : j in [0..r-1]];
      phiK := -3*s3*Evaluate(Sr,s3*xK-rho^2)*(1-Evaluate(Sr,s3*xK-rho^2))*(Evaluate(Sr,s3*xK-rho^2)+rho);
      phi := Polynomial(ChangeUniverse(Eltseq(phiK),Rationals()));
      found := ijk;
    end for;
    end if;
    
    if found ne [] then
      if found eq [1,3,2] then
        phi := phi/(1-phi);
      elif found eq [2,1,3] then
        phi := 1-phi;
      elif found eq [2,3,1] then
        phi := 1/(1-phi);  // may need to switch with the one below
      elif found eq [3,1,2] then
        phi := 1-1/phi;
      elif found eq [3,2,1] then
        phi := 1/phi;
      end if;
      return true, phi;
    end if;
  end if;
  return false, _; 
end function;


hardcode_test := function(sigma);
  bl, phi := IsGenusZeroHardCoded(sigma);
  if not bl then return true; end if;
  Gamma := TriangleSubgroup(sigma);
  phinum := Numerator(phi);
  phiden := Denominator(phi);
  Gamma`TriangleExactBelyiMapLeadingCoefficient := Rationals()!1;
  Gamma`TriangleExactBelyiMapNumeratorCoefficients := Eltseq(phinum);
  Gamma`TriangleExactBelyiMapDenominatorCoefficients := Eltseq(phiden);
  _ := TriangleGenusZeroMakeCurveAndMap(Gamma);
  return BelyiMapSanityCheck(Gamma), phi;
end function;

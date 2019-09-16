declare verbose Refined, 1;

intrinsic RefinedPassportJV(sigma::SeqEnum[GrpPermElt]) -> Any
  {}
  error "not implemented!";
  // setup
  if #sigma eq 2 then
    Append(~sigma, (sigma[2]*sigma[1])^-1);
  end if;
  assert #sigma eq 3;
  d := Degree(Parent(sigma[1]));
  G := sub<Sym(d)|sigma>;
  t0 := Cputime();
  C := ClassMap(G);
  t1 := Cputime();
  vprintf Refined : "Class map computed: %o s\n", t1-t0;
  Csigma := [C(si) : si in sigma];
  vprintf Refined : "Classes of sigma: %o\n", Csigma;
  // passport lemma
  tau := sigma[1..2];
  t0 := Cputime();
  passlemma := PassportLemma(G, tau : generates := true);
  t1 := Cputime();
  vprintf Refined : "Passport lemma: %o s\n", t1-t0;
  vprintf Refined : "#passlemma = %o\n", #passlemma;
  t0 := Cputime();
  passlemma := [sigmap : sigmap in passlemma | [C(sip) : sip in sigmap] eq Csigma[1..2]];
  t1 := Cputime();
  vprintf Refined : "#passlemma matching %o = %o\n", Csigma[1..2], #passlemma;
  // make triples
  sigmaps := [];
  for sigmap in passlemma do
    sigmapoo := (sigmap[2]*sigmap[1])^-1;
    if C(sigmapoo) eq Csigma[3] then
      Append(~sigmaps, sigmap cat [sigmapoo]);
    end if;
  end for;
  vprintf Refined : "#triples matching %o = %o\n", Csigma, #sigmaps;
  // outer automorphisms
  t0 := Cputime();
  AutG := AutomorphismGroup(G);
  t1 := Cputime();
  vprintf Refined : "AutG computed: %o s\n", t1-t0;
  t0 := Cputime();
  AutGFP, mGFP := FPGroup(AutG);
  OutG, mOutG := OuterFPGroup(AutG);
  OutGperm, mOutGperm := PermutationGroup(OutG);
  OutAuts := [mGFP(g@@mOutGperm@@mOutG) : g in OutGperm];
  OutAuts := [outtau : outtau in OutAuts | [C(outtau(si)) : si in sigma] eq Csigma];
  t1 := Cputime();
  vprintf Refined : "outer auts computed: %o s\n", t1-t0;
  // identify via outer auts
  sigmapsAll := sigmaps; // just to save them...
  j := #sigmaps;
  vprintf Refined : "checking %o triples:\n", #sigmaps;
  while j gt 1 do
    foundj := false;
    for outtau in OutAuts do
      outtausigmapj := [outtau(sip) : sip in sigmaps[j]];
      for i := 1 to j-1 do
        if IsConjugate(G, sigmaps[i], outtausigmapj) then
          /* vprintf Refined : "\tremoving triple %o\n", j; */
          Remove(~sigmaps, j);
          foundj := true;
          break;
        end if;
      end for;
      if foundj then
        break;
      end if;
    end for;
    if not foundj then
      vprintf Refined : "\tkeeping triple %o\n", j;
    end if;
    j -:= 1;
  end while;
  // return
  vprintf Refined : "#triples before outer aut identification = %o\n", #sigmapsAll;
  vprintf Refined : "#triples after  outer aut identification = %o\n", #sigmaps;
  return sigmaps, sigmapsAll;
end intrinsic;

intrinsic RefinedPassportLemma(G::GrpPerm, tau::SeqEnum[GrpPermElt] : generates := true) -> Any
  {}
  /* N := NormalizerInAmbient(G); */
  N := G;
  S := SymmetricGroup(Degree(G));
  nus := DoubleCosetRepresentatives(N, Centralizer(N, tau[1]), Centralizer(N, tau[2]));
  pairs := [ [ G ! tau[1], G ! (nu*tau[2]*nu^(-1)) ] : nu in nus ];
  if generates then
      pairs := [ pair : pair in pairs | sub< S | pair > eq G ];
  end if;
  return pairs;
end intrinsic;

intrinsic RefinedPassport(sigma::SeqEnum[GrpPermElt]) -> Any
  {}
  // setup
  if #sigma eq 2 then
    Append(~sigma, (sigma[2]*sigma[1])^-1);
  end if;
  assert #sigma eq 3;
  d := Degree(Parent(sigma[1]));
  G := sub<Sym(d)|sigma>;
  t0 := Cputime();
  C := ClassMap(G);
  t1 := Cputime();
  vprintf Refined : "Class map computed: %o s\n", t1-t0;
  Csigma := [C(si) : si in sigma];
  vprintf Refined : "Classes of sigma: %o\n", Csigma;
  // RefinedPassportLemma
  tau := sigma[1..2];
  t0 := Cputime();
  passlemma := RefinedPassportLemma(G, tau : generates := true);
  t1 := Cputime();
  vprintf Refined : "Passport lemma: %o s\n", t1-t0;
  vprintf Refined : "#passlemma = %o\n", #passlemma;
  t0 := Cputime();
  passlemma := [sigmap : sigmap in passlemma | [C(sip) : sip in sigmap] eq Csigma[1..2]];
  t1 := Cputime();
  vprintf Refined : "#passlemma matching %o = %o\n", Csigma[1..2], #passlemma;
  // extend pairs to triples
  sigmaps := [];
  for sigmap in passlemma do
    sigmapoo := (sigmap[2]*sigmap[1])^-1;
    if C(sigmapoo) eq Csigma[3] then
      Append(~sigmaps, sigmap cat [sigmapoo]);
    end if;
  end for;
  vprintf Refined : "#triples matching %o = %o\n", Csigma, #sigmaps;
  // outer auts
  t0 := Cputime();
  AutG := AutomorphismGroup(G);
  t1 := Cputime();
  vprintf Refined : "AutG computed: %o s\n", t1-t0;
  t0 := Cputime();
  AutGFP, mGFP := FPGroup(AutG);
  OutG, mOutG := OuterFPGroup(AutG);
  OutGperm, mOutGperm := PermutationGroup(OutG);
  OutAuts := [mGFP(g@@mOutGperm@@mOutG) : g in OutGperm];
  OutAuts := [outtau : outtau in OutAuts | [C(outtau(si)) : si in sigma] eq Csigma];
  t1 := Cputime();
  vprintf Refined : "outer auts computed: %o s\n", t1-t0;
  // identify via outer auts
  sigmapsAll := sigmaps; // just to save them...
  j := #sigmaps;
  vprintf Refined : "checking %o triples:\n", #sigmaps;
  while j gt 1 do
    foundj := false;
    for outtau in OutAuts do
      outtausigmapj := [outtau(sip) : sip in sigmaps[j]];
      for i := 1 to j-1 do
        if IsConjugate(G, sigmaps[i], outtausigmapj) then
          /* vprintf Refined : "\tremoving triple %o\n", j; */
          Remove(~sigmaps, j);
          foundj := true;
          break;
        end if;
      end for;
      if foundj then
        break;
      end if;
    end for;
    if not foundj then
      vprintf Refined : "\tkeeping triple %o\n", j;
    end if;
    j -:= 1;
  end while;
  // return
  vprintf Refined : "#triples before outer aut identification = %o\n", #sigmapsAll;
  vprintf Refined : "#triples after  outer aut identification = %o\n", #sigmaps;
  return sigmaps, sigmapsAll;
end intrinsic;

intrinsic IsEquivalent(aut::GrpAutoElt, sigma::SeqEnum[GrpPermElt], sigmap::SeqEnum[GrpPermElt]) -> BoolElt
  {}
  assert #sigma eq 3 and #sigmap eq 3;
  if not aut(sigma[1]) eq sigmap[1] then
    return false;
  else
    if not aut(sigma[2]) eq sigmap[2] then
      return false;
    else
      if not aut(sigma[3]) eq sigmap[3] then
        return false;
      else
        return true;
      end if;
    end if;
  end if;
end intrinsic;

intrinsic IsEquivalent(AutG::GrpAuto, sigma::SeqEnum[GrpPermElt], sigmap::SeqEnum[GrpPermElt]) -> BoolElt
  {}
  for aut in AutG do
    if IsEquivalent(aut, sigma, sigmap) then
      return true;
    end if;
  end for;
  return false;
end intrinsic;

intrinsic RefinedPassportBrutal(sigma::SeqEnum[GrpPermElt]) -> Any
  {}
  // setup
  if #sigma eq 2 then
    Append(~sigma, (sigma[2]*sigma[1])^-1);
  end if;
  assert #sigma eq 3;
  d := Degree(Parent(sigma[1]));
  G := sub<Sym(d)|sigma>;
  t0 := Cputime();
  C := ClassMap(G);
  t1 := Cputime();
  vprintf Refined : "Class map computed: %o s\n", t1-t0;
  Csigma := [C(si) : si in sigma];
  vprintf Refined : "Classes of sigma: %o\n", Csigma;
  Cs := [Class(G, si) : si in sigma];
  vprintf Refined : "#G   = %o\n", #G;
  vprintf Refined : "#c0  = %o\n", #Cs[1];
  vprintf Refined : "#c1  = %o\n", #Cs[2];
  vprintf Refined : "#coo = %o\n", #Cs[3];
  // build triples brutally
  triples := [];
  c0, c1, coo := Explode(Cs);
  for g0 in c0 do
    t0 := Cputime();
    for g1 in c1 do
      goo := (g1*g0)^-1;
      if goo in coo then
        if sub<Sym(d)|[g0,g1,goo]> eq G then
          Append(~triples, [g0,g1,goo]);
        end if;
      end if;
    end for;
    t1 := Cputime();
    vprintf Refined : "i=%o out of %o: %o s\n", Index(SetToSequence(c0), g0), #c0, t1-t0;
  end for;
  // automorphisms of G
  AutG := AutomorphismGroup(G);
  mp, AutGperm := PermutationRepresentation(AutG);
  assert #AutG eq #AutGperm;
  auts := [aut@@mp : aut in AutGperm];
  // pop the stack
  all := triples; // just to save
  vprintf Refined : "checking %o triples:\n", #triples;
  vprintf Refined : "#AutG = %o\n", #AutG;
  j := #triples;
  while j gt 1 do
    foundj := false;
    for aut in auts do
      for i := 1 to j-1 do
        if IsEquivalent(aut, triples[i], triples[j]) then
          Remove(~triples, j);
          foundj := true;
          break;
        end if;
      end for;
      if foundj then
        break;
      end if;
    end for;
    if not foundj then
      vprintf Refined : "\tkeeping triple %o\n", j;
    end if;
    j -:= 1;
  end while;
  return triples;
end intrinsic;

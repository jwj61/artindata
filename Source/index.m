// Make LMFDB Artin database, TD Apr 2013
// Modified JJ Nov 2015


load mmylib;

//Masterdir:="F:\\LMFDB\\Source";
//LMFDBdir:="F:\\LMFDB\\";
Masterdir:="/home/jj/lmfdb/art3/Source";
LMFDBdir:="/home/jj/lmfdb/art3/";

DATAdir:=LMFDBdir*"DATA/";
SrcNFDBdir:=LMFDBdir*"NFDB/";
errfile:=LMFDBdir*"errors.txt";
skipfile:=LMFDBdir*"skip.txt";
DBfile:=LMFDBdir*"artin.dat";
NFgalfile:=LMFDBdir*"nfgal.json";
NFartfile:=LMFDBdir*"artrep.json";
PGLdir:=LMFDBdir*"pgl/";
GL2file:=PGLdir*"gl2.dat";
SD16file:=PGLdir*"sd16.dat";
Wreathdir:=LMFDBdir*"wreath/";
Smalldir:=LMFDBdir*"small/";
Familiesdir:=LMFDBdir*"families/";
MalleKluenersdir:=LMFDBdir*"trans2d/";
OldArtdir:=LMFDBdir*"oldart/";
artlabeldir:=LMFDBdir*"artlabels/";
MaxFrobP:=1000;
FrobHigh:=10000;
FrobHighPrimes:=PrimesUpTo(FrobHigh);
computesign:=true;
PrintFound:=false;
MaxLCf:=10000;
wrcon:=false;


load jjlib;

load galpols;
R<x>:=PR(Q);

try 
  SkipPolys:=[eval d: d in Split(Read(skipfile))];
catch e 
  SkipPolys:=[];
end try;

function Indices2File(dir,file)
  dir:=Sprint(dir);
  file:=Sprint(file);
  if #dir eq 1 then dir:="0"*dir; end if;
  if #file eq 1 then file:="0"*file; end if;
  return DATAdir*dir*"/"*file*".dat";
end function;


function NF2File(K)
  if Type(K) eq RngUPolElt then
    K:=NumberField(K); 
  end if;
  O:=IntegerRing(K);
  D:=Discriminant(O) mod 9973;
  dir:=D div 100;
  file:=D mod 100;
  return Indices2File(dir,file);
end function;


procedure EraseData()
  "Erasing the database";
  System("rmdir /S /Q "*Prune(DATAdir));
  System("mkdir "*Prune(DATAdir));
end procedure;


procedure MakeDirs()
  "Creating 00..99 directories in data/";
  for i:=0 to 99 do
    s:=Sprint(i);
    if #s eq 1 then s:="0"*s; end if;   
    System("mkdir "*DATAdir*s);  
  end for;
end procedure;


function SortPolyList(list)
  list:=SetToSequence(Set(list));
  D:=[[Degree(f),Max([Abs(c): c in Eltseq(f)]),#[x: x in Eltseq(f) | x ne 0]]: f in list];
  Sort(~D,~P);
  list:=PermuteSequence(list,P);
  return list;
end function;



function FrobResData(d,c)
  if Type(d) eq RngIntElt then
    s:=[* c,"\"CYC\"",0,d *];
  elif Type(d[3]) eq Tup then
    s:=[* c,"\"RES\"",d[3],d[4] *];
  else
    s:=[*c,"\"ALT\"",DelSpaces(Eltseq(d[3])),d[4] *];
  end if;
  return s;
end function;


function ArtinLocalData(loc)
  rts, DSn, Frob, ramgps, DtoG, Zpprec:=Explode(loc);
  K:=Parent(rts[1]);
  return <Characteristic(ResidueClassField(IntegerRing(K))),
          #DSn,
          #ramgps[1][1],
          [Eltseq(x): x in MinimizedGenerators(DSn)],
          [Eltseq(c[3]): c in ConjugacyClasses(DSn)],
          Eltseq(Frob),
          DtoG,
          [r[2]: r in ramgps],
          [[Eltseq(g): g in MinimizedGenerators(r[1])]: r in ramgps],         
          Zpprec,
          PrintRelExtEquation(K),
          [PrintRelExtElement(r): r in rts]
          >;
end function;


procedure PolyIndex(f,~I,~j,force)
  cfs:=[Eltseq(cf): cf in Coefficients(f)];
  if not force then
    j:=Position(I,cfs);
    if j ne 0 then return; end if;
  end if;  
  Append(~I,cfs);
  j:=#I;
  return;
end procedure;


function FieldExists(F)
  filename:=NF2File(F);
  D:=Discriminant(IntegerRing(F));
  poly:=DefiningPolynomial(F);
  found:=0;
  fdata:=0;
  try
    for data in Split(Read(filename)) do
      f,Df:=Explode(Split(data," "));
      f:=eval f;
      Df:=eval Df;
      // if (D eq Df) and ((D eq 1) or (f eq poly) or IsIsomorphic(F,NumberField(f))) then
      if (D eq Df) and ((D eq 1) or (f eq poly)) then
        found:=f; fdata:=data; break;
      end if;
    end for;
  catch e
    ; // e;     
  end try;
  if (found ne 0) and PrintFound then
    Sprintf("%o: %o -> %o",
      filename,
      DelSpaces(R!DefiningPolynomial(F)), 
      found eq 0 select "not found" else DelSpaces(f)
    );
  end if;
  return found,fdata;
end function;


function SmallPrimeDivisors(N)
  fct:=Factorisation(Z!N: TrialDivisionLimit:=50000,
    PollardRhoLimit:=200000, MPQSLimit:=0, ECMLimit:=20, Proof:=false);
  P:={f[1]: f in fct};
  return P;
end function;


function HardPrimes(F)
  P:=SmallPrimeDivisors(F`artinrepdata`disc);
  for d in F`artinrepdata`Inv do
    if Type(d) eq RngIntElt then continue; end if;
    fs:=[u[1]: u in d[4]];
    for i:=2 to #fs do
    f:=fs[i];
    if Degree(f) gt 100 then continue; end if;  
    for j:=1 to i-1 do
      g:=fs[j];       
      if Degree(g) gt 100 then continue; end if;  
      res:=Resultant(f,g);
      P:=P join SmallPrimeDivisors(res);
    end for;
    end for;
  end for;
  P:=Sort([p: p in P | p lt 10^50]);
  return P;
end function;


procedure IndexDatabase()
  sizes:=[];
  list:=[];
  ArtArray:=AssociativeArray();
// This will have poly coeffs as string as key and value is a list
// of artrep labels
  "INDEXING THE DATABASE";
  "Reading...";
  for dir,file in [0..99] do
    filename:=Indices2File(dir,file);   
    try 
      data:=Split(Read(filename)); 
    catch e;
      data:=[];
    end try;
    for i:=1 to #data do
      s:=data[i];
      f,D,s1,s2,GSize,Gname,Ggens,CC,cycs,p,QprM,Qprprec,r,complexconjugation,
        Frobs,FrobRes,localdata,ArtKernels,ArtIndices,gorbs,ArtReps
          :=Explode(Split(s," "));
      f:=eval f;
      nfkey := Sprint(Coefficients(f));
      ArtReps := eval ArtReps;
      ArtReps := <z[10] : z in ArtReps>; // Label pairs
      ArtArray[nfkey] := ArtReps;
      D:=eval D;
      GSize:=eval GSize;
      size:=[GSize,Degree(f),Abs(D)];
      Append(~sizes,size);
      Append(~list,[dir,file,i]);
    end for;
  end for;
  "Total entries:",#sizes;
  "Sorting...";
  Sort(~sizes,~P);
  list:=PermuteSequence(list,P);
  "Dumping...";
  Rewrite(DBfile);
  for d in list do
    dir,file,i:=Explode(d);
    filename:=Indices2File(dir,file);   
    data:=Split(Read(filename))[i]; 
    f,D,s1,s2,GSize,Gname,Ggens,CC,cycs,p,QprM,Qprprec,r,complexconjugation,
      Frobs,FrobRes,localdata,ArtKernels,ArtIndices,gorbs,ArtReps
        :=Explode(Split(data," "));
    ArtReps:= eval ArtReps;
    ArtKernels:= eval ArtKernels;
    ArtKernels:= < z eq 0 select Sprint(Coefficients(eval f)) else Sprint(Coefficients(z)) : z in ArtKernels>;
    ArtIndices:= eval ArtIndices;
    mylabels := <ArtArray[ArtKernels[j]][ArtIndices[j]] : j in [1..#ArtKernels]>;
    //data:=&cat[c eq "\"" select "'" else c: c in Eltseq(data)];
    "Got ",mylabels;
    write(DBfile,data: con:=false);        
  end for;
  "Indexing done";
end procedure;


function PruneZeroes(c)
  assert Type(c) eq SeqEnum;
  nz:=[j: j in [1..#c] | c[j] ne 0];
  if IsEmpty(nz) then 
    return [0]; 
  elif nz[#nz] eq #c then 
    return c;
  else 
    return c[[1..nz[#nz]]];
  end if;
end function;

function w(val: q:=true)
  if Type(val) eq MonStgElt
    then s:=(q select "\"" else "")*val*(q select "\"" else "");
  elif Type(val) eq SeqEnum 
    then if IsEmpty(val) 
      then s:="[]"; 
      else s:="["*&*[w(val[i])*(i eq #val select "" else ","): i in [1..#val]]*"]";
    end if;
  elif (Type(val) eq Tup) and (#val eq 3) and (Type(val[1]) eq RngIntElt) 
    and (Type(val[2]) eq RngIntElt) and (Type(val[3]) eq SeqEnum) then
    s:=Sprintf("{\"Order\":%o, \"Size\":\"%o\", \"Representative\": %o}",
      val[1],val[2],w(Eltseq(val[3])));
  elif (Type(val) eq Tup) and (#val eq 5) and (Type(val[1]) eq RngIntElt)
    and (Type(val[2]) eq RngIntElt) and (Type(val[3]) eq RngIntElt) and
    (Type(val[4]) eq RngIntElt) and (Type(val[5]) eq SeqEnum) then
    s:=Sprintf("{\"Dim\":%o, \"Conductor\":\"%o\", \"DBIndex\":%o, \"CharacterField\":%o, \"Character\":%o}",
      val[1],val[2],val[3],val[4],w([PruneZeroes(c): c in val[5]]));
  elif (Type(val) eq Tup) and (#val eq 2) and (Type(val[1]) eq RngUPolElt)
    and (Type(val[2]) eq RngIntElt) then
    s:=Sprintf("{\"RootOf\":%o, \"ConjugacyClass\":%o}",
      //w(ReplaceStringFunc(DelSpaces(val[1]),"a","x")),val[2]);
      w([Sprint(c): c in Coefficients(val[1])]),val[2]);
  elif (Type(val) eq Tup) and (#val eq 2) and (Type(val[1]) eq SeqEnum)
    and (Type(val[2]) eq RngUPolElt) then
    s:=Sprintf("{\"Powers\":%o, \"Resolvent\":%o}",
      w(val[1]),
      //w(ReplaceStringFunc(DelSpaces(val[2]),"a","x")));
      w(Coefficients(val[2])));
  elif Type(val) eq AlgChtrElt then
    s:=w([PruneZeroes(Eltseq(x)): x in Eltseq(val)]);
  elif (Type(val) eq List) and (#val eq 4) and (Type(val[2]) eq MonStgElt) then
    s:=Sprintf("{\"CycleType\":%o, \"Algorithm\":\"%o\", \"Data\":%o, \"Classes\":%o}",
      w(val[1]),val[2],w(val[3]),w(val[4]));
  else s:=Sprint(val);
  end if;
  return s; 
end function;


AddAttribute(FldRat,"outg");
AddAttribute(FldRat,"outa");
Q`outg:="";
Q`outa:="";


procedure wkg(key,val: last:=false, q:=true, fin:=false)
  /* if key cmpeq "" 
    then write(NFgalfile,val: con:=wrcon);
    else write(NFgalfile,w(key: q:=true)*":"*w(val: q:=q)*(last select "" else ","): con:=wrcon);
  end if; */

  Q:=Rationals();
  if key cmpeq "" 
    then Q`outg*:=Sprint(val)*"\n";
    else Q`outg*:=w(key: q:=true)*":"*w(val: q:=q)*(last select "" else ",")*"\n";
  end if;

  if (last and (Random(50) eq 0)) or fin then
    write(NFgalfile,Prune(Q`outg): con:=wrcon);
    Q`outg:="";
  end if;
end procedure;


procedure wka(key,val: last:=false, q:=true, fin:=false)
  Q:=Rationals();
  if key cmpeq "" 
    then Q`outa*:=Sprint(val)*"\n";
    else Q`outa*:=w(key: q:=true)*":"*w(val: q:=q)*(last select "" else ",")*"\n";
  end if;
  if (last and (Random(50) eq 0)) or fin then
    write(NFartfile,Prune(Q`outa): con:=wrcon);
    Q`outa:="";
  end if;

  /* if key cmpeq ""
    then write(NFartfile,val: con:=wrcon);
    else write(NFartfile,w(key: q:=true)*":"*w(val: q:=q)*(last select "" else ","): con:=wrcon);
  end if; */
end procedure;


procedure ExportToJSON(: test:=false)
  C<i>:=ComplexField(9);
  "Erasing nfgal.json and artrep.json";
  Rewrite(NFgalfile);
  Rewrite(NFartfile);
  "Reading database";
  DATA:=Split(Read(DBfile)); 
  "Processing";
  fs:=[];
  NFind:=[];
  NFpolys:=[]; // will hold list of polys
  NFart:=[]; // will hold parallel list of artin reps
  ARind:=[];
  wkg("","[");
  wka("","[");
  for j:=1 to #DATA do
    if test and (j gt 100) then break; end if;
    data:=DATA[j];
    f,D,r1,r2,Gsize,Gname,Ggens,CC,cycs,p,QprM,Qprprec,r,complexconjugation,
      Frobs,FrobRes,localdata,ArtKernels,ArtIndices,ArtReps
        :=Explode(Split(data," "));
    if j mod 100 eq 1 then Sprintf("%o/%o %o",j,#DATA,f); end if;
    f:=eval f;
    Gsize:=eval Gsize;
    Append(~NFpolys,f);
    ArtReps:=eval ArtReps;     
    ArtKernels:=eval ArtKernels;     
    ArtIndices:=eval ArtIndices;     
    alist:=[];
    //NFdata:=[Degree(f),Gsize];
    //Append(~NFind,NFdata);
    //nfindex:=Count(NFind,NFdata);
// Need to group Galois orbits
// explode first one and pull it out?
    indices:=[];
    for ai:=1 to #ArtReps do
      adata:=ArtReps[ai];
      dim,cond,n,char,I,badpr,bads,ind,lvals,alabel,galnt:=Explode(adata);
      sign,lval1,lval2:=Explode([C!x: x in lvals]);
      if ind eq 1 
        then sign:=1;
        elif ind eq -1 then sign:=Round(Real(sign)); 
        else sign:=0;
      end if;
      anf:=ArtKernels[ai]; // 0 or poly for smaller field
      if anf eq 0 then
        ad:=[dim,cond];
        Append(~ARind,ad);
        index:=Count(ARind,ad);
        wka("",Degree(f) eq 1 select "{" else ",{");
        wka("Dim",dim);
        wka("Conductor",Sprint(cond));
        // label instead
        // split the label?
        wka("Label",alabel);
        //wka("DBIndex",index);
        wka("NFGal",Coefficients(f));
        //wka("NFGal",[Sprint(x): x in [Degree(f),Gsize,nfindex]]);
        wka("Galois_nt", galnt);
        wka("CharacterField",n);
        wka("Character",[PruneZeroes(c): c in char]);
        wka("LocalFactors",[[PruneZeroes(c): c in f]: f in I]);
        wka("BadPrimes",[Sprint(p): p in PrimeDivisors(cond)]);
        wka("HardPrimes",[Sprint(p): p in badpr]);
        wka("HardFactors",bads);
        wka("Sign",sign);
        wka("Indicator",ind: last:=true);        
        wka("","}");
      else // ref smaller field
        u:=Position(NFpolys,anf);
        error if u eq 0, "Field "*DelSpaces(anf)* " not found, while processing "*
          DelSpaces(f);
        // Turn this into a label
        index:=NFart[u][ArtIndices[ai]];
      end if;                
      Append(~alist,<dim,cond,index,n,char>);
      Append(~indices,index);
    end for;
    Append(~NFart,indices);
    wkg("",Degree(f) eq 1 select "{" else ",{");
    wkg("TransitiveDegree",Degree(f));
    wkg("Polynomial",Coefficients(f));
    wkg("Size",Sprint(Gsize));
    //wkg("DBIndex",nfindex);
    wkg("G-Gens",Ggens: q:=false);
    wkg("G-Name",Gname: q:=false);
    wkg("QpRts-p",p: q:=false);
    wkg("QpRts-minpoly",Coefficients(eval QprM));
    wkg("QpRts-prec",Qprprec: q:=false);
    wkg("QpRts",[[Sprint(c): c in u]: u in eval r]: q:=false);
    wkg("ConjClasses",eval CC);
    wkg("ComplexConjugation",complexconjugation: q:=false);
    wkg("Frobs",Frobs: q:=false);
    wkg("ArtinReps",alist);
    //wkg("FrobResolvents",FrobRes: last:=true);
    wkg("FrobResolvents",eval FrobRes: last:=true);
    wkg("","}");

   
/*    
"TransitiveDegree":4,
"Polynomial":[1,1,1,1,1],
"Size":"4",
"DBIndex":41,
"G-Gens":[[4,3,1,2],[2,1,4,3]],
"G-Name":"C4",
"QpRts-p":11,
"QpRts-minpoly":[-1,1],
"QpRts-prec":5,
"QpRts":["46709","37107","-56601","-27216"],
"ConjClasses":[{"Order":1, "Size":"1", "Representative": [1,2,3,4]},{"Order":2, "Size":"1", "Representative": [2,1,4,3]},{"Order":4, "Size":"1", "Representative": [4,3,1,2]},{"Order":4, "Size":"1", "Representative": [3,4,2,1]}],
"ComplexConjugation":2,
"Frobs":[3,4,1,3,1,4,3,2,4,2,1,3,1,4,3,4,2,1,3,1,4,2,4,2,3,1,4,3,2,4,3,1,3,2,2,1,3,4,3,4,2,1,1,4,3,2,1,4,3,2,4,2,1,1,3,4,2,1,3,1,4,4,3,1,4,3,1,3,3,2,4,2,3,4,2,4,2,3,1,2,2,1,1,4,2,4,2,3,1,4,3,2,3,1,2,4,2,1,4,1,3,3,4,2,1,3,3,4,2,1,3,4,3,2,1,1,4,3,4,2,1,4,3,4,1,1,2,2,3,4,2,4,1,3,1,2,4,3,3,2,1,1,4,3,2,2,4,3,2,4,3,1,4,3,3,1,2,2,3,1,3,4,3,1,3,4,1,3],
"ArtinReps":[{"Dim":1, "Conductor":"1", "DBIndex":1, "CharacterField":1, 
  "Character":[[1],[1],[1],[1]]},
  {"Dim":1, "Conductor":"5", "DBIndex":2, "CharacterField":4, 
  "Character":[[1],[-1],[0,1],[0,-1]]},
  {"Dim":1, "Conductor":"5", "DBIndex":1, "CharacterField":1, 
  "Character":[[1],[1],[-1],[-1]]},
  {"Dim":1, "Conductor":"5", "DBIndex":3, "CharacterField":4, 
  "Character":[[1],[-1],[0,-1],[0,1]]}],
"FrobResolvents":[{"CycleType":[1,1,1,1], "Algorithm":"CYC", "Data":0, "Classes":1},
   {"CycleType":[2,2], "Algorithm":"CYC", "Data":0, "Classes":2},
   {"CycleType":[4], "Algorithm":"RES", "Data":{"Powers":[1], "Resolvent":[0,0,1]},
   "Classes":[{"RootOf":["1","1"], "ConjugacyClass":3},
     {"RootOf":["-4","1"], "ConjugacyClass":4}]}]
*/
    
  
/*
    f,                                       // Defining polynomial
    PrintFactor(ODisc),                      // Discriminant
    sign1,                                   // r1=#real embeddings
    sign2,                                   // r2=#complex pairs of embeddings
    #G,                                      // Size of G
    GroupName(G),                            // Name of G
    [Eltseq(g): g in Generators(G)],         // Generators of G
    [<c[1],c[2],Eltseq(c[3])>: c in CC],     // CC (Conjugacy Classes)
    D`cycs,                                  // cycle types
    Characteristic(ResidueClassField(Qpr)),  // p
    DelSpaces(QprM),                         // Min poly of Qp unr ext
    Min([Precision(r): r in D`r]),           // p-adic precision
    [PrintRelExtElement(x): x in D`r],          // roots
    gToCC(FrobeniusElement(F,Infinity())),   // CC of complex conjugation
    [gToCC(FrobeniusElement(F,p)): p in PrimesUpTo(MaxFrobP)],
                                             // Frobenius elements
    [FrobResData(D`Inv[i],D`cycs[i]): i in [1..#D`Inv]],
                                             // Frobenius resolvents
    [*ArtinLocalData(loc): loc in D`localdata*],   // Ramification data
    ArtKernels,                              // Indices of base fields of ArtReps
    ArtIndices,                              // and of the representations in
                                             //    their character tables
    [ArtinRepData(a): a in A],               // Artin Representations themselves
*/

    Append(~fs,f);
  end for;
  wkg("","]": fin:=true);
  wka("","]": fin:=true);
  "ExportToJSON done";
end procedure;


procedure PGLExtensions(d: MinN:=11, MaxN:=100, prec:=0)
  filename:=PGLdir*Sprint(d)*".dat";
  Rewrite(filename);
  R<x>:=PR(Q);
  D:=CremonaDatabase();
  list:=[];
  for n:=MinN to MaxN do
  for E in EllipticCurves(D,n) do
    f:=PGLPolynomial(E,d);
    den:=Denominator(VectorContent(Coefficients(f)));
    f*:=den;
    if not IsIrreducible(f) then continue; end if;
    K:=NumberField(f);
    if exists{F: F in list | IsIsomorphic(K,F)} then continue; end if;
    Append(~list,K);

    fred:=PolRedAbs(f);
    
    /*
    O:=IntegerRing(K);
    O:=LLL(O);
    pols:=[MinimalPolynomial(x): x in Basis(O)];
    pols:=[f: f in pols | Degree(f) eq AbsoluteDegree(K)];
    if IsEmpty(pols) then continue; end if;
    B:=[Max([Abs(c): c in Eltseq(f)]): f in pols];  
    mi,j:=Min(B);
    fred:=pols[j];
    */
    
    CremonaReference(E),d,GroupName(GaloisGroup(f)),DelSpaces(fred);
    write(filename,DelSpaces(fred): con:=false);
  end for;
  end for;
end procedure;


procedure CreatePGLExtensions()
  //PGLExtensions(3: MaxN:=200);
  //PGLExtensions(4: MaxN:=200);
  //PGLExtensions(5: MaxN:=200);
  //PGLExtensions(6: MaxN:=100);  // does not work, even for 11a1
  //PGLExtensions(7: MaxN:=100);
  //PGLExtensions(8: MaxN:=100);
  //PGLExtensions(9: MaxN:=100);
  //PGLExtensions(10: MaxN:=100);   // does not work, even for 11a1
  //PGLExtensions(11: MaxN:=100);
end procedure;



procedure CreateSmallCoefficients(d,B)
  maxb:=0;
  for b:=1 to 50 do
    filename:=Sprintf("%od%oB%o.dat",Smalldir,d,b);
    ok,R:=ReadTest(filename);
    if ok and (b ge B) then   
      Sprintf("d=%o B=%o exists",d,b); return;
    end if;
    if ok and (b lt B) then 
      maxb:=b;
    end if;
  end for;
  full:=Factorial(d);
  S:=RSpace(Integers(2*B+1),d);
  "SmallCoefficients:",#S,"polynomials";
  list:=[];
  for H in [maxb+1..B] do
  for x in S do
    cfs:=[Z!c-B: c in Eltseq(x)] cat [1];
    max:=Max([Abs(c): c in cfs]);
    if max ne H then continue; end if;
    f:=Polynomial(Q,cfs);
    if not IsIrreducible(f) then continue; end if;
    G:=GaloisGroup(f);
    if #G eq full then continue; end if;
    GroupName(G),DelSpaces(f);
    Append(~list,f);
  end for;
  end for;
  filename:=Sprintf("%od%oB%o.dat",Smalldir,d,B);
  "writing",filename;
  Rewrite(filename);
  for f in list do
    write(filename, DelSpaces(f): con:=false);
  end for;
end procedure;



function Shorts(L,max)
  d:=2;
  repeat
    S:=ShortVectors(L,d);
    if #S ge max then return S; end if;
    d:=Round(1.25*d+1);
  until false;
end function; 


procedure GenerateWreath(fpol,gpol,MaxD,MaxV,minrepallowed)
  filename:=Sprintf("%o[%o][%o].dat",Wreathdir,fpol,gpol); 
  ok:=ReadTest(filename);
  if ok then filename,"exists"; return; end if;
  for d in [-MaxD..MaxD] do
    if (d in [0,1]) or not IsSquarefree(d) then continue; end if; 
    f:=eval fpol;
    if not IsIrreducible(f) then continue; end if;
    K:=NumberField(f);
    U<y>:=PR(K);
    O:=LLL(IntegerRing(K));
    S:=Shorts(Lattice(O),MaxV);
    V:=[K![Round(c): c in Eltseq(s[1])]: s in S];
    for v in V do
      g:=eval gpol;
      if not IsIrreducible(g) then continue; end if;
      if Degree(g)*Degree(f) gt 16 then continue; end if; 
      poly0:=MinimalPolynomial(AbsoluteField(NumberField(g)).1);
      poly:=PolRedAbs(poly0);    
      G:=GaloisGroup(poly);
      minrep:=Min([10^10] cat [Degree(c): c in CharacterTable(G) | IsFaithful(c)]);
      GroupName(G),minrep,DelSpaces(poly);
      if minrep eq 10^10 then continue; end if;
      if minrep gt minrepallowed then continue; end if;  
      write(filename,poly: con:=false);
    end for;
  end for;
end procedure;


procedure CreateWreathProducts()
  FList:=["x^2-d","x^3-d","x^3+x-d","x^4-d","x^4-x+d"];
  GList:=["y^2-v","y^3-v","y^3+y-v","y^4-v","y^4-v+d"];
  for f in FList, g in GList do
    GenerateWreath(f,g,5,10,8);
  end for;
  GenerateWreath("x^5-d","y^2-v",5,10,8);
  GenerateWreath("x^2-d","y^5-v",5,10,8);
  GenerateWreath("x^6-d","y^2-v",5,10,8);
  GenerateWreath("x^2-d","y^6-v",5,10,8);
  GenerateWreath("x^2-d","y^7-v",5,10,8);
  GenerateWreath("x^2-d","y^8-v",5,10,8);
end procedure;



procedure CreateGL2Extensions()
  if ReadTest(GL2file) then
    GL2file,"exists"; return;
  end if;
  Rewrite(GL2file); 
  D:=CremonaDatabase();
  list:=[];
  for n:=1 to 100 do
  for E in EllipticCurves(D,n) do
    f:=TorsionPolynomial(E,3);
    if not IsIrreducible(f) then continue; end if;
    f:=PolRedAbs(f);
    if f in list then continue; end if;
    Append(~list,f);  
    CremonaReference(E),GroupName(GaloisGroup(f)),DelSpaces(f);
    write(GL2file,DelSpaces(f): con:=false);
  end for;
  end for;
end procedure;


procedure CreateSD16Extensions(:overwrite:=false)
  if not overwrite and ReadTest(SD16file) then
    SD16file,"exists"; return;
  end if;
  Rewrite(SD16file); 
  D:=CremonaDatabase();
  list:=[];
  MaxDisc:=2000000;
  for den in [1..10] do
  for num in [-200..200] do
    if GCD(num,den) ne 1 then continue; end if;
    j:=(num/den)^3;
    E:=MinimalModel(EllipticCurveFromjInvariant(j));
    D:=Z!Discriminant(E);
    twist:=&*[Z|p: p in PrimeDivisors(D) | (p ne 2) and (den mod p ne 0) and Valuation(D,p) eq 6];
    for t in [twist,-twist,2*twist,-2*twist] do
      ET:=MinimalModel(QuadraticTwist(E,t));
      //if Abs(Discriminant(ET)) lt MaxDisc then
      if &*(PrimeDivisors(Z!Discriminant(ET))) lt 100 then
        f:=TorsionPolynomial(ET,3);
        if not IsIrreducible(f) then continue; end if;
        if not (LeadingCoefficient(f) in [1,27]) then continue; end if;
        f:=PolRedAbs(f);
        if f in list then continue; end if;
        Append(~list,f);  
        DelSpaces(aInvariants(ET)),PrintFactor(Discriminant(ET)),DelSpaces(f);
        write(SD16file,DelSpaces(f): con:=false);
      end if;
    end for;
  end for;
  end for;
end procedure;


procedure CreateMalleKluenersFamilies()
// converts htmls from the Malle-Klueners website into dat files
  dir:=MalleKluenersdir;
  count:=0;
  for filename in FindFiles(dir,".html": full:=true) do
    filename;
    outfile:=ReplaceStringFunc(filename,".html",".dat");
    Rewrite(outfile);
    for s in Split(Read(filename),"<>\n") do
      if (#s eq 0) or (s[1] ne "x") then continue; end if;
      write(outfile,DelSpaces(eval s): con:=false);
      count+:=1;
    end for;
  end for;  
  "Total:",count;
end procedure;


procedure CreateFamily(Gname,H,MaxDisc)
  F:=FamilyWithGaloisGroup(Gname);
  if F cmpeq 0 then 
    "Family "*Gname*" not found"; return;
  end if;
  filename:=Familiesdir*Gname*".dat";
  Rewrite(filename);
  S:=SpecializeFamily(F,H);
  "Have",#S,"polynomials";
  total:=0;
  list:=[];
  for i:=1 to #S do
    f:=PolRedAbs(S[i]);
    if f in list then continue; end if;
    Append(~list,f);
    D:=Abs(Discriminant(f));
    size:=Round(Log(10,D));
    if size gt MaxDisc then continue; end if;
    total+:=1;
    guesstot:=Round(total/i*#S);
    Sprintf("%o %o %o (%o)",size,GroupName(GaloisGroup(f)),DelSpaces(f),guesstot);
    write(filename,DelSpaces(f): con:=false);
  end for;
  "Written",total,"polynomials";
end procedure;


procedure CreateFamilies()
  CreateFamily("D8",10,10);
  CreateFamily("D10",3,30);
  CreateFamily("D14",30,50);
  CreateFamily("T(6,3)",30,15);
  CreateFamily("T(6,5)",30,15);
  CreateFamily("T(8,13)",20,25);    
  CreateFamily("T(8,24)",20,15);    
  CreateFamily("T(8,37)",30,40);    
  CreateFamily("T(6,12)",25,20);  
  CreateFamily("T(8,5)",30,60);       //visit
  CreateFamily("T(8,6)",30,35);    
  CreateFamily("T(8,7)",30,80);       //visit
  CreateFamily("T(8,8)",15,130);   
  CreateFamily("T(8,11)",20,50);      //not found
  CreateFamily("T(8,12)",20,45);  
  CreateFamily("T(8,17)",25,30);  
  CreateFamily("T(8,23)",35,40);  
  CreateFamily("T(7,3)",30,90);       //visit
  CreateFamily("T(7,5)",25,40);  
  CreateFamily("T(8,49)",10,60);
  CreateFamily("D16",2,40);
  CreateFamily("PSL(2,11)",4,60);
  for p in PrimesUpTo(50) do
    if p mod 7 eq 1 then
      CreateFamily("27B"*Sprint(p),10,60);
    end if;
  end for;
end procedure;



function PrintLarge(n)
  if n lt 10^8 then return Sprint(n); end if;
  return ReplaceStringFunc(Sprint(RealField(2)!n),"E","\\!\\cdot\\!10^{")*"}";
end function;


procedure Statistics()

  tabfile:=Masterdir*"stats.inc";
  Rewrite(tabfile);

  timessym:="\\times ";
  
  C<i>:=ComplexField(10);

  filename:=LMFDBdir*"artsmall.dat";
  filename:=LMFDBdir*"artin.dat";

  "Statistics: processing",filename;

  None:=10^1000;
  MaxDim:=1000;
  S:=[{* <0,""> *}: dim in [1..MaxDim]];
  MinCond:=[None: dim in [1..MaxDim]];
  MaxCond:=[0: dim in [1..MaxDim]];
  for data in Split(Read(filename)) do
    f,D,r1,r2,Gsize,Gname,Ggens,CC,cycs,p,QprM,Qprprec,r,complexconjugation,
      Frobs,FrobRes,localdata,ArtKernels,ArtIndices,ArtReps :=Explode(Split(data," "));
    ArtKernels:=eval ArtKernels;     
    ArtIndices:=eval ArtIndices;     
    ArtReps:=eval ArtReps;
    for i:=1 to #ArtReps do     
      if ArtKernels[i] ne 0 then continue; end if;
      dim,cond,n,char,I,badpr,bads,ind,lvals:=Explode(ArtReps[i]);
      Include(~(S[dim]),<eval Gsize,Gname>);
      if MinCond[dim] gt cond then MinCond[dim]:=cond; end if;
      if MaxCond[dim] lt cond then MaxCond[dim]:=cond; end if;
    end for;
  end for;
  "Dumping";
  for dim:=1 to MaxDim do
    if MinCond[dim] eq None then continue; end if;
    write(tabfile,Sprintf("$d=%o\\quad N\\in[%o,%o]$\\\\",
      dim,PrintLarge(MinCond[dim]),PrintLarge(MaxCond[dim])
    ));
    X:=S[dim];
    XSet:=Sort(SetToSequence(Set(X)));
    for s in XSet do
      if s[2] eq "" then continue; end if;
      Gname:=ReplaceStringFunc(s[2],"\"","");
      try G:=GP(Gname); assert GroupName(G) eq Gname; Gname:=GroupName(G: tex:=true);
      catch e;
      end try;
      for n:=50 to 1 by -1 do
        ReplaceString(~Gname,"C"*Sprint(n),"C_{"*Sprint(n)*"}");
        ReplaceString(~Gname,"D"*Sprint(n),"D_{"*Sprint(n)*"}");
        ReplaceString(~Gname,"S"*Sprint(n),"S_{"*Sprint(n)*"}");
        ReplaceString(~Gname,"A"*Sprint(n),"A_{"*Sprint(n)*"}");
        ReplaceString(~Gname,"Q"*Sprint(n),"Q_{"*Sprint(n)*"}");
        ReplaceString(~Gname,"L(2,"*Sprint(n),"L_2(\\F_{"*Sprint(n)*"}");
        ReplaceString(~Gname,"L2(F"*Sprint(n),"L_2(\\F_{"*Sprint(n)*"}");
      end for;
      ReplaceString(~Gname,"G20","G_{20}");
      ReplaceString(~Gname,"\\wr","$");
      ReplaceString(~Gname,"wr","$");
      ReplaceString(~Gname,"$","\\wr ");
      ReplaceString(~Gname,"*","\\!\\times\\!");
      ReplaceString(~Gname,":","$");ReplaceString(~Gname,"$","\\!:\\!");
      write(tabfile,Sprintf("\\hbox{$%o\\,(%o)$}%o ",
        Gname,Multiplicity(X,s),Last(XSet) eq s select "\\par\\smallskip" else ",\\ ")
      );
    end for;
  end for;
end procedure;



// run EraseData() to completely Erase the database and directory structure
// run MakeDirs() to create the directory structure
// run System("del /S /Q F:\\LMFDB\\data\\*.dat >nul");quit;  to erase the files


/*
Jones();
Bordeaux(1000);
CreatePGLExtensions();

CreateSmallCoefficients(3,15);
CreateSmallCoefficients(4,5);     // many D8
CreateSmallCoefficients(5,5);     // quite slow
CreateSmallCoefficients(6,2);    
CreateSmallCoefficients(7,2);    
CreateSmallCoefficients(8,1);    
CreateSmallCoefficients(9,1);    
CreateSmallCoefficients(10,1);    

CreateWreathProducts();
CreateMalleKluenersFamilies();
CreateFamilies();
CreateGL2Extensions();
CreateSD16Extensions();
*/


//System("del /S /Q F:\\LMFDB\\data\\*.dat >nul");
//_:=Process(Q);

//_:=Process(x^2+15);
//_:=Process(x^4-2);
//_:=Process(x^2-x+1);
//_:=Process(x^2+21);
//_:=Process(x^4-5*x^2+7);


//SetVerbose("ArtRep",2);

/* 
ProcessNFDBdatFile(LMFDBdir*"nfdb_bordeaux.dat"); // Bordeaux database (done)
*/
/* 
ProcessNFDBdatFile(LMFDBdir*"nfdb_jones.dat");    // Jones database (done)
ProcessDirectory(OldArtdir);
*/

/*
ProcessSmallFields(5);     // Done
ProcessSmallFields(20);    // Done
ProcessSmallFields(40);    // Done    

ProcessDirectory(MalleKluenersdir);               // Done
ProcessBosman();                                  // Done
ProcessNFDBdatFile(Familiesdir*"bookerstr.dat");  // Done

ProcessDirectory(Familiesdir);       // Done [deleted some D16]
ProcessPGLExtensions(: easy:=true);  // Done
ProcessSmallCoefficients();          // Done
ProcessNFDBdatFile(SD16file);        // Done
ProcessNFDBdatFile(GL2file);         // Done
ProcessSmallFields(60);              // Done
*/
/* */
//ProcessWreathProducts();             // unfinished   [polredabs slow?]


/* 
Testing old problems:

_:=Process(x^6-9*x^4-21*x^2+21);
_:=Process(x^3-x^2-2*x+1);
_:=Process(x^4-10*x^2+1);
_:=Process(x^3-2);
_:=Process(x^6-x^5-6*x^4+7*x^3+4*x^2-5*x+1); // testing Minimize for C3*S3
_:=Process(x^8-28); // got a "Too many transforms" error once(!)
_:=Process(x^6-x^5-6*x^4+7*x^3+4*x^2-5*x+1);    // Minimize for C3*S3
_:=Process(x^6-9*x^4-21*x^2+21);                // bug in classmap in artin
_:=Process(NumberField(x^6-9*x^4+43*x^2-43));   // hangs
_:=Process(x^6-x^5-10*x^4-20*x^2-128*x+8);      // old PolR coercion in EulerFactor
_:=Process(x^12 - 27*x^9 + 234*x^6 - 702*x^3 + 702); // freezes. High inertia?
_:=Process(x^8-x^7-7*x^6+5*x^5+15*x^4-7*x^3-10*x^2+2*x+1);  // gamma in L-series
_:=Process(x^6+2*x^4+8*x^3+20*x^2+16*x-8);      // LocalData
_:=Process(x^6-18*x^4+36*x^2+24);               // LocalData
*/


/*
C2^3:PSL(2,7)   (n)
// http://mathoverflow.net/questions/95176/octic-family-with-galois-group-of-order-1344
x^8+3*x^7-15*x^6-29*x^5+79*x^4+61*x^3+29*x+16-n*x^2;  


PSL(2,7)=PSL(3,2)   (a,b,t)   Malle 2000

(x^4 + 2*a*x^3 + 2*(a^2 + a* b - 1)*x^2 + 4*b*x + 2*b^2)*
  (x^3 - (a + b + 3)*x^2 + 2*(b + a + 1)*x + b) + t*x^3 * (x - 2)


PGL(2,7) (a,t)  Malle 2000

x*(x^7 + x^6 + 14*a*x^5 + 7a*x^4 + 49*a^2*x^3 + 14*a^2*x^2 + 49*a^3*x + 7*a^3)
+t*(7*x^2 + x + 1)  


2^3.PSL(3,2)*(a,t)  Malle 2000

x^4*(x^2 + a*x + 2*a)*(x - 2)^2 + t*((a - 5)*(x^2 + x) - 2*a - 2)*(x - 1)^2  

  
3^2.GL(2,3)  (a,t)  Malle 2000

(x^2 + (a - 3)*x - a)^3*(x^3 + (3*a - 12)*x^2 + (3*a^2 - 18*a + 36)*x - a^3)
+t*x^2*(x - 3)

Aut(S_6)   (a,t)  Malle 2000

((x + 8)*(x - 1)^4 + 4*(x + 2)*(x - 1)^2*a + 2*x* a^2)^2 + t*(x^2 + 8)


GL(2,9)  a   Malle 2000

2*(8*a - 81)*(2*a - 27)^2*t - 4*a^4 + 360*a^3 - 8856*a^2 + 86022*a - 295245

PSL_2(p) Shih's theorem

*/



///////////////////////////////////////////////////////////////////////
////////////////////////   EXECUTION PART  ////////////////////////////
///////////////////////////////////////////////////////////////////////



//  Add some more fields if not in the database already 

/*
*/
IndexDatabase();


/*
Process(x^8+2);
Process(x^8-2);
Process(x^8-3);
Process(x^7 - 3*x^6 + 3*x^5 - 3*x^4 + 3*x^3 - 3*x^2 + 10*x - 9);
//ProcessNFDBdatFile(LMFDBdir*"tmp1.dat");    // 
ProcessNFDBdatFile(LMFDBdir*"nfdb_bordeaux.dat"); // Bordeaux database
ProcessDirectory(OldArtdir);

*/

//  Export the database to json, write statistics to stats.inc

// ExportToJSON(: test:=false);

/*
Statistics();
*/


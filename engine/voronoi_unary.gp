\\ voronoi_unary.gp -- Voronoi's algorithm for perfect unary Humbert forms over a
\\ totally real cubic field K, in EXACT arithmetic.
\\
\\ A unary form is a totally positive a in K; its value on x in O_K is Tr(a x^2),
\\ a positive definite ternary Z-form in an integral basis.  A perfect form of
\\ minimum 1 has a in K (it solves a rational linear system), so every vertex of
\\ the Ryshkov polyhedron is represented exactly.  Minimal vectors come from
\\ qfminim on an INTEGRAL Gram matrix (exact Fincke-Pohst).  Equivalence tests
\\ (a ~ a u^2, x ~ x u) are exact membership tests in O_K^x and do not use the
\\ (GRH-conditional) units of bnfinit; those are used only to keep heights small.
\\ Output: n_K (perfect forms), E_K (edges), mu_K (extreme-square classes), all
\\ modulo units, together with an Euler-characteristic check n - E + mu = 0.

vu_nf(f) = my(K = nfinit(subst(f, 'x, 'y))); if (K.sign != [3,0], error("need a totally real cubic")); K;

vu_zk(K) = vector(#K.zk, i, Mod(K.zk[i], K.pol));

vu_elt(K, v) = Mod(K.zk * v, K.pol);

vu_gram(K, a) =
{
  my(z = vu_zk(K), n = #z, G = matrix(n, n));
  for (i = 1, n, for (j = i, n, G[i,j] = G[j,i] = trace(a * z[i] * z[j])));
  G;
}

vu_emb(K, a) = my(p = lift(a)); vector(#K.roots, i, subst(p, 'y, K.roots[i]));


\\ Exact enumeration for a positive definite rational 3x3 Gram matrix G:
\\ all x (one per +-pair) with x~Gx <= B.  Uses the exact rational Gauss
\\ reduction; floating point only widens loop ranges (by 1), never decides.
vu_enum(G, B) =
{
  my(Q = qfgaussred(G), res = List(), r3, r2, r1, c2, c1, x);
  r3 = sqrt(B / Q[3,3] * 1.);
  for (x3 = 0, floor(r3) + 1,
    my(R3 = B - Q[3,3] * x3^2); if (R3 < 0, next);
    c2 = -Q[2,3] * x3; r2 = sqrt(R3 / Q[2,2] * 1.);
    for (x2 = floor(c2 - r2) - 1, ceil(c2 + r2) + 1,
      my(R2 = R3 - Q[2,2] * (x2 + Q[2,3] * x3)^2); if (R2 < 0, next);
      c1 = -(Q[1,2] * x2 + Q[1,3] * x3); r1 = sqrt(R2 / Q[1,1] * 1.);
      for (x1 = floor(c1 - r1) - 1, ceil(c1 + r1) + 1,
        if (x3 == 0 && (x2 < 0 || (x2 == 0 && x1 <= 0)), next);
        x = [x1, x2, x3]~;
        if (x~ * G * x <= B, listput(res, x)))));
  Vec(res);
}

\\ LLL-preconditioned exact minimum / strict sub-level set of a rational form
vu_prep(G) = my(L = qflllgram(G)); [L, L~ * G * L];
vu_qmin(G) =
{
  my(P = vu_prep(G), Gr = P[2], m0 = vecmin(vector(3, i, Gr[i,i])), V = vu_enum(Gr, m0), m, W);
  m = vecmin(vector(#V, i, V[i]~ * Gr * V[i]));
  W = select(v -> v~ * Gr * v == m, V);
  [m, P[1] * Mat(W)];
}
vu_qbelow(G) =
{
  my(P = vu_prep(G), Gr = P[2], V = select(v -> v~ * Gr * v < 1, vu_enum(Gr, 1)));
  [#V, if (#V, P[1] * Mat(V), [;])];
}

vu_ispos(a) = polsturm(minpoly(a), [-oo, 0]) == 0;

\\ [minimum, matrix of minimal vectors (one per +-pair)]
vu_minvecs(K, a) = vu_qmin(vu_gram(K, a));

\\ row of the linear functional rho -> Tr(rho x^2) in the integral basis
vu_crow(K, v) = my(z = vu_zk(K), x = vu_elt(K, v)); vector(#z, k, trace(z[k] * x^2));

vu_rows(K, M) = my(k = #M, R = matrix(k, #K.zk)); for (i = 1, k, my(r = vu_crow(K, M[, i])); for (j = 1, #r, R[i, j] = r[j])); R;

\\ a rational number strictly between lo (rational) and hi (real), small height
vu_between(lo, hi) =
{
  my(B = 16, r);
  while (1, r = bestappr((lo + hi) / 2, B); if (r > lo && r < hi, return(r)); B *= 16; if (B > 10^80, error("vu_between")));
}

\\ Walk from a (minimum 1) along rho (Tr(rho x^2) = 0 on the current minimal vectors,
\\ rho not totally nonnegative) to the first form at which a new vector attains 1.
vu_walk(K, a, rho) =
{
  my(ea = vu_emb(K, a), er = vu_emb(K, rho), tmax = oo, t, Ga = vu_gram(K, a), Gr = vu_gram(K, rho),
     G, d, S, v, tn, cnt = 0);
  for (i = 1, #ea, if (er[i] < 0, tmax = min(tmax, -ea[i] / er[i])));
  if (tmax == oo, error("vu_walk: unbounded direction"));
  t = vu_between(0, tmax);
  while (!vu_ispos(a + t * rho), t = t / 2);
  while (1,
    S = vu_qbelow(Ga + t * Gr);
    if (S[1] > 0, break);
    my(t2 = vu_between(t, tmax)); while (!vu_ispos(a + t2 * rho), t2 = (t + t2) / 2); t = t2;
    cnt++; if (cnt > 400, error("vu_walk: no new vector found"));
  );
  while (S[1] > 0,
    tn = t;
    for (j = 1, #S[2], v = S[2][, j]; tn = min(tn, (v~ * Ga * v - 1) / (-(v~ * Gr * v))));
    t = tn;
    S = vu_qbelow(Ga + t * Gr);
  );
  a + t * rho;
}

\\ Edges of the vertex a of the Ryshkov polyhedron: pairs of adjacent vertices of the
\\ polygon conv{q(x) : x in M(a)}; returns [i, j, rho] with rho the edge direction.
vu_edges(K, a, M) =
{
  my(k = #M, R = vu_rows(K, M), E = List(), ker, rho, c, pos, neg, Q, s, ok);
  Q = vector(k, i, my(e = vu_emb(K, vu_elt(K, M[, i]))); vector(#e, j, e[j]^2));
  for (i = 1, k, for (j = i + 1, k,
    ker = matker(matrix(2, #K.zk, r, cc, R[if (r == 1, i, j), cc]));
    if (#ker != 1, next);
    rho = ker[, 1];
    c = vector(k, m, R[m, ] * rho);
    pos = 0; neg = 0;
    for (m = 1, k, if (c[m] > 0, pos++, if (c[m] < 0, neg++)));
    if (pos && neg, next);
    if (pos == 0 && neg == 0, next);
    if (neg, rho = -rho);
    ok = 1;
    for (m = 1, k, if (m != i && m != j && c[m] == 0,
      s = (Q[m] - Q[i]) * (Q[j] - Q[i])~ / ((Q[j] - Q[i]) * (Q[j] - Q[i])~);
      if (s < 0 || s > 1, ok = 0)));
    if (ok, listput(E, [i, j, vu_elt(K, rho)]));
  ));
  Vec(E);
}

vu_isint(K, x) = denominator(nfalgtobasis(K, x)) == 1;

\\ exact test: b = a u^2 for a unit u ?  returns u or 0
vu_formequiv(K, a, b, Einv) =
{
  my(g = b / a, e, s, cc, u);
  if (norm(g) != 1, return(0));
  if (!vu_isint(K, g), return(0));
  e = vu_emb(K, g);
  forvec(sg = vector(#e - 1, i, [0, 1]),
    s = vector(#e, i, if (i == 1, 1, (-1)^sg[i-1]) * sqrt(e[i]));
    cc = round(Einv * s~);
    u = vu_elt(K, cc);
    if (u^2 == g, return(u)));
  0;
}

\\ exact test: x = z u for a unit u
vu_eltequiv(K, x, z) = my(g = x / z); abs(norm(g)) == 1 && vu_isint(K, g);

\\ size reduction modulo units (pw = 2 for forms, 1 for elements); heuristic only.
\\ vu_reduce2 also returns the unit w with  a = a' * w^pw.
vu_reduce2(K, a, U, B2inv, pw) =
{
  my(e = vu_emb(K, a), L = vector(#e, i, log(abs(e[i]))), m = vecsum(L) / #L, c, w);
  c = round(B2inv * [L[1] - m, L[2] - m]~ / pw);
  w = U[1]^c[1] * U[2]^c[2];
  [a * w^(-pw), w];
}
vu_reduce(K, a, U, B2inv, pw) =
{
  my(e = vu_emb(K, a), L = vector(#e, i, log(abs(e[i]))), m = vecsum(L) / #L, c);
  c = round(B2inv * [L[1] - m, L[2] - m]~ / pw);
  a * U[1]^(-pw * c[1]) * U[2]^(-pw * c[2]);
}

\\ Main routine.  Returns
\\ [D, R(bnf), h(bnf), n_K, E_K, mu_K, n-E+mu, index of found units in bnf units,
\\  degree multiset, norms of extreme squares]
vu_voronoi(f, maxv = 10^6) =
{
  my(K = vu_nf(f), bnf = bnfinit(K, 1), U, Bm, B2inv, Emb, Einv, a, mv, M, R, ker, rho,
     forms = List(), fkey = List(), qi = 1, degs = List(), xs = List(), xkey = List(),
     Ed, b, nb, found, ufound = List(), dsum, ends, idx, ex);
  U = apply(u -> if (type(u) == "t_POLMOD", u, nfbasistoalg(K, u)), bnf.fu);
  Bm = matrix(2, 2, i, j, log(abs(vu_emb(K, U[j])[i])) - vecsum(apply(t -> log(abs(t)), vu_emb(K, U[j]))) / 3);
  B2inv = Bm^(-1);
  Emb = matrix(3, 3, i, k, vu_emb(K, Mod(K.zk[k], K.pol))[i]);
  Einv = Emb^(-1);
  \\ ascent from the trace form to a first perfect form
  a = Mod(1, K.pol); mv = vu_minvecs(K, a); a = a / mv[1]; M = vu_minvecs(K, a)[2];
  R = vu_rows(K, M);
  while (matrank(R) < 3,
    ker = matker(R); rho = vu_elt(K, ker[, 1]);
    if (vecmin(vu_emb(K, rho)) >= 0, rho = -rho);
    a = vu_walk(K, a, rho); M = vu_minvecs(K, a)[2]; R = vu_rows(K, M));
  a = vu_reduce(K, a, U, B2inv, 2);
  listput(forms, a); listput(fkey, norm(a));
  while (qi <= #forms,
    a = forms[qi]; qi++;
    mv = vu_minvecs(K, a);
    if (mv[1] != 1, error("minimum not 1"));
    M = mv[2];
    if (matrank(vu_rows(K, M)) < 3, error("form in queue not perfect"));
    Ed = vu_edges(K, a, M);
    listput(degs, #Ed);
    ends = Set(concat(vector(#Ed, e, [Ed[e][1], Ed[e][2]])));
    if (#ends != #Ed, error("polygon: #vertices != #edges"));
    for (e = 1, #ends,
      my(x = vu_elt(K, M[, ends[e]]), nx = abs(norm(x)), isnew = 1);
      for (j = 1, #xs, if (xkey[j] == nx && vu_eltequiv(K, x, xs[j]), isnew = 0; break));
      if (isnew, listput(xs, vu_reduce(K, x, U, B2inv, 1)); listput(xkey, nx)));
    for (e = 1, #Ed,
      my(rw = vu_reduce2(K, vu_walk(K, a, Ed[e][3]), U, B2inv, 2));
      b = rw[1]; nb = norm(b); found = 0;
      for (j = 1, #forms, if (fkey[j] == nb, my(u = vu_formequiv(K, forms[j], b, Einv)); if (u, found = j; listput(ufound, u * rw[2]); break)));
      if (!found, listput(forms, b); listput(fkey, nb); listput(ufound, rw[2]));
      if (#forms > maxv, error("too many forms")));
  );
  dsum = vecsum(Vec(degs));
  \\ index of the subgroup generated by the side-pairing units inside bnf's unit group
  ex = vector(#ufound, i, my(w = bnfisunit(bnf, ufound[i])); if (#w == 0, error("unit not in bnf group")); [w[1], w[2]]);
  ex = matrix(2, #ufound, i, j, ex[j][i]);
  idx = if (#ufound == 0, 0, my(H = mathnf(ex)); if (#H == 2, abs(matdet(H)), -#H));
  [K.disc, bnf.reg, bnf.no, #forms, dsum / 2, #xs, #forms - dsum / 2 + #xs, idx, vecsort(Vec(degs)), vecsort(Vec(xkey))];
}

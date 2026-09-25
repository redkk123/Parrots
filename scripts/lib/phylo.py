"""Minimal phylogenetic utilities (no R available in this environment).

Covers what Phase 0 needs: reading an ape `phylo` object, Newick I/O, pruning,
covariance matrices under Brownian motion (BM), Pagel's lambda and OU transforms,
GLS/PGLS by maximum likelihood, and BM ancestral-state estimates.
"""
from __future__ import annotations

import math
import re

import numpy as np
from scipy import optimize


class Tree:
    """Rooted tree stored as parent/children arrays. Node 0 is the root."""

    def __init__(self, parent, length, label):
        self.parent = list(parent)          # parent[i] = index of parent, -1 for root
        self.length = np.asarray(length, float)
        self.label = list(label)            # tip label or "" for internal nodes
        self.children = [[] for _ in parent]
        for i, p in enumerate(parent):
            if p >= 0:
                self.children[p].append(i)
        self.root = self.parent.index(-1)

    # ---- construction -------------------------------------------------
    @classmethod
    def from_ape(cls, edge, edge_length, tip_label):
        edge = np.asarray(edge, int)
        ntip = len(tip_label)
        nnode = int(edge.max())
        parent = [-1] * nnode
        length = [0.0] * nnode
        for (a, b), l in zip(edge, edge_length):
            parent[b - 1] = a - 1
            length[b - 1] = float(l)
        label = [str(t) for t in tip_label] + [""] * (nnode - ntip)
        return cls(parent, length, label)

    @classmethod
    def from_newick(cls, text):
        tokens = re.findall(r"\(|\)|,|;|:[^,();]+|[^,();:]+", text.strip())
        parent, length, label, stack, cur = [], [], [], [], None
        for tok in tokens:
            if tok == "(":
                idx = len(parent)
                parent.append(stack[-1] if stack else -1)
                length.append(0.0)
                label.append("")
                stack.append(idx)
                cur = None
            elif tok == ",":
                cur = None
            elif tok == ")":
                cur = stack.pop()
            elif tok == ";":
                break
            elif tok.startswith(":"):
                length[cur] = float(tok[1:])
            else:
                if cur is not None:  # internal-node label after ')': ignored
                    continue
                idx = len(parent)
                parent.append(stack[-1])
                length.append(0.0)
                label.append(tok.strip())
                cur = idx
        return cls(parent, length, label)

    # ---- queries -------------------------------------------------------
    def is_tip(self, i):
        return not self.children[i]

    def tips(self):
        return [i for i in range(len(self.parent)) if self.is_tip(i)]

    def tip_labels(self):
        return [self.label[i] for i in self.tips()]

    def postorder(self):
        out, stack = [], [(self.root, False)]
        while stack:
            n, done = stack.pop()
            if done:
                out.append(n)
            else:
                stack.append((n, True))
                stack.extend((c, False) for c in reversed(self.children[n]))
        return out

    def depths(self):
        d = np.zeros(len(self.parent))
        for n in reversed(self.postorder()):  # preorder
            if self.parent[n] >= 0:
                d[n] = d[self.parent[n]] + self.length[n]
        return d

    def descendants_tips(self):
        desc = [[] for _ in self.parent]
        for n in self.postorder():
            desc[n] = [n] if self.is_tip(n) else [t for c in self.children[n] for t in desc[c]]
        return desc

    def is_ultrametric(self, tol=1e-4):
        d = self.depths()[self.tips()]
        return (d.max() - d.min()) / d.max() < tol

    # ---- transformations ------------------------------------------------
    def prune(self, keep):
        """Return a new tree with only tips whose label is in `keep`; unary nodes are collapsed."""
        keep = set(keep)
        new_parent, new_len, new_label = [], [], []

        def rec(n, acc_len):
            if self.is_tip(n):
                if self.label[n] not in keep:
                    return None
                return ("tip", self.label[n], acc_len + self.length[n])
            kids = [k for k in (rec(c, 0.0) for c in self.children[n]) if k is not None]
            if not kids:
                return None
            if len(kids) == 1:
                kind, payload, l = kids[0]
                return (kind, payload, l + self.length[n] + acc_len)
            return ("node", kids, acc_len + self.length[n])

        def build(item, parent_idx):
            kind, payload, l = item
            idx = len(new_parent)
            new_parent.append(parent_idx)
            new_len.append(l if parent_idx >= 0 else 0.0)
            new_label.append(payload if kind == "tip" else "")
            if kind == "node":
                for k in payload:
                    build(k, idx)

        root = rec(self.root, 0.0)
        build(root, -1)
        return Tree(new_parent, new_len, new_label)

    def to_newick(self, digits=6):
        def rec(n):
            s = self.label[n] if self.is_tip(n) else "(" + ",".join(rec(c) for c in self.children[n]) + ")"
            if self.parent[n] >= 0:
                s += f":{self.length[n]:.{digits}f}"
            return s
        return rec(self.root) + ";"

    # ---- covariance ------------------------------------------------------
    def vcv_all(self):
        """BM covariance among all nodes (tips + internal), ordered by node index."""
        n = len(self.parent)
        sub = [[] for _ in range(n)]
        for v in self.postorder():
            sub[v] = [v] + [x for c in self.children[v] for x in sub[c]]
        C = np.zeros((n, n))
        for v in range(n):
            if self.parent[v] >= 0 and self.length[v] > 0:
                idx = np.array(sub[v])
                C[np.ix_(idx, idx)] += self.length[v]
        return C

    def vcv(self, labels):
        C = self.vcv_all()
        idx = [self.label.index(l) for l in labels]
        return C[np.ix_(idx, idx)]


# ---- models ---------------------------------------------------------------

def lambda_transform(C, lam):
    D = np.diag(np.diag(C))
    return lam * (C - D) + D


def ou_transform(C, alpha):
    """OU covariance with the root state treated as fixed (ultrametric tree, height T)."""
    if alpha < 1e-8:
        return C
    dij = np.add.outer(np.diag(C), np.diag(C)) - 2 * C  # patristic distance
    return np.exp(-alpha * dij) * (1 - np.exp(-2 * alpha * C)) / (2 * alpha)


def gls_loglik(y, X, V, reml=False):
    """Profile log-likelihood of y ~ X under covariance sigma2 * V."""
    n, k = X.shape
    L = np.linalg.cholesky(V)
    Xs = np.linalg.solve(L, X)
    ys = np.linalg.solve(L, y)
    beta, *_ = np.linalg.lstsq(Xs, ys, rcond=None)
    r = ys - Xs @ beta
    df = n - k if reml else n
    s2 = float(r @ r) / df
    logdet = 2 * np.sum(np.log(np.diag(L)))
    ll = -0.5 * (df * math.log(2 * math.pi * s2) + logdet + df)
    if reml:
        ll -= 0.5 * np.linalg.slogdet(Xs.T @ Xs)[1]
    cov_beta = s2 * np.linalg.inv(Xs.T @ Xs)
    return ll, beta, cov_beta, s2


def fit_pgls(y, X, C, model="lambda"):
    """Fit a GLS with BM, Pagel's lambda or OU covariance by maximum likelihood."""
    Cn = C / np.max(np.diag(C))  # unit-height tree; sigma2 absorbs the scale

    def V_of(p):
        if model == "BM":
            return Cn
        if model == "lambda":
            return lambda_transform(Cn, p)
        if model == "OU":
            return ou_transform(Cn, p)
        raise ValueError(model)

    if model == "BM":
        p_hat = None
    elif model == "lambda":
        res = optimize.minimize_scalar(lambda p: -gls_loglik(y, X, V_of(p))[0], bounds=(0.0, 1.0),
                                       method="bounded", options={"xatol": 1e-4})
        p_hat = float(res.x)
    else:  # OU alpha on a log scale, in units of 1/tree-height
        res = optimize.minimize_scalar(lambda la: -gls_loglik(y, X, V_of(math.exp(la)))[0],
                                       bounds=(-6, 5), method="bounded")
        p_hat = float(math.exp(res.x))
    V = V_of(p_hat)
    ll, beta, cov_beta, s2 = gls_loglik(y, X, V)
    k = X.shape[1] + 1 + (0 if model == "BM" else 1)
    return {"model": model, "param": p_hat, "loglik": ll, "aic": 2 * k - 2 * ll, "beta": beta,
            "se": np.sqrt(np.diag(cov_beta)), "sigma2": s2, "V": V, "k": k}


def bm_ancestral(tree, tip_values):
    """ML (= conditional expectation) ancestral states under BM, with conditional SEs.

    tip_values: dict label -> value, covering every tip of `tree`.
    Returns arrays indexed by node: estimate, standard error.
    """
    C = tree.vcv_all()
    tips = tree.tips()
    y = np.array([tip_values[tree.label[t]] for t in tips])
    Ctt = C[np.ix_(tips, tips)]
    one = np.ones(len(tips))
    Ci = np.linalg.inv(Ctt)
    mu = float(one @ Ci @ y) / float(one @ Ci @ one)
    s2 = float((y - mu) @ Ci @ (y - mu)) / (len(y) - 1)
    Cat = C[:, tips]
    est = mu + Cat @ Ci @ (y - mu)
    # conditional variance, plus uncertainty of mu
    cond = np.diag(C) - np.einsum("ij,jk,ik->i", Cat, Ci, Cat)
    w = 1 - Cat @ Ci @ one
    var = s2 * (cond + w ** 2 / float(one @ Ci @ one))
    return est, np.sqrt(np.clip(var, 0, None)), mu, s2

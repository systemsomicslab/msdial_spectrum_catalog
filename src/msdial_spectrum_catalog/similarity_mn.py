"""Versioned, dependency-free molecular-networking similarity kernels.

Low-level kernels accept already processed spectra. Use ``score`` for raw
spectra and a fully specified, canonical preprocessing identity.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from decimal import Decimal
import math
import re

from .similarity import Peaks, entropy_similarity

KERNELS_VERSION = "1"


@dataclass(frozen=True)
class Spectrum:
    precursor_mz: float
    peaks: tuple[tuple[float, float], ...]

    def __post_init__(self):
        # Also make list-based caller input immutable, without sorting it.
        object.__setattr__(self, "peaks", tuple(tuple(p) for p in self.peaks))


@dataclass(frozen=True)
class MnPrep:
    rel: float = 1.0
    abs: float = 0.0
    above: float | None = 0.05
    minmz: float = 0.0
    maxmz: float = 10000.0
    frame: float = 1.0
    pow: float = 0.5
    scale: float = 100.0
    precwin: float = 0.0
    topk: int = 0
    topkwin: float = 50.0


MSDIAL_MN_V1 = MnPrep()


@dataclass(frozen=True)
class KernelResult:
    kernel_id: str
    score: float | None
    score_convention: str
    matched_peak_count: int
    details: dict


@dataclass(frozen=True)
class MsdialMatch:
    e: str
    r: str
    pairs: tuple[tuple[int, int, str], ...]
    used_e: tuple[bool, ...]
    used_r: tuple[bool, ...]


def _sum(values):
    # Python 3.12+ sum(float) uses compensated accumulation; LINQ Sum and
    # matchms score_best_matches add sequentially. Pin that order on 3.10+.
    total = 0.0
    for value in values:
        total += value
    return total

def _divide(numerator, denominator):
    """C# floating point division, including its zero/zero behavior."""
    if denominator == 0:
        return math.nan if numerator == 0 else math.copysign(math.inf, numerator)
    return numerator / denominator


def preprocess(spectrum: Spectrum, prep: MnPrep = MSDIAL_MN_V1) -> Spectrum:
    """Refine, window/top-k, frame, normalise; reject nonfinite input intensities."""
    if any(not math.isfinite(intensity) for _, intensity in spectrum.peaks):
        raise ValueError("preprocess requires finite input intensities")
    base = max((p[1] for p in spectrum.peaks), default=0.0)
    peaks = [(mz, intensity) for mz, intensity in spectrum.peaks
             if not (mz < prep.minmz or mz > prep.maxmz
                     or (prep.above is not None and mz > spectrum.precursor_mz + prep.above)
                     or intensity < prep.abs)
             and intensity >= base * prep.rel * 0.01]
    if prep.precwin > 0:
        peaks = [p for p in peaks if abs(p[0] - spectrum.precursor_mz) > prep.precwin]
    if prep.topk > 0:
        peaks = [p for i, p in enumerate(peaks)
                 if sum(j != i and abs(q[0] - p[0]) <= prep.topkwin and q[1] > p[1]
                        for j, q in enumerate(peaks)) < prep.topk]
    if prep.frame > 0:
        if any(peaks[i][0] < peaks[i - 1][0] for i in range(1, len(peaks))):
            peaks.sort(key=lambda p: p[0])
        framed = []
        last_frame = None
        apex_intensity = 0.0
        for mz, intensity in peaks:
            frame = int(mz / prep.frame)
            if framed and frame == last_frame:
                apex, total = framed[-1]
                total += intensity
                if intensity > apex_intensity:
                    apex_intensity, apex = intensity, mz
                framed[-1] = (apex, total)
            else:
                framed.append((mz, intensity))
                apex_intensity = intensity
                last_frame = frame
        peaks = framed
    if peaks:
        maximum = max(p[1] for p in peaks) ** prep.pow
        peaks = [(mz, _divide(intensity ** prep.pow, maximum) * prep.scale)
                 for mz, intensity in peaks]
    return Spectrum(spectrum.precursor_mz, tuple(peaks))


def msdial_match(a: Spectrum, b: Spectrum, tol: float = 0.025) -> MsdialMatch:
    """SearchMatchedPeaks, strict comparisons and original loop order, no mutation."""
    if a.precursor_mz < b.precursor_mz:
        e, r, ename, rname = a, b, "a", "b"
    else:
        e, r, ename, rname = b, a, "b", "a"
    used_e, used_r = [False] * len(e.peaks), [False] * len(r.peaks)
    pairs = []
    dprec = r.precursor_mz - e.precursor_mz
    for i, (r_mz, r_int) in enumerate(r.peaks):
        best, difference, kind = -1, float.fromhex('0x1.fffffffffffffp+1023'), "product"
        for j, (e_mz, e_int) in enumerate(e.peaks):
            if used_e[j]:
                continue
            if abs(e_mz - r_mz) < tol:
                candidate = "product"
            elif abs((dprec + e_mz) - r_mz) < tol:
                candidate = "neutral_loss"
            else:
                continue
            delta = abs(e_int - r_int)
            if delta < difference:
                best, difference, kind = j, delta, candidate
        if best >= 0:
            used_r[i] = used_e[best] = True
            pairs.append((i, best, kind))
    return MsdialMatch(ename, rname, tuple(pairs), tuple(used_e), tuple(used_r))


def _matched_values(a, b, match, product_only=False):
    e, r = (a, b) if match.e == "a" else (b, a)
    return [(r.peaks[i][1], e.peaks[j][1]) for i, j, kind in match.pairs
            if not product_only or kind == "product"]


def _matched_cosine(values, zero_guard=False):
    product = _sum(r * e for r, e in values)
    sr = _sum(r * r for r, e in values)
    se = _sum(e * e for r, e in values)
    if zero_guard and (sr == 0 or se == 0):
        return 0.0
    return _divide(product, math.sqrt(sr) * math.sqrt(se))


def _bonanza(a, b, match):
    if not match.pairs:
        return 0.0
    product = _sum(r * e for r, e in _matched_values(a, b, match))
    flags_a, flags_b = ((match.used_e, match.used_r) if match.e == "a"
                        else (match.used_r, match.used_e))
    # Keep prop1 then prop2, including the addition order of the C# expression.
    ua = _sum(intensity ** 2 for (_, intensity), used in zip(a.peaks, flags_a) if not used)
    ub = _sum(intensity ** 2 for (_, intensity), used in zip(b.peaks, flags_b) if not used)
    return _divide(product, (product + ua) + ub)


def _msdial_result(base, a, b, tol):
    match = msdial_match(a, b, tol)
    if base == "msdial.bonanza.v1":
        value, convention = _bonanza(a, b, match), "bonanza_ratio"
    else:
        value = (_matched_cosine(_matched_values(a, b, match,
                 product_only=base == "msdial.cosine_matched.legacy"),
                 zero_guard=base == "msdial.cosine_matched.legacy") if match.pairs else 0.0)
        convention = "cosine_matched_only"
    details = {"e": match.e, "r": match.r, "pairs": [list(p) for p in match.pairs],
               "n_product": sum(p[2] == "product" for p in match.pairs),
               "n_neutral_loss": sum(p[2] == "neutral_loss" for p in match.pairs)}
    return KernelResult(format_kernel_id(base, {"prep": "none", "tol": tol}),
                        value, convention, len(match.pairs), details)


def msdial_moddot_legacy(a: Spectrum, b: Spectrum, tol: float = 0.025) -> KernelResult:
    return _msdial_result("msdial.moddot.legacy", a, b, tol)


def msdial_bonanza_v1(a: Spectrum, b: Spectrum, tol: float = 0.025) -> KernelResult:
    return _msdial_result("msdial.bonanza.v1", a, b, tol)


def msdial_all_v1(a: Spectrum, b: Spectrum, tol: float = 0.025) -> tuple[float, int, float, float]:
    """C# All result: (bonanza, n, zero-guarded moddot, product-only cosine)."""
    match = msdial_match(a, b, tol)
    if not match.pairs:
        return 0.0, 0, 0.0, 0.0
    return (_bonanza(a, b, match), len(match.pairs),
            _matched_cosine(_matched_values(a, b, match), True),
            _matched_cosine(_matched_values(a, b, match, True), True))


def _hungarian(weights):
    """Maximum-weight assignment of a square matrix, O(n**3), potentials method."""
    n = len(weights)
    u, v, p, way = [0.0] * (n + 1), [0.0] * (n + 1), [0] * (n + 1), [0] * (n + 1)
    for i in range(1, n + 1):
        p[0] = i
        j0 = 0
        minimum, used = [math.inf] * (n + 1), [False] * (n + 1)
        while True:
            used[j0] = True
            i0, delta, j1 = p[j0], math.inf, 0
            for j in range(1, n + 1):
                if not used[j]:
                    current = -weights[i0 - 1][j - 1] - u[i0] - v[j]
                    if current < minimum[j]:
                        minimum[j], way[j] = current, j0
                    if minimum[j] < delta:
                        delta, j1 = minimum[j], j
            for j in range(n + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minimum[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while True:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            if j0 == 0:
                break
    return [(p[j] - 1, j - 1) for j in range(1, n + 1)]


def _optimal(candidates):
    edges = {}
    adjacency = {}
    for i, j, w in candidates:
        edges[i, j] = max(w, edges.get((i, j), -math.inf))
        adjacency.setdefault((0, i), set()).add((1, j))
        adjacency.setdefault((1, j), set()).add((0, i))
    visited, accepted = set(), []
    for start in sorted(adjacency):
        if start in visited:
            continue
        stack, left, right = [start], [], []
        visited.add(start)
        while stack:
            side, index = stack.pop()
            (left if side == 0 else right).append(index)
            for neighbour in sorted(adjacency[side, index]):
                if neighbour not in visited:
                    visited.add(neighbour)
                    stack.append(neighbour)
        left.sort()
        right.sort()
        n = max(len(left), len(right))
        weights = [[0.0] * n for _ in range(n)]
        for x, i in enumerate(left):
            for y, j in enumerate(right):
                weights[x][y] = edges.get((i, j), 0.0)
        for x, y in _hungarian(weights):
            if x < len(left) and y < len(right) and (left[x], right[y]) in edges:
                accepted.append((left[x], right[y], edges[left[x], right[y]]))
    return accepted


def modified_cosine(a: Spectrum, b: Spectrum, *, tol: float, assign: str,
                    mz_power: float = 0.0, intensity_power: float = 1.0) -> KernelResult:
    """matchms 0.34.0 windows/order/norms, greedy or component-wise optimal."""
    identity = format_kernel_id("gnps.modified_cosine", {
        "prep": "none", "tol": tol, "assign": assign, "mzpow": mz_power, "ipow": intensity_power})
    ap = sorted(a.peaks, key=lambda p: p[0])
    bp = sorted(b.peaks, key=lambda p: p[0])
    wa = [mz ** mz_power * intensity ** intensity_power for mz, intensity in ap]
    wb = [mz ** mz_power * intensity ** intensity_power for mz, intensity in bp]
    # Degenerate MN normalization can produce NaN. Do not feed nonfinite
    # costs to Hungarian (no finite augmenting path exists in that case).
    if any(not math.isfinite(w) for w in wa + wb):
        return KernelResult(identity, math.nan, "cosine", 0, {"pairs": [], "nonfinite_weights": True})
    candidates = []
    shift = a.precursor_mz - b.precursor_mz
    for offset in ([0.0] if abs(shift) <= tol else [0.0, shift]):
        lowest = 0
        for i, (mz, _) in enumerate(ap):
            for j in range(lowest, len(bp)):
                shifted_mz = bp[j][0] + offset
                if shifted_mz < mz - tol:
                    lowest = j + 1
                elif shifted_mz > mz + tol:
                    break
                else:
                    candidates.append((i, j, wa[i] * wb[j]))
    if assign == "greedy":
        used_a, used_b, accepted = set(), set(), []
        for i, j, w in reversed(sorted(candidates, key=lambda p: p[2])):
            if i not in used_a and j not in used_b:
                used_a.add(i)
                used_b.add(j)
                accepted.append((i, j, w))
    else:
        accepted = _optimal(candidates)
    denominator = math.sqrt(_sum(w * w for w in wa)) * math.sqrt(_sum(w * w for w in wb))
    value = _sum(w for _, _, w in accepted) / denominator if denominator else 0.0
    return KernelResult(identity, value, "cosine", len(accepted),
                        {"pairs": [[i, j] for i, j, _ in accepted], "shift": shift})


def entropy_li2021(a: Spectrum, b: Spectrum, bin: float = 0.05) -> KernelResult:
    identity = format_kernel_id("entropy.li2021", {"prep": "none", "bin": bin, "weighting": "none"})
    pa: Peaks = list(a.peaks)
    pb: Peaks = list(b.peaks)
    comparable = all(all(math.isfinite(intensity) for _, intensity in peaks)
                     and _sum(intensity for _, intensity in peaks) > 0 for peaks in (pa, pb))
    value = entropy_similarity(pa, pb, bin) if comparable else None
    return KernelResult(identity, value, "entropy_similarity", 0, {})


_BASE_KEYS = {
    "msdial.moddot.legacy": {"tol"},
    "msdial.bonanza.v1": {"tol"},
    "msdial.cosine_matched.legacy": {"tol"},
    "gnps.modified_cosine": {"tol", "assign", "mzpow", "ipow"},
    "entropy.li2021": {"bin", "weighting"},
}
_PREP_DEFAULTS = {"prep." + k: v for k, v in asdict(MSDIAL_MN_V1).items()}
_DECIMAL = re.compile(r"-?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)\Z")


def _validate(base, params):
    if base not in _BASE_KEYS:
        raise ValueError("unknown kernel base")
    if params.get("prep") not in ("none", "mn"):
        raise ValueError("prep must be none or mn")
    required = _BASE_KEYS[base] | {"prep"}
    if params["prep"] == "mn":
        required |= _PREP_DEFAULTS.keys()
    if set(params) != required:
        raise ValueError("unknown or missing kernel keys")
    for key, value in params.items():
        if key == "prep":
            continue
        if key in ("assign", "weighting"):
            if value not in (("greedy", "optimal") if key == "assign" else ("none",)):
                raise ValueError("invalid " + key)
        elif key == "prep.above" and value is None:
            continue
        else:
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError("numeric value required: " + key)
            try:
                finite = math.isfinite(value)
            except OverflowError:
                finite = False
            if not finite or value < 0 or (key in ("tol", "bin") and value == 0):
                raise ValueError("invalid numeric value: " + key)
            if key == "prep.topk" and value != int(value):
                raise ValueError("topk must be an integer")
    if params["prep"] == "mn" and params["prep.minmz"] > params["prep.maxmz"]:
        raise ValueError("minmz exceeds maxmz")


def _literal(value):
    if value is None:
        return "none"
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, (int, float)):
        if value == 0:
            return "0"
        result = format(Decimal(str(value)), "f")
        return result.rstrip("0").rstrip(".") if "." in result else result
    return value


def format_kernel_id(base: str, params: dict) -> str:
    _validate(base, params)
    return base + "?" + "&".join(k + "=" + _literal(params[k]) for k in sorted(params))


def parse_kernel_id(identity: str) -> tuple[str, dict]:
    if not isinstance(identity, str):
        raise ValueError("identity must be a string")
    if "?" not in identity:
        if identity in ("msdial.moddot.legacy", "msdial.bonanza.v1"):
            return identity, {"prep": "mn", **_PREP_DEFAULTS, "tol": 0.025}
        if identity == "entropy.li2021":
            return identity, {"prep": "none", "bin": 0.05, "weighting": "none"}
        raise ValueError("unknown or non-aliased bare kernel")
    base, query = identity.split("?", 1)
    params = {}
    for part in query.split("&"):
        if part.count("=") != 1:
            raise ValueError("invalid key=value")
        key, value = part.split("=")
        if key in params:
            raise ValueError("duplicate key")
        if key in ("prep", "assign", "weighting"):
            params[key] = value
        elif value == "none":
            params[key] = None
        elif _DECIMAL.fullmatch(value):
            params[key] = float(value)
        else:
            raise ValueError("invalid numeric spelling")
    _validate(base, params)
    if "prep.topk" in params:
        params["prep.topk"] = int(params["prep.topk"])
    return base, params


def canonical_kernel_id(identity: str) -> str:
    return format_kernel_id(*parse_kernel_id(identity))


def score(kernel_id: str, a_raw: Spectrum, b_raw: Spectrum) -> KernelResult:
    base, params = parse_kernel_id(kernel_id)
    identity = format_kernel_id(base, params)
    if params["prep"] == "mn":
        prep = MnPrep(**{k[5:]: v for k, v in params.items() if k.startswith("prep.")})
        a, b = preprocess(a_raw, prep), preprocess(b_raw, prep)
    else:
        a = Spectrum(a_raw.precursor_mz, tuple(sorted(a_raw.peaks, key=lambda p: p[0])))
        b = Spectrum(b_raw.precursor_mz, tuple(sorted(b_raw.peaks, key=lambda p: p[0])))
    if base.startswith("msdial."):
        result = _msdial_result(base, a, b, params["tol"])
    elif base == "gnps.modified_cosine":
        result = modified_cosine(a, b, tol=params["tol"], assign=params["assign"],
                                 mz_power=params["mzpow"], intensity_power=params["ipow"])
    else:
        result = entropy_li2021(a, b, params["bin"])
    return KernelResult(identity, result.score, result.score_convention,
                        result.matched_peak_count, result.details)


def msdial_edge_passes(result: KernelResult, dprec: float, *, max_dprec: float = 400,
                       min_match: int = 2, cutoff_pct: float = 50) -> bool:
    if abs(dprec) > max_dprec:
        return False
    if result.matched_peak_count < min_match:
        return False
    if result.score is None or result.score < cutoff_pct * 0.01:
        return False
    return True

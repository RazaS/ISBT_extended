#!/usr/bin/env python3
import base64
import csv
import json
import os
import re
from html import unescape
from datetime import datetime
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from io import BytesIO
from pathlib import Path
from urllib.error import URLError
from urllib.parse import quote, urlparse
from urllib.request import Request, urlopen

MPL_CACHE_DIR = os.environ.get("MPLCONFIGDIR", "/tmp/matplotlib-cache")
os.environ.setdefault("MPLCONFIGDIR", MPL_CACHE_DIR)
try:
    os.makedirs(MPL_CACHE_DIR, exist_ok=True)
except Exception:
    pass

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from dna_features_viewer import GraphicFeature, GraphicRecord

    VIZ_AVAILABLE = True
    VIZ_IMPORT_ERROR = ""
except Exception as exc:
    plt = None
    GraphicFeature = None
    GraphicRecord = None
    VIZ_AVAILABLE = False
    VIZ_IMPORT_ERROR = clean_text(exc) if "clean_text" in globals() else str(exc)


def clean_text(value):
    if value is None:
        return ""
    text = str(value)
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    text = text.replace("_x000D_", "")
    return text


def read_csv_rows(path):
    rows = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            cleaned = {k: clean_text(v) for k, v in row.items()}
            rows.append(cleaned)
    return rows


def to_int(value):
    txt = clean_text(value).strip().replace(",", "")
    if not txt:
        return None
    try:
        return int(float(txt))
    except ValueError:
        return None


REF_FASTA_PATH = clean_text(os.environ.get("REF_FASTA_PATH", "")).strip()
VIZ_SNV_CONTEXT_BP = 14
_FASTA_READER_CACHE = None
ENSEMBL_SEQUENCE_URL = "https://rest.ensembl.org/sequence/region/human/"
SNV_SEQUENCE_CACHE_PATH = Path(__file__).resolve().parent.parent / "variant_sequence_cache.json"
_SNV_SEQUENCE_CACHE = None


def load_snv_sequence_cache():
    global _SNV_SEQUENCE_CACHE
    if _SNV_SEQUENCE_CACHE is not None:
        return _SNV_SEQUENCE_CACHE
    try:
        if SNV_SEQUENCE_CACHE_PATH.exists():
            with SNV_SEQUENCE_CACHE_PATH.open("r", encoding="utf-8") as handle:
                data = json.load(handle)
                if isinstance(data, dict):
                    _SNV_SEQUENCE_CACHE = data
                    return _SNV_SEQUENCE_CACHE
    except Exception:
        pass
    _SNV_SEQUENCE_CACHE = {}
    return _SNV_SEQUENCE_CACHE


def save_snv_sequence_cache():
    cache = load_snv_sequence_cache()
    try:
        with SNV_SEQUENCE_CACHE_PATH.open("w", encoding="utf-8") as handle:
            json.dump(cache, handle, ensure_ascii=False)
    except Exception:
        return False
    return True


def short_label(value, max_len=42):
    txt = clean_text(value).strip()
    if len(txt) <= max_len:
        return txt
    return txt[: max_len - 3] + "..."


def parse_positions_from_text(value):
    txt = clean_text(value)
    pair = re.search(r"(\d+)\D+(\d+)", txt)
    if pair:
        a = int(pair.group(1))
        b = int(pair.group(2))
        return (min(a, b), max(a, b))
    single = re.search(r"c\.(-?\d+)", txt)
    if single:
        center = int(single.group(1))
        return (center, center + 1)
    return None


def load_fasta_reader():
    global _FASTA_READER_CACHE
    if not REF_FASTA_PATH:
        return None, "Reference FASTA not configured."
    fasta_path = Path(REF_FASTA_PATH).expanduser()
    if not fasta_path.exists():
        return None, f"Reference FASTA not found: {fasta_path}"
    fai_path = Path(str(fasta_path) + ".fai")
    if not fai_path.exists():
        return None, f"FASTA index (.fai) not found: {fai_path}"

    cached = _FASTA_READER_CACHE
    if cached and cached.get("path") == str(fasta_path):
        return cached, ""

    index = {}
    try:
        with fai_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                chrom = parts[0]
                index[chrom] = {
                    "length": int(parts[1]),
                    "offset": int(parts[2]),
                    "line_bases": int(parts[3]),
                    "line_width": int(parts[4]),
                }
    except Exception as exc:
        return None, "Failed reading FASTA index: " + clean_text(exc)

    _FASTA_READER_CACHE = {"path": str(fasta_path), "index": index}
    return _FASTA_READER_CACHE, ""


def resolve_chromosome(chrom_value, index):
    raw = clean_text(chrom_value).strip()
    if not raw:
        return ""
    base = raw[3:] if raw.lower().startswith("chr") else raw
    candidates = [raw, base, "chr" + base, raw.upper(), base.upper(), "chr" + base.upper()]
    lower_map = {k.lower(): k for k in index.keys()}
    for cand in candidates:
        if cand in index:
            return cand
        mapped = lower_map.get(cand.lower())
        if mapped:
            return mapped
    return ""


def fetch_reference_window_local(chrom_value, start_1based, end_1based):
    reader, err = load_fasta_reader()
    if err:
        return "", "", err, ""

    index = reader.get("index", {})
    chrom = resolve_chromosome(chrom_value, index)
    if not chrom:
        return "", "", "Chromosome not found in FASTA index: " + clean_text(chrom_value), ""

    meta = index.get(chrom, {})
    chrom_len = int(meta.get("length", 0))
    start = max(1, int(start_1based))
    end = min(chrom_len, int(end_1based))
    if end < start:
        return "", "", "Invalid coordinate window.", ""

    try:
        seq_chunks = []
        pos = start
        with Path(reader["path"]).open("rb") as handle:
            while pos <= end:
                line_bases = int(meta["line_bases"])
                line_width = int(meta["line_width"])
                base_offset = int(meta["offset"])
                line_idx = (pos - 1) // line_bases
                pos_in_line = (pos - 1) % line_bases
                take = min(line_bases - pos_in_line, end - pos + 1)
                byte_offset = base_offset + line_idx * line_width + pos_in_line
                handle.seek(byte_offset)
                raw = handle.read(take).decode("ascii", errors="ignore")
                seq_chunks.append(re.sub(r"[^A-Za-z]", "", raw))
                pos += take
        return "".join(seq_chunks).upper(), chrom, "", "local_fasta"
    except Exception as exc:
        return "", "", "Failed reading FASTA: " + clean_text(exc), ""


def normalize_chrom_for_ensembl(chrom_value):
    raw = clean_text(chrom_value).strip()
    if raw.lower().startswith("chr"):
        raw = raw[3:]
    return raw


def fetch_reference_window_remote(chrom_value, start_1based, end_1based):
    chrom = normalize_chrom_for_ensembl(chrom_value)
    if not chrom:
        return "", "", "Missing chromosome for remote sequence lookup.", ""
    start = max(1, int(start_1based))
    end = max(start, int(end_1based))
    region = f"{chrom}:{start}..{end}:1"
    url = ENSEMBL_SEQUENCE_URL + quote(region)
    try:
        request = Request(
            url,
            headers={
                "Accept": "application/json",
                "Content-Type": "application/json",
                "User-Agent": "geno-mock-python/1.0 (+https://github.com/RazaS/geno_mock_python)",
            },
        )
        with urlopen(request, timeout=20) as response:
            payload = json.loads(response.read().decode("utf-8"))
        seq = clean_text(payload.get("seq", "")).upper()
        if not seq:
            return "", "", "Remote sequence API returned empty sequence.", ""
        return seq, chrom, "", "ensembl_rest"
    except Exception as exc:
        return "", "", "Remote sequence lookup failed: " + clean_text(exc), ""


def fetch_reference_window(chrom_value, start_1based, end_1based):
    errors = []
    if REF_FASTA_PATH:
        seq, chrom, err, source = fetch_reference_window_local(chrom_value, start_1based, end_1based)
        if not err and seq:
            return seq, chrom, "", source
        if err:
            errors.append(err)

    seq, chrom, err, source = fetch_reference_window_remote(chrom_value, start_1based, end_1based)
    if not err and seq:
        return seq, chrom, "", source
    if err:
        errors.append(err)
    return "", "", " | ".join(errors) if errors else "Sequence lookup failed.", ""


def get_ref_alt_bases(row):
    ref = clean_text(row.get("Ref_allele_curated") or row.get("Ref_allele") or row.get("Ref") or "").strip().upper()
    alt = clean_text(row.get("Alt_allele_curated") or row.get("Alt_allele") or row.get("Alt") or "").strip().upper()
    return ref, alt


def build_snv_sequence_panel(row):
    start = to_int(row.get("Hg38_start") or row.get("Start"))
    end = to_int(row.get("Hg38_end") or row.get("End"))
    chrom = clean_text(row.get("Chromosome") or row.get("chromosome") or "").strip()
    if start is None or end is None or not chrom:
        return {
            "status": "unavailable",
            "message": "Sequence unavailable: missing genomic coordinates/chromosome.",
            "rows": [],
        }
    if end < start:
        start, end = end, start

    ref, alt = get_ref_alt_bases(row)
    if not (len(ref) == 1 and len(alt) == 1 and re.fullmatch(r"[ACGTN]", ref) and re.fullmatch(r"[ACGTN]", alt) and start == end):
        return {
            "status": "not_applicable",
            "message": "Sequence table shown only for SNVs with single-base REF/ALT.",
            "rows": [],
        }

    cache = load_snv_sequence_cache()
    cache_key = "|".join([clean_text(chrom).upper(), str(start), ref, alt, str(VIZ_SNV_CONTEXT_BP)])
    cached_panel = cache.get(cache_key)
    if isinstance(cached_panel, dict):
        return dict(cached_panel)

    window_start = max(1, start - VIZ_SNV_CONTEXT_BP)
    window_end = start + VIZ_SNV_CONTEXT_BP
    seq, resolved_chrom, err, seq_source = fetch_reference_window(chrom, window_start, window_end)
    if err:
        return {"status": "unavailable", "message": err, "rows": []}
    if not seq:
        return {"status": "unavailable", "message": "Reference sequence not returned.", "rows": []}

    center_idx = start - window_start
    if center_idx < 0 or center_idx >= len(seq):
        return {"status": "unavailable", "message": "Variant position outside sequence window.", "rows": []}

    left = seq[:center_idx]
    center_ref = seq[center_idx]
    right = seq[center_idx + 1 :]
    warnings = []
    if center_ref != ref:
        warnings.append("Reference mismatch: FASTA has " + center_ref + ", row REF is " + ref + ".")

    rows = [
        {"label": "Reference (5'→3')", "five_prime": left, "center": center_ref, "three_prime": right},
        {"label": "Alternate (5'→3')", "five_prime": left, "center": alt, "three_prime": right},
    ]
    meta = (
        resolved_chrom
        + ":"
        + str(window_start)
        + "-"
        + str(window_end)
        + " (variant at "
        + str(start)
        + ")"
        + (" | source: " + seq_source if seq_source else "")
    )
    panel = {
        "status": "ok",
        "message": "SNV local sequence window ±" + str(VIZ_SNV_CONTEXT_BP) + " bp",
        "meta": meta,
        "rows": rows,
        "warnings": warnings,
    }
    cache[cache_key] = panel
    save_snv_sequence_cache()
    return panel


def normalize_row_for_snv_panel(row):
    if not isinstance(row, dict):
        return {}
    out = dict(row)

    chrom = first_nonempty_ci(out, ["Chromosome", "chromosome", "grch38_chr", "grch37_chr", "chr"])
    start = first_nonempty_ci(out, ["Hg38_start", "hg38_start", "Start", "start", "grch38_pos", "grch37_pos", "position"])
    end = first_nonempty_ci(out, ["Hg38_end", "hg38_end", "End", "end", "grch38_pos", "grch37_pos", "position"])
    ref = first_nonempty_ci(out, ["Ref_allele_curated", "Ref_allele", "Ref", "ref", "grch38_ref", "grch37_ref"])
    alt = first_nonempty_ci(out, ["Alt_allele_curated", "Alt_allele", "Alt", "alt", "grch38_alt", "grch37_alt"])

    if chrom:
        out["Chromosome"] = chrom
    if start:
        out["Hg38_start"] = start
    if end:
        out["Hg38_end"] = end
    if ref:
        out["Ref_allele_curated"] = ref
    if alt:
        out["Alt_allele_curated"] = alt
    return out


def build_isbt_snv_sequence_panel(row):
    candidates = []
    raw_rows = row.get("__raw_variant_rows")
    if isinstance(raw_rows, list) and raw_rows:
        candidates.extend([item for item in raw_rows if isinstance(item, dict)])
    if isinstance(row, dict):
        candidates.append(row)

    if not candidates:
        return {
            "status": "unavailable",
            "message": "Sequence unavailable: no variant row data supplied.",
            "rows": [],
        }

    best_panel = None
    for candidate in candidates:
        panel = build_snv_sequence_panel(normalize_row_for_snv_panel(candidate))
        if panel.get("status") == "ok":
            return panel
        if best_panel is None:
            best_panel = panel

    return best_panel or {
        "status": "not_applicable",
        "message": "Sequence table shown only for SNVs with single-base REF/ALT.",
        "rows": [],
    }


def build_exon_index(rows):
    index = {}
    for row in rows:
        gene = clean_text(row.get("Gene", "")).strip()
        start = to_int(row.get("Start"))
        end = to_int(row.get("End"))
        if not gene or start is None or end is None:
            continue
        if end < start:
            start, end = end, start
        index.setdefault(gene, []).append(
            {
                "start": start,
                "end": end,
                "label": clean_text(row.get("Exon_name") or row.get("Exon_ID") or "Exon"),
                "orientation": clean_text(row.get("Orientation", "")).strip(),
            }
        )
    for gene in index:
        index[gene].sort(key=lambda r: (r["start"], r["end"]))
    return index


def build_intron_features(sorted_exons, window_start=None, window_end=None):
    features = []
    for idx in range(len(sorted_exons) - 1):
        left = sorted_exons[idx]
        right = sorted_exons[idx + 1]
        start = int(left["end"]) + 1
        end = int(right["start"]) - 1
        if end <= start:
            continue
        if window_start is not None and end < window_start:
            continue
        if window_end is not None and start > window_end:
            continue
        features.append(
            {
                "start": start,
                "end": end,
                "label": f"Intron {idx + 1}",
                "color": "#adb5bd",
                "strand": +1,
            }
        )
    return features


def build_variant_viz_bundle(row):
    start = to_int(row.get("Hg38_start") or row.get("Start"))
    end = to_int(row.get("Hg38_end") or row.get("End"))
    if start is None and end is None:
        parsed = parse_positions_from_text(row.get("Nucleotide_change") or row.get("DNA Change") or row.get("HGVS Transcript"))
        if parsed:
            start, end = parsed
    if start is None and end is not None:
        start = end
    if end is None and start is not None:
        end = start + 1
    if start is None or end is None:
        return {"features": [], "xlim": None, "notes": []}
    if end < start:
        start, end = end, start
    features = []
    label = short_label(row.get("Nucleotide_change") or row.get("DNA Change") or row.get("Variant_id") or "Variant")
    features.append({"start": start, "end": max(end, start + 1), "label": label, "color": "#d1495b", "strand": +1})

    gene = clean_text(row.get("Gene", "")).strip()
    window_start = max(0, start - 20000)
    window_end = end + 20000
    exons = EXON_INDEX_BY_GENE.get(gene, [])
    exons_in_window = []
    for exon in exons:
        if exon["end"] < window_start or exon["start"] > window_end:
            continue
        exons_in_window.append(exon)
        strand = -1 if exon.get("orientation") == "-" else +1
        features.append(
            {
                "start": exon["start"],
                "end": max(exon["end"], exon["start"] + 1),
                "label": short_label(exon["label"], 28),
                "color": "#2a9d8f",
                "strand": strand,
            }
        )

    if exons_in_window:
        features.extend(build_intron_features(exons_in_window, window_start, window_end))

    notes = [
        "Variant coordinate: " + str(start) + "-" + str(max(end, start + 1)),
        "Nucleotide change: " + clean_text(row.get("Nucleotide_change") or row.get("DNA Change") or "-"),
    ]
    if clean_text(row.get("Exon_Intron")):
        notes.append("Exon/Intron annotation: " + clean_text(row.get("Exon_Intron")))
    return {"features": features, "xlim": (window_start, window_end), "notes": notes}


def build_exon_viz_bundle(row):
    start = to_int(row.get("Start"))
    end = to_int(row.get("End"))
    if start is None or end is None:
        return {"features": [], "xlim": None, "notes": []}
    if end < start:
        start, end = end, start
    orientation = clean_text(row.get("Orientation", "")).strip()
    strand = -1 if orientation == "-" else 1
    gene = clean_text(row.get("Gene", "")).strip()
    features = []
    exons = EXON_INDEX_BY_GENE.get(gene, [])
    for exon in exons:
        exon_strand = -1 if exon.get("orientation") == "-" else +1
        is_selected = exon["start"] == start and exon["end"] == end
        features.append(
            {
                "start": exon["start"],
                "end": max(exon["end"], exon["start"] + 1),
                "label": short_label(exon["label"], 28),
                "color": "#e76f51" if is_selected else "#2a9d8f",
                "strand": exon_strand,
            }
        )
    features.extend(build_intron_features(exons))
    notes = [
        "Exon coordinate: " + str(start) + "-" + str(max(end, start + 1)),
        "Gene: " + gene,
    ]
    return {"features": features, "xlim": (max(0, start - 5000), end + 5000), "notes": notes}


def build_isbt_viz_bundle(row):
    raw_rows = row.get("__raw_variant_rows")
    if not isinstance(raw_rows, list) or not raw_rows:
        raw_rows = [row]
    features = []
    synthetic_cursor = 1
    notes = []
    for raw in raw_rows[:80]:
        dna = clean_text(raw.get("DNA Change") or raw.get("Nucleotide_change") or raw.get("HGVS Transcript") or raw.get("input"))
        parsed = parse_positions_from_text(dna)
        if parsed:
            start, end = parsed
        else:
            start, end = synthetic_cursor, synthetic_cursor + 1
            synthetic_cursor += 2
        label = short_label(dna or raw.get("isbt_allele") or raw.get("variant_id") or "ISBT change")
        features.append({"start": start, "end": max(end, start + 1), "label": label, "color": "#457b9d"})
        exon_intr = clean_text(raw.get("Exon/Intron") or raw.get("exon_intron") or "")
        if exon_intr and exon_intr != "-":
            features.append(
                {
                    "start": start,
                    "end": max(end, start + 1),
                    "label": short_label(exon_intr, 26),
                    "color": "#f4a261",
                    "strand": -1,
                }
            )
        if dna:
            notes.append(dna)
    xlim = None
    numeric_starts = [int(f["start"]) for f in features if isinstance(f.get("start"), int)]
    numeric_ends = [int(f["end"]) for f in features if isinstance(f.get("end"), int)]
    if numeric_starts and numeric_ends:
        xlim = (max(0, min(numeric_starts) - 20), max(numeric_ends) + 20)
    return {"features": features, "xlim": xlim, "notes": notes[:10]}


def encode_dna_feature_plot(features, title, xlim=None):
    if not features:
        raise ValueError("No plottable coordinates for this row.")
    min_pos = min(int(f["start"]) for f in features)
    max_pos = max(int(f["end"]) for f in features)
    offset = 0
    if min_pos < 0:
        offset = abs(min_pos) + 10

    graphic_features = []
    for feature in features:
        start = int(feature["start"]) + offset
        end = int(feature["end"]) + offset
        if end <= start:
            end = start + 1
        graphic_features.append(
            GraphicFeature(
                start=start,
                end=end,
                strand=int(feature.get("strand", +1)),
                color=feature.get("color", "#457b9d"),
                label=clean_text(feature.get("label", "")),
            )
        )

    seq_len = max(int(max_pos + offset + 10), 30)
    graphic_record = GraphicRecord(sequence_length=seq_len, features=graphic_features)
    fig_height = min(10, max(3, 1.8 + len(graphic_features) * 0.12))
    fig, ax = plt.subplots(1, 1, figsize=(12, fig_height))
    graphic_record.plot(ax=ax)
    if xlim and len(xlim) == 2:
        try:
            x0 = int(xlim[0]) + offset
            x1 = int(xlim[1]) + offset
            if x1 > x0:
                ax.set_xlim(x0, x1)
        except Exception:
            pass
    ax.set_title(clean_text(title))
    ax.set_xlabel("Position")
    fig.tight_layout()

    png_buffer = BytesIO()
    svg_buffer = BytesIO()
    fig.savefig(png_buffer, format="png", dpi=220)
    fig.savefig(svg_buffer, format="svg")
    plt.close(fig)
    return {
        "png_base64": base64.b64encode(png_buffer.getvalue()).decode("ascii"),
        "svg_base64": base64.b64encode(svg_buffer.getvalue()).decode("ascii"),
    }


def build_row_viz_payload(source, row):
    if not VIZ_AVAILABLE:
        return {
            "ok": False,
            "error": "Visualization backend unavailable. Install dna-features-viewer + matplotlib. " + VIZ_IMPORT_ERROR,
        }
    if not isinstance(row, dict):
        return {"ok": False, "error": "Invalid row payload."}

    src = clean_text(source).strip().lower()
    sequence_panel = {
        "status": "not_applicable",
        "message": "Sequence table is available for SNV variant rows.",
        "rows": [],
    }
    if src == "exon":
        bundle = build_exon_viz_bundle(row)
        features = bundle["features"]
        title = f"Exon Viz: {clean_text(row.get('Exon_name') or row.get('Exon_ID') or row.get('Gene') or '')}"
    elif src == "variant":
        bundle = build_variant_viz_bundle(row)
        features = bundle["features"]
        title = f"Variant Viz: {clean_text(row.get('Variant_id') or row.get('Nucleotide_change') or '')}"
        sequence_panel = build_snv_sequence_panel(row)
    elif src == "isbt":
        bundle = build_isbt_viz_bundle(row)
        features = bundle["features"]
        title = f"ISBT Viz: {clean_text(row.get('ISBT_Allele') or row.get('Allele_id') or row.get('Group') or '')}"
        sequence_panel = build_isbt_snv_sequence_panel(row)
    else:
        return {"ok": False, "error": "Unsupported viz source."}

    try:
        encoded = encode_dna_feature_plot(features, title, bundle.get("xlim"))
    except Exception as exc:
        return {"ok": False, "error": clean_text(exc)}

    return {
        "ok": True,
        "title": title,
        "image_png_base64": encoded.get("png_base64", ""),
        "image_svg_base64": encoded.get("svg_base64", ""),
        "notes": bundle.get("notes", []),
        "sequence_panel": sequence_panel,
    }


ROOT_DIR = Path(__file__).resolve().parent.parent
GENE_CSV = ROOT_DIR / "Gene_table.csv"
ALLELE_CSV = ROOT_DIR / "Allele_table.csv"
VARIANT_CSV = ROOT_DIR / "Variant_table.csv"
EXON_CSV = ROOT_DIR / "Exon_table.csv"
BRIDGE_CSV = ROOT_DIR / "Bridge_table.csv"

gene_table = read_csv_rows(GENE_CSV)
allele_table = read_csv_rows(ALLELE_CSV)
variant_table = read_csv_rows(VARIANT_CSV)
exon_table = read_csv_rows(EXON_CSV)
bridge_table = read_csv_rows(BRIDGE_CSV)
EXON_INDEX_BY_GENE = build_exon_index(exon_table)


def prefetch_variant_snv_sequences(max_rows=0):
    downloaded = 0
    skipped = 0
    attempted = 0
    for row in variant_table:
        start = to_int(row.get("Hg38_start") or row.get("Start"))
        end = to_int(row.get("Hg38_end") or row.get("End"))
        chrom = clean_text(row.get("Chromosome") or row.get("chromosome") or "").strip()
        ref, alt = get_ref_alt_bases(row)
        if not chrom or start is None or end is None:
            skipped += 1
            continue
        if not (start == end and len(ref) == 1 and len(alt) == 1 and re.fullmatch(r"[ACGTN]", ref) and re.fullmatch(r"[ACGTN]", alt)):
            skipped += 1
            continue
        attempted += 1
        panel = build_snv_sequence_panel(row)
        if panel.get("status") == "ok":
            downloaded += 1
        if max_rows and attempted >= max_rows:
            break
    return {"attempted": attempted, "downloaded_or_cached": downloaded, "skipped": skipped}

for row in allele_table:
    row["Allele_id"] = clean_text(row.get("Allele_id", ""))
for row in variant_table:
    row["Variant_id"] = clean_text(row.get("Variant_id", ""))
for row in bridge_table:
    row["Allele_id"] = clean_text(row.get("Allele_id", ""))
    row["Variant_id"] = clean_text(row.get("Variant_id", ""))

bridge_by_allele = {}
bridge_by_variant = {}
for row in bridge_table:
    aid = row.get("Allele_id", "")
    vid = row.get("Variant_id", "")
    if aid and vid:
        bridge_by_allele.setdefault(aid, []).append(vid)
        bridge_by_variant.setdefault(vid, []).append(aid)


def unique_keep_order(items):
    out = []
    seen = set()
    for item in items:
        if item not in seen:
            out.append(item)
            seen.add(item)
    return out


for row in allele_table:
    aid = row.get("Allele_id", "")
    linked = unique_keep_order(bridge_by_allele.get(aid, []))
    row["Variants"] = "\n".join(linked)

all_genes_set = set()
for source in (gene_table, allele_table, variant_table, exon_table):
    for row in source:
        gene = clean_text(row.get("Gene", "")).strip()
        if gene:
            all_genes_set.add(gene)
all_genes = sorted(all_genes_set)

gene_view_cols = list(gene_table[0].keys()) if gene_table else []
if "Gene_id" in gene_view_cols and "Gene" in gene_view_cols:
    gene_view_cols = [col for col in gene_view_cols if col != "Gene"]
    gene_id_idx = gene_view_cols.index("Gene_id")
    gene_view_cols.insert(gene_id_idx + 1, "Gene")
allele_view_cols = list(allele_table[0].keys()) if allele_table else []
if "Variants" not in allele_view_cols:
    allele_view_cols.append("Variants")
variant_view_cols = list(variant_table[0].keys()) if variant_table else []
exon_view_cols = list(exon_table[0].keys()) if exon_table else []

ISBT_SYSTEMS_URL = "https://api-blooddatabase.isbtweb.org/system"
ISBT_ANTIGEN_URL = "https://api-blooddatabase.isbtweb.org/antigen/search?system_symbol="
ISBT_ALLELE_URL = "https://api-blooddatabase.isbtweb.org/allele/search?system_symbol="
ISBT_RELEASE_SEARCH_URL = "https://api-blooddatabase.isbtweb.org/release/search?page=1&limit=12"
ISBT_VARIANT_VERSION_HISTORY_URL = "https://api-blooddatabase.isbtweb.org/variant/version-history"
ISBT_CACHE_PATH = ROOT_DIR / "isbt_cache.json"
ISBT_CACHE_MAX_DAYS = 30
ISBT_FLAG_COLS = ["sv_allele", "null_allele", "mod_allele", "partial_allele", "weak_allele", "el_allele"]
ISBT_GROUPED_VIEW_COLS = [
    "Group",
    "Allele_id",
    "ISBT_Allele",
    "ISBT_Phenotype",
    "Gene_name",
    "DNA Change",
    "Exon/Intron",
    "HGVS Transcript",
    "sv_allele",
    "null_allele",
    "mod_allele",
    "partial_allele",
    "weak_allele",
    "el_allele",
]
ISBT_VARIANT_VIEW_COLS = [
    "Group",
    "Allele_id",
    "ISBT_Allele",
    "Gene_name",
    "Variant_id",
    "DNA Change",
    "Exon/Intron",
    "HGVS Transcript",
    "sv_allele",
    "null_allele",
    "mod_allele",
    "partial_allele",
    "weak_allele",
    "el_allele",
]


def now_utc_epoch():
    return int(datetime.utcnow().timestamp())


def extract_records(payload):
    if isinstance(payload, list):
        return [item for item in payload if isinstance(item, dict)]
    if isinstance(payload, dict):
        for key in ("results", "data", "items", "systems", "antigens", "alleles", "variants"):
            value = payload.get(key)
            if isinstance(value, list):
                return [item for item in value if isinstance(item, dict)]
        return [payload]
    return []


def strip_html(value):
    txt = clean_text(value)
    if not txt:
        return ""
    no_tags = re.sub(r"<[^>]+>", " ", txt)
    plain = unescape(no_tags)
    return re.sub(r"\s+", " ", plain).strip()


def build_release_history_rows():
    payload = fetch_json(ISBT_RELEASE_SEARCH_URL)
    records = []
    if isinstance(payload, dict):
        data_rows = payload.get("data")
        if isinstance(data_rows, list):
            records = [row for row in data_rows if isinstance(row, dict)]
    elif isinstance(payload, list):
        records = [row for row in payload if isinstance(row, dict)]

    releases = []
    for row in records:
        release_version = to_int(row.get("releaseVersion"))
        suggestion_count = to_int(row.get("suggestionCount"))
        releases.append(
            {
                "id": to_int(row.get("id")),
                "name": clean_text(row.get("name")),
                "releaseVersion": release_version if release_version is not None else clean_text(row.get("releaseVersion")),
                "appliedAt": clean_text(row.get("appliedAt")),
                "updatedAt": clean_text(row.get("updatedAt")),
                "suggestionCount": suggestion_count if suggestion_count is not None else clean_text(row.get("suggestionCount")),
                "alleleCount": to_int(row.get("alleleCount")),
                "variantCount": to_int(row.get("variantCount")),
                "antigenCount": to_int(row.get("antigenCount")),
                "geneCount": to_int(row.get("geneCount")),
                "systemCount": to_int(row.get("systemCount")),
                "notesSummary": strip_html(row.get("notes")),
            }
        )

    releases.sort(
        key=lambda row: (
            to_int(row.get("releaseVersion")) if to_int(row.get("releaseVersion")) is not None else -1,
            clean_text(row.get("appliedAt")),
            clean_text(row.get("updatedAt")),
        ),
        reverse=True,
    )
    return releases


def flatten_record(record, prefix=""):
    out = {}
    if not isinstance(record, dict):
        return out
    for key, value in record.items():
        clean_key = clean_text(key).strip()
        if not clean_key:
            continue
        full_key = f"{prefix}_{clean_key}" if prefix else clean_key
        if isinstance(value, dict):
            nested = flatten_record(value, full_key)
            out.update(nested)
        elif isinstance(value, list):
            if all(not isinstance(item, (dict, list)) for item in value):
                out[full_key] = ", ".join(clean_text(item) for item in value)
            else:
                out[full_key] = json.dumps(value, ensure_ascii=False)
        else:
            out[full_key] = clean_text(value)
    return out


def first_nonempty_ci(record, candidate_keys):
    if not isinstance(record, dict):
        return ""
    for wanted in candidate_keys:
        wanted_low = wanted.lower()
        for key, value in record.items():
            if clean_text(key).lower() == wanted_low:
                txt = clean_text(value).strip()
                if txt:
                    return txt
    return ""


def to_bool_flag_text(value):
    txt = clean_text(value).strip()
    low = txt.lower()
    if low in ("1", "true", "t", "yes", "y"):
        return "Yes"
    if low in ("0", "false", "f", "no", "n"):
        return "No"
    return txt


def normalize_dna_change(raw_value):
    txt = clean_text(raw_value).strip()
    if not txt:
        return "-"
    if ":" in txt:
        txt = txt.split(":", 1)[1].strip()
    txt = re.sub(r"^(c\.[0-9_+\-*]+)[A-Za-z]+=$", r"\1=", txt)
    return txt or "-"


def normalize_exon_label(raw_value):
    txt = clean_text(raw_value).strip()
    if not txt:
        return ""
    low = txt.lower()
    if low.startswith("exon"):
        rest = txt[4:].strip()
        return f"Exon {rest}" if rest else "Exon"
    if re.fullmatch(r"-?\d+", txt):
        return f"Exon {txt}"
    return f"Exon {txt}"


def normalize_intron_label(raw_value):
    txt = clean_text(raw_value).strip()
    if not txt:
        return ""
    low = txt.lower()
    if low.startswith("intron"):
        rest = txt[6:].strip()
        return f"Intron {rest}" if rest else "Intron"
    return f"Intron {txt}"


def normalize_exon_intron(exon_value, intron_value):
    exon_label = normalize_exon_label(exon_value)
    intron_label = normalize_intron_label(intron_value)
    if exon_label and intron_label:
        return f"{exon_label} / {intron_label}"
    if exon_label:
        return exon_label
    if intron_label:
        return intron_label
    return "-"


def parse_c_position(dna_change):
    txt = clean_text(dna_change)
    match = re.search(r"c\.(-?\d+)", txt)
    if match:
        try:
            return int(match.group(1))
        except ValueError:
            return 10**9
    return 10**9


def parse_isbt_number(record):
    txt = first_nonempty_ci(record, ["isbt_number", "system_number", "number"])
    match = re.search(r"\d+", txt)
    if match:
        try:
            return int(match.group(0))
        except ValueError:
            return 10**9
    return 10**9


def fetch_json(url, timeout_sec=25):
    request = Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "geno-mock-python/1.0 (+https://github.com/RazaS/geno_mock_python)",
        },
    )
    with urlopen(request, timeout=timeout_sec) as response:
        body = response.read().decode("utf-8")
    return json.loads(body)


def variant_like_record(record):
    check_fields = [
        "input",
        "dna_change",
        "nucleotide_change",
        "hgvs",
        "hgvs_transcript",
        "exon",
        "intron",
        "variant_id",
    ]
    return bool(first_nonempty_ci(record, check_fields))


def dedupe_rows(rows, key_fields):
    seen = set()
    out = []
    for row in rows:
        key = tuple(clean_text(row.get(field, "")) for field in key_fields)
        if key in seen:
            continue
        seen.add(key)
        out.append(row)
    return out


def variant_matches_group(variant_row, selected_group):
    group = clean_text(selected_group).upper().strip()
    symbol = first_nonempty_ci(variant_row, ["system_symbol"]).upper()
    if group in ("RHD", "RHCE"):
        isbt_allele = first_nonempty_ci(variant_row, ["isbt_allele"]).upper()
        gene_name = first_nonempty_ci(variant_row, ["gene_name", "gene"]).upper()
        return isbt_allele.startswith(group + "*") or gene_name == group
    return symbol == group


def variant_target_groups(variant_row):
    symbol = first_nonempty_ci(variant_row, ["system_symbol"]).upper()
    if not symbol:
        return []
    if symbol == "RH":
        groups = [group for group in ("RHD", "RHCE") if variant_matches_group(variant_row, group)]
        return groups
    return [symbol]


def to_isbt_variant_view_row(raw_row, group):
    out = dict(raw_row)
    out["Group"] = group
    out["__group"] = group
    out["Allele_id"] = first_nonempty_ci(raw_row, ["allele_id", "Allele_id"])
    out["ISBT_Allele"] = first_nonempty_ci(raw_row, ["isbt_allele", "allele_name", "ISBT_Allele"])
    out["Gene_name"] = first_nonempty_ci(raw_row, ["gene_name", "gene", "Gene_name"])
    out["Variant_id"] = first_nonempty_ci(raw_row, ["variant_id", "id", "Variant_id"])
    out["DNA Change"] = first_nonempty_ci(raw_row, ["DNA Change", "dna_change", "input", "nucleotide_change", "hgvs"]) or "-"
    out["Exon/Intron"] = first_nonempty_ci(raw_row, ["Exon/Intron", "exon_intron"]) or "-"
    out["HGVS Transcript"] = first_nonempty_ci(raw_row, ["HGVS Transcript", "hgvs_transcript", "input", "hgvs"]) or "-"
    for flag_col in ISBT_FLAG_COLS:
        out[flag_col] = to_bool_flag_text(first_nonempty_ci(raw_row, [flag_col]))
    return out


def make_isbt_variant_rows(variants_rows):
    variant_rows = []
    for raw in variants_rows:
        for group in variant_target_groups(raw):
            variant_rows.append(to_isbt_variant_view_row(raw, group))

    by_variant_id = {}
    order = []

    def nonempty_score(row):
        fields = ["DNA Change", "Exon/Intron", "HGVS Transcript", "Gene_name", "Allele_id", "ISBT_Allele"]
        return sum(1 for field in fields if clean_text(row.get(field, "")).strip() and clean_text(row.get(field, "")) != "-")

    for row in variant_rows:
        variant_id = clean_text(row.get("Variant_id", "")).strip()
        if not variant_id:
            variant_id = f"missing_variant_id:{len(order) + 1}"
            row["Variant_id"] = variant_id

        if variant_id not in by_variant_id:
            by_variant_id[variant_id] = {
                "row": dict(row),
                "allele_ids": unique_keep_order([clean_text(row.get("Allele_id", ""))]),
                "allele_names": unique_keep_order([clean_text(row.get("ISBT_Allele", ""))]),
                "groups": unique_keep_order([clean_text(row.get("Group", ""))]),
                "score": nonempty_score(row),
            }
            order.append(variant_id)
            continue

        entry = by_variant_id[variant_id]
        candidate_score = nonempty_score(row)
        if candidate_score > entry["score"]:
            entry["row"] = dict(row)
            entry["score"] = candidate_score

        entry["allele_ids"] = unique_keep_order(entry["allele_ids"] + [clean_text(row.get("Allele_id", ""))])
        entry["allele_names"] = unique_keep_order(entry["allele_names"] + [clean_text(row.get("ISBT_Allele", ""))])
        entry["groups"] = unique_keep_order(entry["groups"] + [clean_text(row.get("Group", ""))])

    out_rows = []
    for variant_id in order:
        entry = by_variant_id[variant_id]
        row = dict(entry["row"])
        allele_ids = [value for value in entry["allele_ids"] if value]
        allele_names = [value for value in entry["allele_names"] if value]
        groups = [value for value in entry["groups"] if value]

        if allele_ids:
            row["Allele_id"] = ", ".join(allele_ids)
        if allele_names:
            row["ISBT_Allele"] = ", ".join(allele_names)
        if groups:
            row["Group"] = ", ".join(groups)
            row["__group"] = groups[0]
        out_rows.append(row)

    out_rows.sort(
        key=lambda row: (
            clean_text(row.get("Group", "")),
            parse_c_position(row.get("DNA Change", "")),
            clean_text(row.get("DNA Change", "")),
            clean_text(row.get("Variant_id", "")),
        )
    )
    return out_rows


def make_isbt_grouped_rows(variants_rows, allele_lookup=None):
    allele_lookup = allele_lookup if isinstance(allele_lookup, dict) else {}
    by_group = {}
    for row in variants_rows:
        for group in variant_target_groups(row):
            by_group.setdefault(group, []).append(row)

    grouped_rows = []
    for group, group_rows in by_group.items():
        by_allele = {}
        for row in group_rows:
            allele_id = first_nonempty_ci(row, ["allele_id"])
            isbt_allele = first_nonempty_ci(row, ["isbt_allele"])
            if not allele_id:
                allele_id = isbt_allele or f"{group}:{len(by_allele) + 1}"
            row_key = f"{group}|{allele_id}"
            by_allele.setdefault(row_key, []).append(row)

        for row_key, rows in by_allele.items():
            rows_sorted = sorted(
                rows,
                key=lambda r: (
                    parse_c_position(first_nonempty_ci(r, ["DNA Change", "dna_change", "input", "nucleotide_change", "hgvs"])),
                    first_nonempty_ci(r, ["DNA Change", "dna_change", "input", "nucleotide_change", "hgvs"]),
                ),
            )

            dna_lines = []
            exon_intron_lines = []
            hgvs_lines = []
            for raw in rows_sorted:
                dna_lines.append(first_nonempty_ci(raw, ["DNA Change", "dna_change", "input", "nucleotide_change", "hgvs"]) or "-")
                exon_intron_lines.append(first_nonempty_ci(raw, ["Exon/Intron", "exon_intron"]) or "-")
                hgvs_lines.append(first_nonempty_ci(raw, ["HGVS Transcript", "hgvs_transcript", "input", "hgvs"]) or "-")

            first_row = rows_sorted[0]
            allele_meta = allele_lookup.get(clean_text(first_nonempty_ci(first_row, ["allele_id"])), {})
            grouped_row = {
                "Group": group,
                "Allele_id": first_nonempty_ci(first_row, ["allele_id"]),
                "ISBT_Allele": first_nonempty_ci(first_row, ["isbt_allele", "allele_name"]),
                "ISBT_Phenotype": first_nonempty_ci(first_row, ["isbt_phenotype"]) or first_nonempty_ci(
                    allele_meta, ["isbt_phenotype"]
                ),
                "Gene_name": first_nonempty_ci(first_row, ["gene_name", "gene"]),
                "DNA Change": "\n".join(dna_lines),
                "Exon/Intron": "\n".join(exon_intron_lines),
                "HGVS Transcript": "\n".join(hgvs_lines),
                "__group": group,
                "__row_key": row_key,
                "__raw_variant_rows": rows_sorted,
                "__variant_view_rows": [to_isbt_variant_view_row(raw, group) for raw in rows_sorted],
            }
            for flag_col in ISBT_FLAG_COLS:
                grouped_row[flag_col] = first_nonempty_ci(first_row, [flag_col])
            grouped_rows.append(grouped_row)

    grouped_rows.sort(
        key=lambda row: (
            clean_text(row.get("Group", "")),
            clean_text(row.get("ISBT_Allele", "")),
            clean_text(row.get("Allele_id", "")),
        )
    )
    return grouped_rows


def parse_iso_datetime_utc(raw_value):
    txt = clean_text(raw_value).strip()
    if not txt:
        return None
    if txt.endswith("Z"):
        txt = txt[:-1] + "+00:00"
    try:
        return datetime.fromisoformat(txt)
    except ValueError:
        return None


def make_isbt_update_row(variant_row, group, event_at, operation="", actor="", version=""):
    return {
        "Group": clean_text(group),
        "Variant_id": first_nonempty_ci(variant_row, ["variant_id", "id", "Variant_id"]),
        "Allele_id": first_nonempty_ci(variant_row, ["allele_id", "Allele_id"]),
        "ISBT_Allele": first_nonempty_ci(variant_row, ["isbt_allele", "allele_name", "ISBT_Allele"]),
        "Gene_name": first_nonempty_ci(variant_row, ["gene_name", "gene", "Gene_name"]),
        "DNA Change": first_nonempty_ci(variant_row, ["DNA Change", "dna_change", "input", "nucleotide_change", "hgvs"]) or "-",
        "Exon/Intron": first_nonempty_ci(variant_row, ["Exon/Intron", "exon_intron"]) or "-",
        "HGVS Transcript": first_nonempty_ci(variant_row, ["HGVS Transcript", "hgvs_transcript", "input", "hgvs"]) or "-",
        "eventAt": clean_text(event_at),
        "operation": clean_text(operation).lower() or "update",
        "actor": clean_text(actor),
        "version": clean_text(version) or first_nonempty_ci(variant_row, ["version"]),
    }


def build_isbt_table_update_rows(variant_rows, max_events=40, use_history=True):
    if not isinstance(variant_rows, list) or not variant_rows:
        return [], []

    def row_event_epoch(row):
        dt = parse_iso_datetime_utc(first_nonempty_ci(row, ["updatedAt", "createdAt"]))
        return int(dt.timestamp()) if dt else 0

    sorted_rows = sorted(
        variant_rows,
        key=lambda row: (
            row_event_epoch(row),
            clean_text(first_nonempty_ci(row, ["variant_id", "id"])),
        ),
        reverse=True,
    )

    updates = []
    warnings = []

    if use_history:
        variant_candidates = []
        seen_variant_ids = set()
        for row in sorted_rows:
            variant_id = first_nonempty_ci(row, ["variant_id", "id"])
            if not variant_id or variant_id in seen_variant_ids:
                continue
            seen_variant_ids.add(variant_id)
            variant_candidates.append((variant_id, row))
            if len(variant_candidates) >= 30:
                break

        for variant_id, row in variant_candidates:
            history_url = (
                f"{ISBT_VARIANT_VERSION_HISTORY_URL}"
                f"?variant_id={quote(variant_id)}&all_versions=true&offset=0&limit=3"
            )
            try:
                history_payload = fetch_json(history_url)
                history_rows = extract_records(history_payload)
            except Exception as exc:
                warnings.append(f"variant_version_history:{variant_id}:{clean_text(exc)}")
                history_rows = []

            if not history_rows:
                history_rows = [{}]

            for history_row in history_rows[:2]:
                variant_data = history_row.get("variant_data") if isinstance(history_row, dict) else {}
                merged = dict(row)
                if isinstance(variant_data, dict):
                    for key, value in variant_data.items():
                        if isinstance(value, (dict, list)):
                            continue
                        merged[key] = clean_text(value)

                event_at = first_nonempty_ci(history_row, ["created_at", "createdAt", "updatedAt"]) or first_nonempty_ci(
                    merged, ["updatedAt", "createdAt"]
                )
                operation = first_nonempty_ci(history_row, ["operation"]) or "update"
                actor = first_nonempty_ci(history_row, ["createdBy", "updatedBy"])
                version = first_nonempty_ci(variant_data, ["version"]) if isinstance(variant_data, dict) else ""
                for group in variant_target_groups(merged):
                    updates.append(make_isbt_update_row(merged, group, event_at, operation=operation, actor=actor, version=version))

            if len(updates) >= max_events * 3:
                break

    if not updates:
        for row in sorted_rows[: max_events * 2]:
            event_at = first_nonempty_ci(row, ["updatedAt", "createdAt"])
            for group in variant_target_groups(row):
                updates.append(
                    make_isbt_update_row(
                        row,
                        group,
                        event_at,
                        operation="update",
                        actor="",
                        version=first_nonempty_ci(row, ["version"]),
                    )
                )

    updates = dedupe_rows(updates, ["Group", "Variant_id", "operation", "eventAt", "version"])

    def update_event_epoch(update_row):
        dt = parse_iso_datetime_utc(update_row.get("eventAt"))
        return int(dt.timestamp()) if dt else 0

    updates.sort(
        key=lambda row: (
            update_event_epoch(row),
            clean_text(row.get("Group", "")),
            clean_text(row.get("Variant_id", "")),
        ),
        reverse=True,
    )
    return updates[:max_events], warnings


def build_isbt_dataset():
    systems_payload = fetch_json(ISBT_SYSTEMS_URL)
    systems_raw = extract_records(systems_payload)

    systems = []
    symbols = []
    for record in systems_raw:
        flat = flatten_record(record)
        symbol = first_nonempty_ci(flat, ["system_symbol", "symbol"]).upper()
        if not symbol:
            continue
        flat["system_symbol"] = symbol
        systems.append(flat)
        symbols.append(symbol)

    systems = sorted(
        dedupe_rows(systems, ["system_symbol"]),
        key=lambda row: (parse_isbt_number(row), clean_text(row.get("system_symbol", ""))),
    )
    symbols = [clean_text(row.get("system_symbol", "")).upper() for row in systems]

    antigens = []
    alleles = []
    variants = []
    releases = []
    pull_errors = []

    for symbol in symbols:
        if not symbol:
            continue

        try:
            antigen_payload = fetch_json(ISBT_ANTIGEN_URL + quote(symbol))
            antigen_rows = extract_records(antigen_payload)
        except Exception as exc:
            pull_errors.append(f"antigen:{symbol}:{clean_text(exc)}")
            antigen_rows = []
        for record in antigen_rows:
            flat = flatten_record(record)
            flat["system_symbol"] = symbol
            antigens.append(flat)

        try:
            allele_payload = fetch_json(ISBT_ALLELE_URL + quote(symbol))
            allele_rows = extract_records(allele_payload)
        except Exception as exc:
            pull_errors.append(f"allele:{symbol}:{clean_text(exc)}")
            allele_rows = []
        for record in allele_rows:
            if not isinstance(record, dict):
                continue
            flat_allele = flatten_record(record)
            flat_allele["system_symbol"] = symbol

            allele_id = first_nonempty_ci(flat_allele, ["allele_id", "id"])
            if not allele_id:
                allele_id = f"{symbol}:{len(alleles) + 1}"
            flat_allele["allele_id"] = allele_id
            flat_allele["isbt_allele"] = first_nonempty_ci(flat_allele, ["isbt_allele", "allele_name", "name"])
            flat_allele["isbt_phenotype"] = first_nonempty_ci(flat_allele, ["isbt_phenotype"])
            flat_allele["gene_name"] = first_nonempty_ci(flat_allele, ["gene_name", "gene", "gene_symbol"])

            for flag_col in ISBT_FLAG_COLS:
                flat_allele[flag_col] = to_bool_flag_text(first_nonempty_ci(flat_allele, [flag_col]))

            alleles.append(flat_allele)

            nested_variants = []
            for key in ("variants", "variant", "variant_list", "allele_variants"):
                value = record.get(key)
                if isinstance(value, list):
                    nested_variants = [item for item in value if isinstance(item, dict)]
                    if nested_variants:
                        break

            if nested_variants:
                for variant_record in nested_variants:
                    flat_variant = flatten_record(variant_record)
                    flat_variant["system_symbol"] = symbol
                    flat_variant["allele_id"] = allele_id
                    if not first_nonempty_ci(flat_variant, ["isbt_allele"]):
                        flat_variant["isbt_allele"] = flat_allele.get("isbt_allele", "")
                    if not first_nonempty_ci(flat_variant, ["isbt_phenotype"]):
                        flat_variant["isbt_phenotype"] = flat_allele.get("isbt_phenotype", "")
                    if not first_nonempty_ci(flat_variant, ["gene_name", "gene"]):
                        flat_variant["gene_name"] = flat_allele.get("gene_name", "")
                    variants.append(flat_variant)
            elif variant_like_record(flat_allele):
                variants.append(dict(flat_allele))

    allele_by_id = {}
    for row in alleles:
        allele_id = clean_text(row.get("allele_id", ""))
        if not allele_id:
            continue
        if allele_id not in allele_by_id:
            allele_by_id[allele_id] = dict(row)
            continue
        existing = allele_by_id[allele_id]
        for key, value in row.items():
            if value and not clean_text(existing.get(key, "")).strip():
                existing[key] = value
    alleles = list(allele_by_id.values())

    allele_flags_by_id = {
        clean_text(row.get("allele_id", "")): {flag: clean_text(row.get(flag, "")) for flag in ISBT_FLAG_COLS}
        for row in alleles
    }

    normalized_variants = []
    for idx, row in enumerate(variants, start=1):
        out = dict(row)
        allele_id = first_nonempty_ci(out, ["allele_id"])
        if not allele_id:
            allele_id = f"missing:{idx}"
        out["allele_id"] = allele_id
        out["isbt_allele"] = first_nonempty_ci(out, ["isbt_allele", "allele_name"])
        if not first_nonempty_ci(out, ["isbt_phenotype"]):
            out["isbt_phenotype"] = first_nonempty_ci(allele_by_id.get(allele_id, {}), ["isbt_phenotype"])
        out["gene_name"] = first_nonempty_ci(out, ["gene_name", "gene", "gene_symbol"])
        out["variant_id"] = first_nonempty_ci(out, ["variant_id", "id"]) or f"v:{idx}"

        dna_change = normalize_dna_change(first_nonempty_ci(out, ["input", "dna_change", "nucleotide_change", "hgvs", "hgvs_c"]))
        exon_value = first_nonempty_ci(out, ["exon", "exon_number", "exon_no"])
        intron_value = first_nonempty_ci(out, ["intron", "intron_number", "intron_no"])
        exon_intron = normalize_exon_intron(exon_value, intron_value)
        hgvs_transcript = clean_text(first_nonempty_ci(out, ["input", "hgvs_transcript", "hgvs"])).strip() or "-"

        out["DNA Change"] = dna_change
        out["Exon/Intron"] = exon_intron
        out["HGVS Transcript"] = hgvs_transcript

        allele_flags = allele_flags_by_id.get(allele_id, {})
        for flag_col in ISBT_FLAG_COLS:
            if not first_nonempty_ci(out, [flag_col]):
                out[flag_col] = allele_flags.get(flag_col, "")
            out[flag_col] = to_bool_flag_text(out.get(flag_col, ""))

        normalized_variants.append(out)

    variants = dedupe_rows(
        normalized_variants,
        ["system_symbol", "allele_id", "variant_id", "DNA Change", "Exon/Intron", "HGVS Transcript"],
    )
    antigens = dedupe_rows(antigens, ["system_symbol", "antigen_id", "id", "name"])

    group_symbols = [clean_text(row.get("system_symbol", "")).upper() for row in systems if clean_text(row.get("system_symbol", ""))]
    groups = []
    for symbol in group_symbols:
        if symbol == "RH":
            groups.extend(["RHD", "RHCE"])
        else:
            groups.append(symbol)
    groups = unique_keep_order(groups)

    grouped_rows = make_isbt_grouped_rows(variants, allele_lookup=allele_by_id)
    variant_rows = make_isbt_variant_rows(variants)
    table_updates, revision_warnings = build_isbt_table_update_rows(variants, max_events=40, use_history=True)
    pull_errors.extend(revision_warnings)
    raw_export_columns = sorted(
        {key for row in (variants + variant_rows) for key in row.keys() if not key.startswith("__")}
    )

    try:
        releases = build_release_history_rows()
    except Exception as exc:
        pull_errors.append(f"release_search:{clean_text(exc)}")

    meta = {
        "updated_at": datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC"),
        "updated_epoch": now_utc_epoch(),
        "source_urls": [
            ISBT_SYSTEMS_URL,
            ISBT_ANTIGEN_URL + "<SYMBOL>",
            ISBT_ALLELE_URL + "<SYMBOL>",
            ISBT_RELEASE_SEARCH_URL,
            ISBT_VARIANT_VERSION_HISTORY_URL + "?variant_id=<ID>&all_versions=true",
        ],
        "row_counts": {
            "systems": len(systems),
            "antigens": len(antigens),
            "alleles": len(alleles),
            "variants": len(variants),
            "grouped_rows": len(grouped_rows),
            "variant_rows": len(variant_rows),
            "table_updates": len(table_updates),
            "releases": len(releases),
        },
        "pull_errors": pull_errors,
    }

    return {
        "metadata": meta,
        "systems": systems,
        "antigens": antigens,
        "alleles": alleles,
        "variants": variants,
        "releases": releases,
        "variant_rows": variant_rows,
        "table_updates": table_updates,
        "groups": groups,
        "grouped_rows": grouped_rows,
        "columns": {
            "grouped_view": ISBT_GROUPED_VIEW_COLS,
            "variant_view": ISBT_VARIANT_VIEW_COLS,
            "raw_export": raw_export_columns,
        },
    }


def empty_isbt_dataset(error_message=""):
    return {
        "metadata": {
            "updated_at": "Not available",
            "updated_epoch": 0,
            "source_urls": [
                ISBT_SYSTEMS_URL,
                ISBT_ANTIGEN_URL + "<SYMBOL>",
                ISBT_ALLELE_URL + "<SYMBOL>",
                ISBT_RELEASE_SEARCH_URL,
                ISBT_VARIANT_VERSION_HISTORY_URL + "?variant_id=<ID>&all_versions=true",
            ],
            "row_counts": {
                "systems": 0,
                "antigens": 0,
                "alleles": 0,
                "variants": 0,
                "grouped_rows": 0,
                "variant_rows": 0,
                "table_updates": 0,
                "releases": 0,
            },
            "pull_errors": [],
            "error": clean_text(error_message),
        },
        "systems": [],
        "antigens": [],
        "alleles": [],
        "variants": [],
        "releases": [],
        "variant_rows": [],
        "table_updates": [],
        "groups": [],
        "grouped_rows": [],
        "columns": {"grouped_view": ISBT_GROUPED_VIEW_COLS, "variant_view": ISBT_VARIANT_VIEW_COLS, "raw_export": []},
    }


def ensure_isbt_dataset_shape(dataset):
    if not isinstance(dataset, dict):
        return empty_isbt_dataset("Invalid ISBT dataset shape")

    dataset.setdefault("systems", [])
    dataset.setdefault("antigens", [])
    dataset.setdefault("alleles", [])
    dataset.setdefault("variants", [])
    dataset.setdefault("releases", [])
    dataset.setdefault("groups", [])
    dataset.setdefault("grouped_rows", [])

    allele_lookup = {}
    for row in dataset.get("alleles", []):
        allele_id = clean_text(row.get("allele_id", ""))
        if allele_id and allele_id not in allele_lookup:
            allele_lookup[allele_id] = dict(row)

    # Always regenerate grouped_rows so newly added derived columns are available
    # immediately even when the cache was created by an older app version.
    dataset["grouped_rows"] = make_isbt_grouped_rows(dataset.get("variants", []), allele_lookup=allele_lookup)
    if not isinstance(dataset.get("groups"), list):
        dataset["groups"] = []
    # Always regenerate variant_rows to keep one-row-per-Variant_id guarantees
    # even when older cache payloads already contain stale/duplicated variant_rows.
    dataset["variant_rows"] = make_isbt_variant_rows(dataset.get("variants", []))
    if not isinstance(dataset.get("table_updates"), list):
        updates, _ = build_isbt_table_update_rows(dataset.get("variants", []), max_events=40, use_history=False)
        dataset["table_updates"] = updates

    columns = dataset.setdefault("columns", {})
    columns["grouped_view"] = ISBT_GROUPED_VIEW_COLS
    if not isinstance(columns.get("variant_view"), list):
        columns["variant_view"] = ISBT_VARIANT_VIEW_COLS
    if not isinstance(columns.get("raw_export"), list):
        variants_for_export = list(dataset.get("variants", [])) + list(dataset.get("variant_rows", []))
        columns["raw_export"] = sorted(
            {key for row in variants_for_export for key in row.keys() if not key.startswith("__")}
        )

    meta = dataset.setdefault("metadata", {})
    meta.setdefault("updated_at", "Not available")
    meta.setdefault("updated_epoch", 0)
    source_urls = meta.setdefault(
        "source_urls",
        [
            ISBT_SYSTEMS_URL,
            ISBT_ANTIGEN_URL + "<SYMBOL>",
            ISBT_ALLELE_URL + "<SYMBOL>",
            ISBT_RELEASE_SEARCH_URL,
            ISBT_VARIANT_VERSION_HISTORY_URL + "?variant_id=<ID>&all_versions=true",
        ],
    )
    if isinstance(source_urls, list):
        version_history_url = ISBT_VARIANT_VERSION_HISTORY_URL + "?variant_id=<ID>&all_versions=true"
        if version_history_url not in source_urls:
            source_urls.append(version_history_url)
    meta.setdefault("pull_errors", [])
    row_counts = meta.setdefault("row_counts", {})
    row_counts["systems"] = len(dataset.get("systems", []))
    row_counts["antigens"] = len(dataset.get("antigens", []))
    row_counts["alleles"] = len(dataset.get("alleles", []))
    row_counts["variants"] = len(dataset.get("variants", []))
    row_counts["grouped_rows"] = len(dataset.get("grouped_rows", []))
    row_counts["variant_rows"] = len(dataset.get("variant_rows", []))
    row_counts["table_updates"] = len(dataset.get("table_updates", []))
    row_counts["releases"] = len(dataset.get("releases", []))
    return dataset


def load_isbt_dataset():
    cached = None
    if ISBT_CACHE_PATH.exists():
        try:
            with ISBT_CACHE_PATH.open("r", encoding="utf-8") as handle:
                cached = json.load(handle)
        except Exception:
            cached = None

    if cached:
        had_variant_rows = isinstance(cached.get("variant_rows"), list)
        had_table_updates = isinstance(cached.get("table_updates"), list)
        had_variant_view_cols = isinstance(cached.get("columns", {}).get("variant_view"), list)
        cache_shape_updated = False
        updated_epoch = int(cached.get("metadata", {}).get("updated_epoch", 0) or 0)
        if updated_epoch > 0:
            age_days = (now_utc_epoch() - updated_epoch) / 86400.0
            if age_days <= ISBT_CACHE_MAX_DAYS:
                cached_releases = cached.get("releases")
                if not isinstance(cached_releases, list) or not cached_releases:
                    try:
                        cached["releases"] = build_release_history_rows()
                        cached_meta = cached.setdefault("metadata", {})
                        cached_counts = cached_meta.setdefault("row_counts", {})
                        cached_counts["releases"] = len(cached.get("releases", []))
                        cache_shape_updated = True
                    except Exception as exc:
                        cached_meta = cached.setdefault("metadata", {})
                        cached_meta["release_feed_error"] = clean_text(exc)
                cached = ensure_isbt_dataset_shape(cached)
                if (not had_variant_rows) or (not had_table_updates) or (not had_variant_view_cols) or cache_shape_updated:
                    try:
                        with ISBT_CACHE_PATH.open("w", encoding="utf-8") as handle:
                            json.dump(cached, handle, ensure_ascii=False)
                    except OSError as exc:
                        cached_meta = cached.setdefault("metadata", {})
                        cached_meta["cache_write_error"] = clean_text(exc)
                return cached

    try:
        fresh = build_isbt_dataset()
    except (URLError, TimeoutError, json.JSONDecodeError, OSError, ValueError) as exc:
        if cached:
            cached_meta = cached.setdefault("metadata", {})
            cached_meta["refresh_error"] = clean_text(exc)
            return ensure_isbt_dataset_shape(cached)
        return ensure_isbt_dataset_shape(empty_isbt_dataset(clean_text(exc)))

    fresh = ensure_isbt_dataset_shape(fresh)

    # Cache write failures should not discard freshly pulled data.
    try:
        with ISBT_CACHE_PATH.open("w", encoding="utf-8") as handle:
            json.dump(fresh, handle, ensure_ascii=False)
    except OSError as exc:
        fresh_meta = fresh.setdefault("metadata", {})
        fresh_meta["cache_write_error"] = clean_text(exc)
    return fresh


isbt_data = load_isbt_dataset()

snv_prefetch_summary = {}
if clean_text(os.environ.get("PREFETCH_SNV_SEQUENCE_CACHE", "0")).strip() == "1":
    try:
        snv_prefetch_summary = prefetch_variant_snv_sequences()
    except Exception as exc:
        snv_prefetch_summary = {"error": clean_text(exc)}

INIT_PAYLOAD = {
    "gene_table": gene_table,
    "allele_table": allele_table,
    "variant_table": variant_table,
    "exon_table": exon_table,
    "bridge_by_allele": bridge_by_allele,
    "bridge_by_variant": bridge_by_variant,
    "all_genes": all_genes,
    "counts": {
        "gene": len(gene_table),
        "allele": len(allele_table),
        "variant": len(variant_table),
        "exon": len(exon_table),
        "bridge": len(bridge_table),
        "isbt_grouped": len(isbt_data.get("grouped_rows", [])),
        "isbt_variant_rows": len(isbt_data.get("variant_rows", [])),
    },
    "columns": {
        "gene_table": gene_view_cols,
        "allele_view": allele_view_cols,
        "variant_view": variant_view_cols,
        "exon_view": exon_view_cols,
    },
    "isbt": isbt_data,
    "snv_sequence_cache": {
        "enabled": True,
        "cache_path": str(SNV_SEQUENCE_CACHE_PATH),
        "cached_entries": len(load_snv_sequence_cache()),
        "prefetch_summary": snv_prefetch_summary,
    },
}

FEEDBACK_LOG = []

INDEX_HTML = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>Geno Database (Python)</title>
  <style>
    html, body {
      margin: 0;
      padding: 0;
      height: 100%;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
      background: #f3f4f6;
      color: #1f2937;
    }
    * { box-sizing: border-box; }
    .app-shell {
      min-height: 100vh;
      padding: 12px;
      display: flex;
      flex-direction: column;
      gap: 12px;
    }
    .app-root {
      display: flex;
      gap: 12px;
      min-height: calc(100vh - 36px);
    }
    .panel-a {
      flex: 0 0 80%;
      min-width: 0;
      display: flex;
      flex-direction: column;
      gap: 12px;
    }
    .panel-b {
      flex: 0 0 20%;
      min-width: 0;
      min-height: 0;
      height: 100%;
      display: flex;
      flex-direction: column;
      gap: 12px;
      overflow: hidden;
    }
    .app-root.panel-b-collapsed .panel-a {
      flex: 1 1 100%;
    }
    .app-root.panel-b-collapsed .panel-b {
      display: none;
    }
    .section-box {
      border: 1px solid #c8c8c8;
      border-radius: 8px;
      background: #fff;
      padding: 10px;
      min-width: 0;
    }
    .a1-box {
      min-height: 78px;
      display: flex;
      align-items: center;
    }
    .a2-box {
      flex: 1 1 auto;
      min-height: 0;
      min-width: 0;
      display: flex;
      flex-direction: column;
      overflow: hidden;
    }
    .mode-row {
      display: flex;
      gap: 8px;
      width: 100%;
      align-items: center;
      flex-wrap: wrap;
    }
    .sidebar-toggle-wrap {
      margin-left: auto;
      display: flex;
      gap: 8px;
      align-items: center;
    }
    .toggle-btn {
      border: 1px solid #b8b8b8;
      border-radius: 6px;
      background: #f7f7f7;
      padding: 10px 14px;
      font-size: 15px;
      cursor: pointer;
      min-width: 120px;
    }
    .toggle-btn.active {
      background: #dbeaf8;
      border-color: #6a9ecf;
      font-weight: 600;
    }
    .download-actions {
      display: flex;
      gap: 8px;
      align-items: center;
      flex-wrap: wrap;
    }
    .download-note {
      margin-top: 6px;
      font-size: 0.92em;
      color: #555;
    }
    .b1-box {
      flex: 0 0 26%;
      min-height: 120px;
      max-height: 220px;
      display: flex;
      flex-direction: column;
      overflow: hidden;
    }
    .b1-box h4 {
      margin: 0 0 6px 0;
      font-size: 14px;
    }
    #latest-updates {
      font-size: 12px;
      line-height: 1.3;
      overflow-y: auto;
      min-height: 0;
      padding-right: 4px;
    }
    #latest-updates ul {
      margin: 0;
      padding-left: 16px;
    }
    #latest-updates li {
      margin-bottom: 4px;
    }
    .b2-box {
      flex: 1 1 auto;
      min-height: 0;
      overflow-y: auto;
    }
    .a2-allele-layout {
      display: flex;
      gap: 12px;
      flex: 1 1 auto;
      height: 100%;
      min-height: 0;
      overflow: hidden;
    }
    .a2-allele-a {
      flex: 0 0 24%;
      display: grid;
      grid-template-rows: auto minmax(0, 1fr);
      height: 100%;
      min-height: 0;
      overflow: hidden;
    }
    .gene-search-wrap {
      flex: 0 0 auto;
      min-height: 0;
    }
    .a2-allele-b {
      flex: 1 1 76%;
      min-width: 0;
      height: 100%;
      min-height: 0;
      display: flex;
      flex-direction: column;
      overflow: hidden;
    }
    .gene-search-wrap label,
    .detail-controls-grid label,
    .isbt-controls-grid label {
      display: block;
      margin-bottom: 4px;
      font-weight: 600;
    }
    .input {
      width: 100%;
      border: 1px solid #c5c9cf;
      border-radius: 6px;
      padding: 8px 10px;
      font-size: 14px;
      background: #fff;
    }
    .gene-list-wrap {
      margin-top: 6px;
      border: 1px solid #e2e2e2;
      border-radius: 6px;
      background: #fbfbfb;
      overflow: auto;
      display: flex;
      flex-direction: column;
      flex: 1 1 auto;
      min-height: 0;
    }
    .detail-controls {
      flex: 0 0 auto;
      margin-bottom: 8px;
    }
    .detail-controls-grid {
      display: grid;
      grid-template-columns: minmax(220px, 1fr) 140px 140px 170px;
      gap: 8px;
      align-items: end;
    }
    .a2-allele-layout.gene-selector-collapsed .a2-allele-a {
      display: none;
    }
    .a2-allele-layout.gene-selector-collapsed .a2-allele-b {
      flex: 1 1 100%;
    }
    .a2-isbt-layout {
      display: flex;
      gap: 12px;
      flex: 1 1 auto;
      height: 100%;
      min-height: 0;
      overflow: hidden;
    }
    .a2-isbt-a {
      flex: 0 0 24%;
      display: grid;
      grid-template-rows: auto minmax(0, 1fr);
      height: 100%;
      min-height: 0;
      overflow: hidden;
    }
    .a2-isbt-b {
      flex: 1 1 76%;
      min-width: 0;
      height: 100%;
      min-height: 0;
      display: flex;
      flex-direction: column;
      overflow: hidden;
    }
    .isbt-controls-grid {
      display: grid;
      grid-template-columns: minmax(220px, 1fr) minmax(230px, 1fr) minmax(220px, 1fr) 190px;
      gap: 8px;
      align-items: end;
    }
    .subtable-toggle {
      display: flex;
      gap: 6px;
      align-items: center;
      flex-wrap: wrap;
    }
    .btn.active {
      background: #dbeaf8;
      border-color: #6a9ecf;
      font-weight: 600;
    }
    .a2-isbt-layout.group-selector-collapsed .a2-isbt-a {
      display: none;
    }
    .a2-isbt-layout.group-selector-collapsed .a2-isbt-b {
      flex: 1 1 100%;
    }
    .detail-table-wrap {
      flex: 1 1 auto;
      min-height: 0;
      overflow: hidden;
    }
    .table-shell {
      border: 1px solid #d9dce1;
      border-radius: 6px;
      background: #fff;
      overflow-x: auto;
      overflow-y: auto;
      width: 100%;
      min-width: 0;
    }
    .table-top-scroll {
      border: 1px solid #d9dce1;
      border-radius: 6px;
      background: #fff;
      overflow-x: auto;
      overflow-y: hidden;
      height: 14px;
      margin-bottom: 6px;
    }
    .table-top-scroll-inner {
      height: 1px;
    }
    table.data-table {
      width: max-content;
      min-width: 100%;
      border-collapse: collapse;
      font-size: 14px;
      background: #fff;
    }
    table.data-table th,
    table.data-table td {
      border-bottom: 1px solid #e5e7eb;
      padding: 8px 10px;
      text-align: left;
      vertical-align: top;
      white-space: nowrap;
    }
    table.data-table thead th {
      position: sticky;
      top: 0;
      z-index: 2;
      background: #f8fafc;
      font-weight: 700;
      border-bottom: 1px solid #d1d5db;
    }
    table.data-table td.multiline {
      white-space: pre-line;
    }
    table.data-table th.allele-name-col,
    table.data-table td.allele-name-col {
      width: 25ch;
      min-width: 25ch;
      max-width: 25ch;
      white-space: normal;
      overflow-wrap: anywhere;
      word-break: break-word;
    }
    table.data-table th.tight-wrap-col,
    table.data-table td.tight-wrap-col {
      width: 16ch;
      min-width: 16ch;
      max-width: 16ch;
      white-space: normal;
      overflow-wrap: anywhere;
      word-break: break-word;
    }
    table.data-table tr.selected-row {
      background: #dbeafe;
    }
    table.data-table tr:hover {
      background: #f3f6fb;
    }
    .check-cell {
      width: 34px;
      text-align: center;
    }
    .viz-cell {
      width: 72px;
      text-align: center;
    }
    .viz-btn {
      border: 1px solid #1d4ed8;
      border-radius: 5px;
      background: #eff6ff;
      color: #1e3a8a;
      font-size: 12px;
      padding: 4px 8px;
      cursor: pointer;
    }
    .viz-btn:hover {
      background: #dbeafe;
    }
    .hidden {
      display: none !important;
    }
    .last-feedback {
      margin-top: 10px;
      border-top: 1px solid #e6e6e6;
      padding-top: 8px;
      font-size: 0.92em;
    }
    .muted {
      color: #6b7280;
    }
    .modal-overlay {
      position: fixed;
      inset: 0;
      background: rgba(31, 41, 55, 0.45);
      display: none;
      align-items: center;
      justify-content: center;
      z-index: 1000;
      padding: 18px;
    }
    .modal-overlay.open {
      display: flex;
    }
    .modal-box {
      width: min(1200px, 95vw);
      max-height: 90vh;
      background: #fff;
      border-radius: 10px;
      border: 1px solid #d1d5db;
      display: flex;
      flex-direction: column;
      overflow: hidden;
    }
    .modal-head {
      padding: 12px 14px;
      border-bottom: 1px solid #e5e7eb;
      font-weight: 700;
    }
    .modal-body {
      padding: 12px;
      min-height: 0;
      overflow: hidden;
      flex: 1 1 auto;
    }
    .modal-footer {
      padding: 10px 12px;
      border-top: 1px solid #e5e7eb;
      display: flex;
      gap: 8px;
      justify-content: flex-end;
      flex-wrap: wrap;
    }
    .help-content {
      font-size: 14px;
      line-height: 1.4;
    }
    .help-content ul {
      margin: 0;
      padding-left: 18px;
    }
    .help-content li {
      margin-bottom: 8px;
    }
    .btn {
      border: 1px solid #b8b8b8;
      border-radius: 6px;
      background: #f7f7f7;
      padding: 8px 12px;
      font-size: 14px;
      cursor: pointer;
    }
    .btn-primary {
      background: #2563eb;
      color: #fff;
      border-color: #1d4ed8;
    }
    .status-ok {
      color: #047857;
      font-size: 13px;
      margin-top: 8px;
    }
    .viz-image-wrap {
      width: 100%;
      height: 100%;
      overflow: auto;
      border: 1px solid #e5e7eb;
      border-radius: 8px;
      background: #f8fafc;
      padding: 8px;
    }
    .viz-image-wrap img {
      width: auto;
      max-width: none;
      height: auto;
      display: block;
      background: #fff;
      border-radius: 6px;
      transform-origin: top left;
    }
    .viz-toolbar {
      display: flex;
      gap: 8px;
      align-items: center;
      margin-right: auto;
    }
    .viz-scale-label {
      font-size: 12px;
      color: #475569;
      min-width: 44px;
      text-align: right;
    }
    .viz-notes {
      margin-bottom: 8px;
      font-size: 13px;
      color: #334155;
      border: 1px solid #e2e8f0;
      background: #f8fafc;
      border-radius: 6px;
      padding: 6px 8px;
      line-height: 1.3;
    }
    .viz-seq-panel {
      margin-top: 10px;
      border: 1px solid #dbe3ee;
      border-radius: 6px;
      background: #ffffff;
      padding: 8px;
    }
    .viz-seq-meta {
      font-size: 12px;
      color: #334155;
      margin-bottom: 6px;
    }
    .viz-seq-table-wrap {
      max-height: 180px;
      overflow: auto;
      border: 1px solid #e2e8f0;
      border-radius: 5px;
      background: #f8fafc;
    }
    .viz-seq-table {
      width: 100%;
      border-collapse: collapse;
      font-family: "SFMono-Regular", Menlo, Monaco, Consolas, "Liberation Mono", monospace;
      font-size: 12px;
    }
    .viz-seq-table th,
    .viz-seq-table td {
      border-bottom: 1px solid #e2e8f0;
      padding: 6px 8px;
      text-align: left;
      vertical-align: top;
      white-space: nowrap;
    }
    .viz-seq-center {
      font-weight: 700;
      color: #7f1d1d;
      background: #fee2e2;
      border-radius: 3px;
    }
    .viz-seq-warn {
      margin-top: 6px;
      font-size: 12px;
      color: #92400e;
    }
    .viz-error {
      color: #991b1b;
      font-size: 14px;
      padding: 10px;
      border: 1px solid #fecaca;
      background: #fef2f2;
      border-radius: 6px;
    }
    @media (max-width: 1100px) {
      .app-root {
        flex-direction: column;
      }
      .panel-a, .panel-b {
        flex: 1 1 auto;
      }
      .a2-allele-layout {
        flex-direction: column;
      }
      .a2-allele-a, .a2-allele-b, .a2-isbt-a, .a2-isbt-b {
        max-height: none;
      }
      .detail-controls-grid, .isbt-controls-grid {
        grid-template-columns: 1fr;
      }
      .gene-list-wrap {
        max-height: 40vh;
      }
    }
  </style>
</head>
<body>
  <div class="app-shell">
    <div class="app-root">
      <div class="panel-a">
        <div class="section-box a1-box">
          <div class="mode-row">
            <strong style="font-size: 18px;">ISBT Table</strong>
            <div class="sidebar-toggle-wrap">
              <button id="nav-help-btn" class="btn">Navigation Help</button>
              <button id="toggle-panel-b" class="btn">Hide Sidebar</button>
            </div>
          </div>
        </div>
        <div class="section-box a2-box">
          <div id="isbt-mode">
            <div id="isbt-layout" class="a2-isbt-layout">
              <div class="a2-isbt-a">
                <div class="gene-search-wrap">
                  <label for="isbt-group-search">Search groups</label>
                  <input id="isbt-group-search" class="input" type="text" placeholder="Type to filter ISBT groups" />
                </div>
                <div class="gene-list-wrap">
                  <div id="isbt-group-selector-table"></div>
                </div>
              </div>
              <div class="a2-isbt-b">
                <div class="detail-controls">
                  <div class="isbt-controls-grid">
                    <div>
                      <label for="isbt-table-search">Search ISBT table</label>
                      <input id="isbt-table-search" class="input" type="text" placeholder="Filter grouped allele rows" />
                    </div>
                    <div>
                      <label>ISBT subtable</label>
                      <div class="subtable-toggle">
                        <button id="isbt-view-allele-btn" class="btn active" type="button">Allele</button>
                        <button id="isbt-view-variant-btn" class="btn" type="button">Variant</button>
                      </div>
                    </div>
                    <div id="isbt-last-updated" class="download-note"></div>
                    <button id="toggle-isbt-selector" class="btn">Hide Group Selector</button>
                  </div>
                </div>
                <div class="detail-table-wrap">
                  <div id="isbt-detail-table"></div>
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>

      <div class="panel-b">
        <div class="section-box">
          <div class="download-actions">
            <button id="download-main-btn" class="btn">Download CSV</button>
            <button id="clear-main-btn" class="btn">Clear All</button>
          </div>
          <div id="download-context" class="download-note"></div>
        </div>
        <div class="section-box b1-box">
          <h4>Latest updates</h4>
          <div id="latest-updates"></div>
        </div>
        <div class="section-box b2-box">
          <h4>Feedback form</h4>
          <label for="feedback-state">State</label>
          <input id="feedback-state" class="input" type="text" readonly />
          <div style="height: 8px;"></div>
          <label for="feedback-email">Email</label>
          <input id="feedback-email" class="input" type="text" />
          <div style="height: 8px;"></div>
          <label for="feedback-comment">Comment</label>
          <textarea id="feedback-comment" class="input" rows="5"></textarea>
          <div style="height: 10px;"></div>
          <button id="feedback-submit" class="btn btn-primary">Submit</button>
          <div id="feedback-flash" class="status-ok"></div>
          <div id="feedback-status" class="last-feedback muted">No submissions yet.</div>
        </div>
      </div>
    </div>
  </div>

  <div id="relation-modal" class="modal-overlay">
    <div class="modal-box">
      <div id="modal-title" class="modal-head"></div>
      <div class="modal-body">
        <div id="modal-table"></div>
      </div>
      <div class="modal-footer">
        <button id="modal-download-btn" class="btn">Download CSV</button>
        <button id="modal-clear-btn" class="btn">Clear All</button>
        <button id="modal-close-btn" class="btn">Close</button>
      </div>
    </div>
  </div>

  <div id="viz-modal" class="modal-overlay">
    <div class="modal-box">
      <div id="viz-title" class="modal-head">Visualization</div>
      <div class="modal-body">
        <div id="viz-content" class="viz-image-wrap"></div>
      </div>
      <div class="modal-footer">
        <div class="viz-toolbar">
          <button id="viz-zoom-out-btn" class="btn">-</button>
          <button id="viz-zoom-in-btn" class="btn">+</button>
          <button id="viz-zoom-reset-btn" class="btn">Reset</button>
          <span id="viz-scale-label" class="viz-scale-label">100%</span>
        </div>
        <button id="viz-close-btn" class="btn">Close</button>
      </div>
    </div>
  </div>

  <div id="help-modal" class="modal-overlay">
    <div class="modal-box" style="width: min(900px, 95vw);">
      <div class="modal-head">Navigation Instructions</div>
      <div class="modal-body">
        <div class="help-content">
          <ul>
            <li>Use the left <strong>Group</strong> selector to choose an ISBT system group, then use the top <strong>ISBT subtable</strong> buttons to switch between <strong>Allele</strong> and <strong>Variant</strong> tables.</li>
            <li>Search supports <strong>AND/OR</strong> logic: use commas for <strong>OR</strong> and semicolons for <strong>AND</strong>. Example: <code>ABO,RH;Exon 7</code> means <em>(ABO OR RH) AND Exon 7</em>.</li>
            <li>Double-click allele or variant rows to open the linked popup table (allele to variants, variant to alleles). Single-click is also supported.</li>
            <li>Use row checkboxes to select entries, then click <strong>Download CSV</strong> to export selected rows from the currently active table. Use <strong>Clear All</strong> to clear selections.</li>
          </ul>
        </div>
      </div>
      <div class="modal-footer">
        <button id="help-close-btn" class="btn">Close</button>
      </div>
    </div>
  </div>

  <script>
    let DATA = null;

    const state = {
      mode: "IT",
      selectedGene: "",
      geneSearch: "",
      granularity: "Allele",
      tableSearch: "",
      selectedIsbtGroup: "",
      isbtGroupSearch: "",
      isbtTableSearch: "",
      vizScale: 1.0,
      panelBCollapsed: false,
      alleleSelectorCollapsed: false,
      isbtSelectorCollapsed: false,
      activeTable: "isbtDetail",
      checked: {
        gene: new Set(),
        exon: new Set(),
        detail: new Set(),
        isbtDetail: new Set(),
        modalVariants: new Set(),
        modalAlleles: new Set()
      },
      modal: {
        open: false,
        type: null,
        rows: [],
        title: ""
      }
    };

    function text(v) {
      if (v === null || v === undefined) return "";
      return String(v);
    }

    function ciContains(haystack, needle) {
      return text(haystack).toLowerCase().includes(text(needle).toLowerCase());
    }

    function parseSearchGroups(rawQuery) {
      const query = text(rawQuery).trim();
      if (!query) return [];
      return query
        .split(";")
        .map((andPart) => andPart.trim())
        .filter(Boolean)
        .map((andPart) => andPart.split(",")
          .map((orTerm) => orTerm.trim())
          .filter(Boolean));
    }

    function matchesSearchGroups(haystack, rawQuery) {
      const groups = parseSearchGroups(rawQuery);
      if (groups.length === 0) return true;
      const source = text(haystack).toLowerCase();
      return groups.every((orGroup) => orGroup.some((term) => source.includes(term.toLowerCase())));
    }

    function uniqueOrdered(list) {
      const out = [];
      const seen = new Set();
      list.forEach((item) => {
        const key = text(item);
        if (!seen.has(key)) {
          seen.add(key);
          out.push(key);
        }
      });
      return out;
    }

    function normalizeRow(row, columns) {
      const out = {};
      columns.forEach((col) => {
        out[col] = row[col] === undefined || row[col] === null ? "" : row[col];
      });
      return out;
    }

    function getFilteredGenes() {
      if (!state.geneSearch.trim()) return DATA.all_genes.slice();
      return DATA.all_genes.filter((g) => matchesSearchGroups(g, state.geneSearch));
    }

    function ensureSelectedGene() {
      const genes = getFilteredGenes();
      if (genes.length === 0) {
        state.selectedGene = "";
        return;
      }
      if (!genes.includes(state.selectedGene)) {
        state.selectedGene = genes[0];
        state.checked.detail.clear();
      }
    }

    function detailBaseRows() {
      if (!state.selectedGene) return [];
      const source = state.granularity === "Allele" ? DATA.allele_table : DATA.variant_table;
      const rows = source.filter((r) => text(r.Gene) === state.selectedGene);
      if (!state.tableSearch.trim()) return rows;
      return rows.filter((row) => {
        const joined = Object.values(row).map(text).join(" || ");
        return matchesSearchGroups(joined, state.tableSearch);
      });
    }

    function getFilteredIsbtGroups() {
      const groups = (DATA.isbt && Array.isArray(DATA.isbt.groups)) ? DATA.isbt.groups : [];
      if (!state.isbtGroupSearch.trim()) return groups.slice();
      return groups.filter((group) => matchesSearchGroups(group, state.isbtGroupSearch));
    }

    function ensureSelectedIsbtGroup() {
      const groups = getFilteredIsbtGroups();
      if (groups.length === 0) {
        state.selectedIsbtGroup = "";
        return;
      }
      if (!groups.includes(state.selectedIsbtGroup)) {
        state.selectedIsbtGroup = groups[0];
        state.checked.isbtDetail.clear();
      }
    }

    function isbtDetailBaseRows() {
      if (!state.selectedIsbtGroup || !DATA.isbt) return [];
      const source = state.granularity === "Allele"
        ? (Array.isArray(DATA.isbt.grouped_rows) ? DATA.isbt.grouped_rows : [])
        : (Array.isArray(DATA.isbt.variant_rows) ? DATA.isbt.variant_rows : []);
      const rows = source.filter((row) => {
        const rowGroup = text(row.__group || row.Group);
        return rowGroup === state.selectedIsbtGroup;
      });
      if (!state.isbtTableSearch.trim()) return rows;
      return rows.filter((row) => {
        const joined = Object.values(row)
          .filter((value) => !Array.isArray(value) && typeof value !== "object")
          .map(text)
          .join(" || ");
        return matchesSearchGroups(joined, state.isbtTableSearch);
      });
    }

    function isbtCurrentDetailColumns() {
      const isbtCols = (DATA.isbt && DATA.isbt.columns) ? DATA.isbt.columns : {};
      if (state.granularity === "Allele") {
        return Array.isArray(isbtCols.grouped_view) ? isbtCols.grouped_view : [];
      }
      return Array.isArray(isbtCols.variant_view) ? isbtCols.variant_view : [];
    }

    function isbtNormalizeVariantRow(row, fallbackGroup) {
      const group = text(row.__group || row.Group || fallbackGroup);
      const out = {};
      out.Group = group;
      out.__group = group;
      out.Allele_id = text(row.Allele_id || row.allele_id || "");
      out.ISBT_Allele = text(row.ISBT_Allele || row.isbt_allele || row.allele_name || "");
      out.Gene_name = text(row.Gene_name || row.gene_name || row.gene || "");
      out.Variant_id = text(row.Variant_id || row.variant_id || row.id || "");
      out["DNA Change"] = text(row["DNA Change"] || row.dna_change || row.input || row.nucleotide_change || row.hgvs || "-") || "-";
      out["Exon/Intron"] = text(row["Exon/Intron"] || row.exon_intron || "-") || "-";
      out["HGVS Transcript"] = text(row["HGVS Transcript"] || row.hgvs_transcript || row.input || row.hgvs || "-") || "-";
      ["sv_allele", "null_allele", "mod_allele", "partial_allele", "weak_allele", "el_allele"].forEach((flag) => {
        out[flag] = text(row[flag] || "");
      });
      return Object.assign({}, row, out);
    }

    function isbtVariantRowsFromAlleleRow(alleleRow) {
      const group = text(alleleRow.__group || alleleRow.Group || "");
      const explicitRows = Array.isArray(alleleRow.__variant_view_rows) ? alleleRow.__variant_view_rows : [];
      if (explicitRows.length > 0) {
        return explicitRows.map((row) => isbtNormalizeVariantRow(row, group));
      }
      const rawRows = Array.isArray(alleleRow.__raw_variant_rows) ? alleleRow.__raw_variant_rows : [];
      return rawRows.map((row) => isbtNormalizeVariantRow(row, group));
    }

    function isbtAlleleRowsFromVariantRow(variantRow) {
      const groupedRows = (DATA.isbt && Array.isArray(DATA.isbt.grouped_rows)) ? DATA.isbt.grouped_rows : [];
      const group = text(variantRow.__group || variantRow.Group || "");
      const variantId = text(variantRow.Variant_id || variantRow.variant_id || variantRow.id || "");
      const alleleId = text(variantRow.Allele_id || variantRow.allele_id || "");
      const isbtAllele = text(variantRow.ISBT_Allele || variantRow.isbt_allele || "");
      return groupedRows.filter((row) => {
        if (group && text(row.__group) !== group) return false;
        const rawRows = Array.isArray(row.__raw_variant_rows) ? row.__raw_variant_rows : [];
        if (variantId) {
          const hasVariant = rawRows.some((raw) => {
            return text(raw.variant_id || raw.Variant_id || raw.id || "") === variantId;
          });
          if (hasVariant) return true;
        }
        const rowAlleleId = text(row.Allele_id || "");
        const rowIsbtAllele = text(row.ISBT_Allele || "");
        if (alleleId && rowAlleleId === alleleId) return true;
        if (isbtAllele && rowIsbtAllele === isbtAllele) return true;
        return false;
      });
    }

    function mapRowsByIds(sourceRows, idField, ids) {
      const byId = new Map();
      sourceRows.forEach((row) => byId.set(text(row[idField]), row));
      const out = [];
      ids.forEach((id) => {
        const row = byId.get(text(id));
        if (row) out.push(row);
      });
      return out;
    }

    function setActiveButtons() {
      const alleleBtn = document.getElementById("isbt-view-allele-btn");
      const variantBtn = document.getElementById("isbt-view-variant-btn");
      if (alleleBtn) alleleBtn.classList.toggle("active", state.granularity === "Allele");
      if (variantBtn) variantBtn.classList.toggle("active", state.granularity === "Variant");
    }

    function currentVisibleStateText() {
      return "Mode: ISBT"
        + " | Group: " + (state.selectedIsbtGroup || "")
        + " | View: ISBT_" + (state.granularity === "Allele" ? "Allele_Table" : "Variant_Table")
        + " | GroupSearch: " + state.isbtGroupSearch
        + " | ISBTTableSearch: " + state.isbtTableSearch;
    }

    function updateStateField() {
      document.getElementById("feedback-state").value = currentVisibleStateText();
    }

    function activeTableInfo() {
      if (state.activeTable === "isbtDetail") {
        const label = state.granularity === "Allele" ? "ISBT Allele Table" : "ISBT Variant Table";
        return { label: label, checked: state.checked.isbtDetail };
      }
      if (state.activeTable === "modalVariants") {
        const label = state.modal.type === "isbtVariants" ? "ISBT Variants Popup Table" : "Variants Popup Table";
        return { label: label, checked: state.checked.modalVariants };
      }
      if (state.activeTable === "modalAlleles") {
        const label = state.modal.type === "isbtAlleles" ? "ISBT Alleles Popup Table" : "Alleles Popup Table";
        return { label: label, checked: state.checked.modalAlleles };
      }
      const fallbackLabel = state.granularity === "Allele" ? "ISBT Allele Table" : "ISBT Variant Table";
      return { label: fallbackLabel, checked: state.checked.isbtDetail };
    }

    function updateDownloadContext() {
      const info = activeTableInfo();
      document.getElementById("download-context").textContent =
        "Active table: " + info.label + " | Checked rows: " + info.checked.size;
    }

    function applyLayoutState() {
      const appRoot = document.querySelector(".app-root");
      const panelBtn = document.getElementById("toggle-panel-b");
      const alleleLayout = document.getElementById("allele-layout");
      const alleleBtn = document.getElementById("toggle-gene-selector");
      const isbtLayout = document.getElementById("isbt-layout");
      const isbtBtn = document.getElementById("toggle-isbt-selector");

      if (appRoot) {
        appRoot.classList.toggle("panel-b-collapsed", state.panelBCollapsed);
      }
      if (panelBtn) {
        panelBtn.textContent = state.panelBCollapsed ? "Show Sidebar" : "Hide Sidebar";
      }
      if (alleleLayout) {
        alleleLayout.classList.toggle("gene-selector-collapsed", state.alleleSelectorCollapsed);
      }
      if (alleleBtn) {
        alleleBtn.textContent = state.alleleSelectorCollapsed ? "Show Gene Selector" : "Hide Gene Selector";
      }
      if (isbtLayout) {
        isbtLayout.classList.toggle("group-selector-collapsed", state.isbtSelectorCollapsed);
      }
      if (isbtBtn) {
        isbtBtn.textContent = state.isbtSelectorCollapsed ? "Show Group Selector" : "Hide Group Selector";
      }
    }

    function closeVizModal() {
      document.getElementById("viz-modal").classList.remove("open");
      document.getElementById("viz-content").innerHTML = "";
      state.vizScale = 1.0;
      const label = document.getElementById("viz-scale-label");
      if (label) label.textContent = "100%";
    }

    function openHelpModal() {
      document.getElementById("help-modal").classList.add("open");
    }

    function closeHelpModal() {
      document.getElementById("help-modal").classList.remove("open");
    }

    function applyVizScale() {
      const img = document.getElementById("viz-image");
      const label = document.getElementById("viz-scale-label");
      if (img) {
        img.style.transform = "scale(" + state.vizScale + ")";
      }
      if (label) {
        label.textContent = Math.round(state.vizScale * 100) + "%";
      }
    }

    function setVizScale(nextScale) {
      const clamped = Math.max(0.25, Math.min(8.0, nextScale));
      state.vizScale = clamped;
      applyVizScale();
    }

    function renderSequencePanel(sequencePanel) {
      const panel = document.createElement("div");
      panel.className = "viz-seq-panel";

      const meta = document.createElement("div");
      meta.className = "viz-seq-meta";
      const status = text(sequencePanel && sequencePanel.status ? sequencePanel.status : "unavailable");
      const message = text(sequencePanel && sequencePanel.message ? sequencePanel.message : "Sequence unavailable.");
      const metaLine = text(sequencePanel && sequencePanel.meta ? sequencePanel.meta : "");
      meta.textContent = message + (metaLine ? (" | " + metaLine) : "");
      panel.appendChild(meta);

      if (status !== "ok") {
        return panel;
      }

      const rows = Array.isArray(sequencePanel.rows) ? sequencePanel.rows : [];
      if (rows.length === 0) {
        return panel;
      }

      const wrap = document.createElement("div");
      wrap.className = "viz-seq-table-wrap";
      const table = document.createElement("table");
      table.className = "viz-seq-table";
      const thead = document.createElement("thead");
      thead.innerHTML = "<tr><th>Version</th><th>5'</th><th>Variant</th><th>3'</th></tr>";
      table.appendChild(thead);
      const tbody = document.createElement("tbody");
      rows.forEach((row) => {
        const tr = document.createElement("tr");
        const tdLabel = document.createElement("td");
        tdLabel.textContent = text(row.label || "");
        tr.appendChild(tdLabel);
        const td5 = document.createElement("td");
        td5.textContent = text(row.five_prime || "");
        tr.appendChild(td5);
        const tdCenter = document.createElement("td");
        tdCenter.className = "viz-seq-center";
        tdCenter.textContent = text(row.center || "");
        tr.appendChild(tdCenter);
        const td3 = document.createElement("td");
        td3.textContent = text(row.three_prime || "");
        tr.appendChild(td3);
        tbody.appendChild(tr);
      });
      table.appendChild(tbody);
      wrap.appendChild(table);
      panel.appendChild(wrap);

      const warnings = Array.isArray(sequencePanel.warnings) ? sequencePanel.warnings : [];
      warnings.forEach((warning) => {
        const warn = document.createElement("div");
        warn.className = "viz-seq-warn";
        warn.textContent = "Warning: " + text(warning);
        panel.appendChild(warn);
      });
      return panel;
    }

    async function requestVizForRow(vizSource, row) {
      const modal = document.getElementById("viz-modal");
      const titleNode = document.getElementById("viz-title");
      const content = document.getElementById("viz-content");
      state.vizScale = 1.0;
      titleNode.textContent = "Rendering visualization...";
      content.innerHTML = "<div class='muted'>Generating plot...</div>";
      modal.classList.add("open");

      try {
        const res = await fetch("/api/viz", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ source: vizSource, row: row || {} })
        });
        const data = await res.json();
        if (!data.ok) {
          titleNode.textContent = "Visualization unavailable";
          content.innerHTML = "<div class='viz-error'>" + text(data.error || "Unknown error") + "</div>";
          return;
        }
        titleNode.textContent = text(data.title || "Visualization");
        const notes = Array.isArray(data.notes) ? data.notes : [];
        const sequencePanel = data.sequence_panel || null;
        const img = document.createElement("img");
        img.id = "viz-image";
        img.alt = "DNA feature visualization";
        const svgB64 = text(data.image_svg_base64 || "");
        const pngB64 = text(data.image_png_base64 || "");
        img.src = svgB64
          ? ("data:image/svg+xml;base64," + svgB64)
          : ("data:image/png;base64," + pngB64);
        content.innerHTML = "";
        if (notes.length > 0) {
          const notesDiv = document.createElement("div");
          notesDiv.className = "viz-notes";
          notesDiv.textContent = notes.join(" | ");
          content.appendChild(notesDiv);
        }
        content.appendChild(img);
        content.appendChild(renderSequencePanel(sequencePanel));
        applyVizScale();
      } catch (err) {
        titleNode.textContent = "Visualization unavailable";
        content.innerHTML = "<div class='viz-error'>Failed to generate visualization.</div>";
      }
    }

    function renderTable(containerId, config) {
      const container = document.getElementById(containerId);
      container.innerHTML = "";

      const rows = config.rows || [];
      const rawRows = config.rawRows || rows;
      const columns = config.columns || [];
      const checkedSet = config.checkedSet;
      const maxHeight = config.maxHeight || "52vh";
      const multilineCols = new Set(config.multilineCols || []);
      const truncateTooltipCols = new Set(config.truncateTooltipCols || []);
      const truncateLength = Number.isInteger(config.truncateLength) ? config.truncateLength : 10;
      const wrapCols = new Set(["Allele_name", "Reference_No_PMID", "Accession_number", "Phenotype"]);
      (config.wrapCols || []).forEach((col) => wrapCols.add(col));
      const tightWrapCols = new Set(config.tightWrapCols || []);
      const showTopScrollbar = !!config.topScrollbar;
      const vizSource = text(config.vizSource || "").trim();
      const showViz = !!vizSource;
      const rowKey = typeof config.rowKey === "function"
        ? config.rowKey
        : ((row, idx) => String(idx + 1));

      const validKeys = new Set(rows.map((row, idx) => text(rowKey(row, idx))));
      Array.from(checkedSet).forEach((key) => {
        if (!validKeys.has(text(key))) checkedSet.delete(key);
      });

      const shell = document.createElement("div");
      shell.className = "table-shell";
      shell.style.maxHeight = maxHeight;
      shell.style.height = maxHeight;
      let topScroll = null;
      let topScrollInner = null;
      if (showTopScrollbar) {
        topScroll = document.createElement("div");
        topScroll.className = "table-top-scroll";
        topScrollInner = document.createElement("div");
        topScrollInner.className = "table-top-scroll-inner";
        topScroll.appendChild(topScrollInner);
      }

      if (rows.length === 0) {
        const empty = document.createElement("div");
        empty.className = "muted";
        empty.style.padding = "10px";
        empty.textContent = "No rows to display.";
        shell.appendChild(empty);
        if (topScroll) container.appendChild(topScroll);
        container.appendChild(shell);
        return;
      }

      const table = document.createElement("table");
      table.className = "data-table";

      const thead = document.createElement("thead");
      const htr = document.createElement("tr");
      const hCheck = document.createElement("th");
      hCheck.className = "check-cell";
      htr.appendChild(hCheck);
      if (showViz) {
        const hViz = document.createElement("th");
        hViz.className = "viz-cell";
        hViz.textContent = "Viz";
        htr.appendChild(hViz);
      }
      columns.forEach((col) => {
        const th = document.createElement("th");
        th.textContent = col;
        if (wrapCols.has(col)) th.classList.add("allele-name-col");
        if (tightWrapCols.has(col)) th.classList.add("tight-wrap-col");
        htr.appendChild(th);
      });
      thead.appendChild(htr);
      table.appendChild(thead);

      const tbody = document.createElement("tbody");
      rows.forEach((row, idx) => {
        const key = text(rowKey(row, idx));
        const rawRow = idx < rawRows.length ? rawRows[idx] : row;
        const tr = document.createElement("tr");

        const checkTd = document.createElement("td");
        checkTd.className = "check-cell";
        const cb = document.createElement("input");
        cb.type = "checkbox";
        cb.checked = checkedSet.has(key);
        cb.addEventListener("change", (ev) => {
          ev.stopPropagation();
          if (cb.checked) {
            checkedSet.add(key);
          } else {
            checkedSet.delete(key);
          }
          state.activeTable = config.tableKey;
          updateDownloadContext();
        });
        checkTd.appendChild(cb);
        tr.appendChild(checkTd);

        if (showViz) {
          const vizTd = document.createElement("td");
          vizTd.className = "viz-cell";
          const vizBtn = document.createElement("button");
          vizBtn.type = "button";
          vizBtn.className = "viz-btn";
          vizBtn.textContent = "Viz";
          vizBtn.addEventListener("click", (ev) => {
            ev.stopPropagation();
            requestVizForRow(vizSource, rawRow);
          });
          vizTd.appendChild(vizBtn);
          tr.appendChild(vizTd);
        }

        columns.forEach((col) => {
          const td = document.createElement("td");
          const fullText = text(row[col]);
          if (truncateTooltipCols.has(col) && fullText.length > truncateLength) {
            td.textContent = "...";
            td.title = fullText;
          } else {
            td.textContent = fullText;
          }
          if (multilineCols.has(col)) td.classList.add("multiline");
          if (wrapCols.has(col)) td.classList.add("allele-name-col");
          if (tightWrapCols.has(col)) td.classList.add("tight-wrap-col");
          tr.appendChild(td);
        });

        tr.addEventListener("click", (ev) => {
          if (ev.target && ev.target.tagName === "INPUT") return;
          if (typeof config.onRowClick === "function") config.onRowClick(row, idx, key, rawRow);
        });
        tr.addEventListener("dblclick", (ev) => {
          if (ev.target && ev.target.tagName === "INPUT") return;
          if (typeof config.onRowDblClick === "function") config.onRowDblClick(row, idx, key, rawRow);
        });
        tbody.appendChild(tr);
      });
      table.appendChild(tbody);

      shell.appendChild(table);
      if (topScroll && topScrollInner) {
        let syncing = false;
        const syncTopBarWidth = () => {
          topScrollInner.style.width = table.scrollWidth + "px";
          topScroll.scrollLeft = shell.scrollLeft;
        };
        topScroll.addEventListener("scroll", () => {
          if (syncing) return;
          syncing = true;
          shell.scrollLeft = topScroll.scrollLeft;
          syncing = false;
        });
        shell.addEventListener("scroll", () => {
          if (syncing) return;
          syncing = true;
          topScroll.scrollLeft = shell.scrollLeft;
          syncing = false;
        });
        syncTopBarWidth();
        if (window.ResizeObserver) {
          const observer = new ResizeObserver(syncTopBarWidth);
          observer.observe(shell);
          observer.observe(table);
        }
        container.appendChild(topScroll);
      }
      container.appendChild(shell);
    }

    function renderGeneSelectorTable() {
      const container = document.getElementById("gene-selector-table");
      container.innerHTML = "";

      const genes = getFilteredGenes();
      const shell = document.createElement("div");
      shell.className = "table-shell";
      shell.style.height = "100%";
      shell.style.maxHeight = "100%";

      if (genes.length === 0) {
        const empty = document.createElement("div");
        empty.className = "muted";
        empty.style.padding = "10px";
        empty.textContent = "No genes match the current search.";
        shell.appendChild(empty);
        container.appendChild(shell);
        return;
      }

      const table = document.createElement("table");
      table.className = "data-table";
      const thead = document.createElement("thead");
      const htr = document.createElement("tr");
      const th = document.createElement("th");
      th.textContent = "Gene";
      htr.appendChild(th);
      thead.appendChild(htr);
      table.appendChild(thead);

      const tbody = document.createElement("tbody");
      genes.forEach((gene) => {
        const tr = document.createElement("tr");
        if (gene === state.selectedGene) tr.classList.add("selected-row");
        const td = document.createElement("td");
        td.textContent = gene;
        tr.appendChild(td);
        tr.addEventListener("click", () => {
          if (state.selectedGene !== gene) {
            state.selectedGene = gene;
            state.checked.detail.clear();
            renderAlleleMode();
            updateStateField();
            updateDownloadContext();
          }
        });
        tbody.appendChild(tr);
      });
      table.appendChild(tbody);
      shell.appendChild(table);
      container.appendChild(shell);
    }

    function renderIsbtGroupSelectorTable() {
      const container = document.getElementById("isbt-group-selector-table");
      container.innerHTML = "";

      const groups = getFilteredIsbtGroups();
      const shell = document.createElement("div");
      shell.className = "table-shell";
      shell.style.height = "100%";
      shell.style.maxHeight = "100%";

      if (groups.length === 0) {
        const empty = document.createElement("div");
        empty.className = "muted";
        empty.style.padding = "10px";
        const isbtMeta = (DATA.isbt && DATA.isbt.metadata) ? DATA.isbt.metadata : {};
        const loadError = text(isbtMeta.error || isbtMeta.refresh_error || "");
        empty.textContent = loadError
          ? ("ISBT data unavailable: " + loadError)
          : "No ISBT groups match the current search.";
        shell.appendChild(empty);
        container.appendChild(shell);
        return;
      }

      const table = document.createElement("table");
      table.className = "data-table";
      const thead = document.createElement("thead");
      const htr = document.createElement("tr");
      const th = document.createElement("th");
      th.textContent = "Group";
      htr.appendChild(th);
      thead.appendChild(htr);
      table.appendChild(thead);

      const tbody = document.createElement("tbody");
      groups.forEach((group) => {
        const tr = document.createElement("tr");
        if (group === state.selectedIsbtGroup) tr.classList.add("selected-row");
        const td = document.createElement("td");
        td.textContent = group;
        tr.appendChild(td);
        tr.addEventListener("click", () => {
          if (state.selectedIsbtGroup !== group) {
            state.selectedIsbtGroup = group;
            state.checked.isbtDetail.clear();
            renderIsbtMode();
            updateStateField();
            updateDownloadContext();
          }
        });
        tbody.appendChild(tr);
      });

      table.appendChild(tbody);
      shell.appendChild(table);
      container.appendChild(shell);
    }

    function setActiveMainTableForMode() {
      state.mode = "IT";
      state.activeTable = "isbtDetail";
    }

    function openRelationModal(type, title, rows) {
      state.modal.open = true;
      state.modal.type = type;
      state.modal.rows = rows.slice();
      state.modal.title = title;
      if (type === "variants" || type === "isbtVariants") {
        state.checked.modalVariants.clear();
        state.activeTable = "modalVariants";
      } else if (type === "alleles" || type === "isbtAlleles") {
        state.checked.modalAlleles.clear();
        state.activeTable = "modalAlleles";
      }
      document.getElementById("modal-title").textContent = title;
      document.getElementById("relation-modal").classList.add("open");
      renderModalTable();
      updateDownloadContext();
    }

    function closeRelationModal() {
      state.modal.open = false;
      state.modal.type = null;
      state.modal.rows = [];
      state.modal.title = "";
      document.getElementById("relation-modal").classList.remove("open");
      setActiveMainTableForMode();
      updateDownloadContext();
    }

    function renderModalTable() {
      if (!state.modal.open) {
        document.getElementById("modal-table").innerHTML = "";
        return;
      }
      if (state.modal.type === "variants") {
        const rawRows = state.modal.rows.slice();
        const displayRows = rawRows.map((r) => normalizeRow(r, DATA.columns.variant_view));
        renderTable("modal-table", {
          rows: displayRows,
          rawRows: rawRows,
          columns: DATA.columns.variant_view,
          rowKey: (row, idx) => String(idx + 1),
          checkedSet: state.checked.modalVariants,
          tableKey: "modalVariants",
          truncateTooltipCols: ["Ref_allele_curated", "Alt_allele_curated"],
          truncateLength: 10,
          vizSource: "variant",
          maxHeight: "45vh"
        });
      } else if (state.modal.type === "isbtVariants") {
        const columns = (DATA.isbt && DATA.isbt.columns && Array.isArray(DATA.isbt.columns.variant_view))
          ? DATA.isbt.columns.variant_view
          : [];
        const rawRows = state.modal.rows.slice();
        const displayRows = rawRows.map((r) => normalizeRow(r, columns));
        renderTable("modal-table", {
          rows: displayRows,
          rawRows: rawRows,
          columns: columns,
          rowKey: (row, idx) => String(idx + 1),
          checkedSet: state.checked.modalVariants,
          tableKey: "modalVariants",
          vizSource: "isbt",
          maxHeight: "45vh"
        });
      } else if (state.modal.type === "alleles") {
        const displayRows = state.modal.rows.map((r) => normalizeRow(r, DATA.columns.allele_view));
        renderTable("modal-table", {
          rows: displayRows,
          columns: DATA.columns.allele_view,
          rowKey: (row, idx) => String(idx + 1),
          checkedSet: state.checked.modalAlleles,
          tableKey: "modalAlleles",
          multilineCols: ["Variants"],
          maxHeight: "45vh"
        });
      } else if (state.modal.type === "isbtAlleles") {
        const columns = (DATA.isbt && DATA.isbt.columns && Array.isArray(DATA.isbt.columns.grouped_view))
          ? DATA.isbt.columns.grouped_view
          : [];
        const displayRows = state.modal.rows.map((r) => normalizeRow(r, columns));
        renderTable("modal-table", {
          rows: displayRows,
          rawRows: state.modal.rows,
          columns: columns,
          rowKey: (row, idx) => String(idx + 1),
          checkedSet: state.checked.modalAlleles,
          tableKey: "modalAlleles",
          multilineCols: ["DNA Change", "Exon/Intron", "HGVS Transcript"],
          vizSource: "isbt",
          maxHeight: "45vh"
        });
      }
    }

    function renderGeneMode() {
      const displayRows = DATA.gene_table.map((r) => normalizeRow(r, DATA.columns.gene_table));
      renderTable("gene-mode", {
        rows: displayRows,
        columns: DATA.columns.gene_table,
        rowKey: (row, idx) => String(idx + 1),
        checkedSet: state.checked.gene,
        tableKey: "gene",
        maxHeight: "100%"
      });
    }

    function renderExonMode() {
      const displayRows = DATA.exon_table.map((r) => normalizeRow(r, DATA.columns.exon_view));
      renderTable("exon-mode", {
        rows: displayRows,
        rawRows: DATA.exon_table,
        columns: DATA.columns.exon_view,
        rowKey: (row, idx) => String(idx + 1),
        checkedSet: state.checked.exon,
        tableKey: "exon",
        vizSource: "exon",
        maxHeight: "100%"
      });
    }

    function renderAlleleMode() {
      renderGeneSelectorTable();

      const detailRows = detailBaseRows();
      if (state.granularity === "Allele") {
        const displayRows = detailRows.map((r) => normalizeRow(r, DATA.columns.allele_view));
        renderTable("detail-table", {
          rows: displayRows,
          columns: DATA.columns.allele_view,
          rowKey: (row, idx) => String(idx + 1),
          checkedSet: state.checked.detail,
          tableKey: "detail",
          multilineCols: ["Variants"],
          maxHeight: "100%",
          onRowDblClick: (row) => {
            const alleleId = text(row.Allele_id);
            const linked = uniqueOrdered(DATA.bridge_by_allele[alleleId] || []);
            const rows = mapRowsByIds(DATA.variant_table, "Variant_id", linked);
            const label = text(row.Allele_name) || alleleId;
            openRelationModal("variants", "Variants associated with allele " + label, rows);
          }
        });
      } else {
        const displayRows = detailRows.map((r) => normalizeRow(r, DATA.columns.variant_view));
        renderTable("detail-table", {
          rows: displayRows,
          rawRows: detailRows,
          columns: DATA.columns.variant_view,
          rowKey: (row, idx) => String(idx + 1),
          checkedSet: state.checked.detail,
          tableKey: "detail",
          truncateTooltipCols: ["Ref_allele_curated", "Alt_allele_curated"],
          truncateLength: 10,
          vizSource: "variant",
          maxHeight: "100%",
          onRowClick: (row) => {
            const variantId = text(row.Variant_id);
            const linked = uniqueOrdered(DATA.bridge_by_variant[variantId] || []);
            const rows = mapRowsByIds(DATA.allele_table, "Allele_id", linked);
            const note = text(row.Nucleotide_change);
            const label = note ? (variantId + " (" + note + ")") : variantId;
            openRelationModal("alleles", "Alleles associated with variant " + label, rows);
          }
        });
      }
    }

    function renderIsbtMode() {
      ensureSelectedIsbtGroup();
      setActiveButtons();
      renderIsbtGroupSelectorTable();

      const lastUpdated = document.getElementById("isbt-last-updated");
      const tableSearchInput = document.getElementById("isbt-table-search");
      const isbtMeta = (DATA.isbt && DATA.isbt.metadata) ? DATA.isbt.metadata : {};
      const updatedText = text(isbtMeta.updated_at || "Not available");
      const groupedCount = Number((isbtMeta.row_counts || {}).grouped_rows || 0).toLocaleString();
      const variantCount = Number((isbtMeta.row_counts || {}).variant_rows || 0).toLocaleString();
      const metaNotes = [];
      if (isbtMeta.refresh_error) metaNotes.push("refresh error");
      if (isbtMeta.cache_write_error) metaNotes.push("cache not persisted");
      if (lastUpdated) {
        const noteSuffix = metaNotes.length > 0 ? (" | Notes: " + metaNotes.join(", ")) : "";
        lastUpdated.textContent = "Last pulled: " + updatedText + " | Allele rows: " + groupedCount + " | Variant rows: " + variantCount + noteSuffix;
      }
      if (tableSearchInput) {
        tableSearchInput.placeholder = state.granularity === "Allele"
          ? "Filter grouped allele rows"
          : "Filter variant rows";
      }

      const detailRows = isbtDetailBaseRows();
      const columns = isbtCurrentDetailColumns();
      const displayRows = detailRows.map((row) => normalizeRow(row, columns));

      const tableConfig = {
        rows: displayRows,
        rawRows: detailRows,
        columns: columns,
        rowKey: (row, idx) => String(idx + 1),
        checkedSet: state.checked.isbtDetail,
        tableKey: "isbtDetail",
        vizSource: "isbt",
        maxHeight: "100%"
      };

      if (state.granularity === "Allele") {
        tableConfig.multilineCols = ["DNA Change", "Exon/Intron", "HGVS Transcript"];
        tableConfig.topScrollbar = true;
        const openAlleleVariantPopup = (row, idx, key, rawRow) => {
          const source = rawRow || row;
          const rows = isbtVariantRowsFromAlleleRow(source);
          const label = text(source.ISBT_Allele || source.Allele_id || source.__row_key || "");
          openRelationModal("isbtVariants", "Variants associated with allele " + label, rows);
        };
        tableConfig.onRowClick = openAlleleVariantPopup;
        tableConfig.onRowDblClick = openAlleleVariantPopup;
      } else {
        tableConfig.wrapCols = ["Allele_id", "ISBT_Allele"];
        tableConfig.tightWrapCols = ["Allele_id", "ISBT_Allele"];
        tableConfig.topScrollbar = true;
        const openVariantAllelePopup = (row, idx, key, rawRow) => {
          const source = rawRow || row;
          const rows = isbtAlleleRowsFromVariantRow(source);
          const variantId = text(source.Variant_id || source.variant_id || "");
          const dna = text(source["DNA Change"] || source.dna_change || source.input || "");
          const label = dna ? (variantId + " (" + dna + ")") : variantId;
          openRelationModal("isbtAlleles", "Alleles associated with variant " + label, rows);
        };
        tableConfig.onRowClick = openVariantAllelePopup;
        tableConfig.onRowDblClick = openVariantAllelePopup;
      }

      renderTable("isbt-detail-table", tableConfig);
    }

    function refreshMainPanels() {
      applyLayoutState();
      setActiveButtons();
      setActiveMainTableForMode();
      renderIsbtMode();
      updateStateField();
      updateDownloadContext();
    }

    function clearAllCheckedRows() {
      state.checked.gene.clear();
      state.checked.exon.clear();
      state.checked.detail.clear();
      state.checked.isbtDetail.clear();
      state.checked.modalVariants.clear();
      state.checked.modalAlleles.clear();
      refreshMainPanels();
      if (state.modal.open) renderModalTable();
    }

    function toCsv(rows, columns) {
      const cols = columns.slice();
      const esc = (val) => {
        const s = text(val).replace(/"/g, '""');
        if (s.includes(",") || s.includes("\\n") || s.includes('"')) {
          return '"' + s + '"';
        }
        return s;
      };
      const lines = [];
      lines.push(cols.map(esc).join(","));
      rows.forEach((row) => {
        lines.push(cols.map((c) => esc(row[c])).join(","));
      });
      return lines.join("\\r\\n");
    }

    function saveCsv(rows, columns, filename) {
      const csvText = toCsv(rows, columns);
      const blob = new Blob([csvText], { type: "text/csv;charset=utf-8;" });
      const url = URL.createObjectURL(blob);
      const link = document.createElement("a");
      link.href = url;
      link.download = filename;
      document.body.appendChild(link);
      link.click();
      link.remove();
      URL.revokeObjectURL(url);
    }

    function checkedRowsByIndex(sourceRows, checkedSet) {
      const idx = Array.from(checkedSet)
        .map((v) => parseInt(v, 10))
        .filter((n) => Number.isInteger(n) && n >= 1 && n <= sourceRows.length)
        .sort((a, b) => a - b);
      return idx.map((n) => sourceRows[n - 1]);
    }

    function expandIsbtGroupedRows(groupedRows) {
      const rawRows = [];
      groupedRows.forEach((row) => {
        const variants = Array.isArray(row.__raw_variant_rows) ? row.__raw_variant_rows : [];
        variants.forEach((raw) => rawRows.push(raw));
      });
      return rawRows;
    }

    function downloadMainSelection() {
      let rows = [];
      let columns = [];
      let filename = "checked_rows.csv";
      const dateTag = new Date().toISOString().slice(0, 10).replace(/-/g, "");

      if (state.activeTable === "gene") {
        columns = DATA.columns.gene_table;
        rows = checkedRowsByIndex(DATA.gene_table, state.checked.gene).map((r) => normalizeRow(r, columns));
        filename = "gene_table_checked_rows_" + dateTag + ".csv";
      } else if (state.activeTable === "exon") {
        columns = DATA.columns.exon_view;
        rows = checkedRowsByIndex(DATA.exon_table, state.checked.exon).map((r) => normalizeRow(r, columns));
        filename = "exon_table_checked_rows_" + dateTag + ".csv";
      } else if (state.activeTable === "detail") {
        if (state.granularity === "Allele") {
          columns = DATA.columns.allele_view;
          rows = checkedRowsByIndex(detailBaseRows(), state.checked.detail).map((r) => normalizeRow(r, columns));
          filename = "allele_detail_table_checked_rows_" + dateTag + ".csv";
        } else {
          columns = DATA.columns.variant_view;
          rows = checkedRowsByIndex(detailBaseRows(), state.checked.detail).map((r) => normalizeRow(r, columns));
          filename = "variant_detail_table_checked_rows_" + dateTag + ".csv";
        }
      } else if (state.activeTable === "modalVariants") {
        columns = state.modal.type === "isbtVariants"
          ? ((DATA.isbt && DATA.isbt.columns && Array.isArray(DATA.isbt.columns.variant_view)) ? DATA.isbt.columns.variant_view : [])
          : DATA.columns.variant_view;
        rows = checkedRowsByIndex(state.modal.rows, state.checked.modalVariants).map((r) => normalizeRow(r, columns));
        filename = "variants_popup_checked_rows_" + dateTag + ".csv";
      } else if (state.activeTable === "modalAlleles") {
        columns = state.modal.type === "isbtAlleles"
          ? ((DATA.isbt && DATA.isbt.columns && Array.isArray(DATA.isbt.columns.grouped_view)) ? DATA.isbt.columns.grouped_view : [])
          : DATA.columns.allele_view;
        rows = checkedRowsByIndex(state.modal.rows, state.checked.modalAlleles).map((r) => normalizeRow(r, columns));
        filename = "alleles_popup_checked_rows_" + dateTag + ".csv";
      } else if (state.activeTable === "isbtDetail") {
        columns = (DATA.isbt && DATA.isbt.columns && Array.isArray(DATA.isbt.columns.raw_export))
          ? DATA.isbt.columns.raw_export
          : [];
        const detailRows = checkedRowsByIndex(isbtDetailBaseRows(), state.checked.isbtDetail);
        if (state.granularity === "Allele") {
          const rawRows = expandIsbtGroupedRows(detailRows);
          rows = rawRows.map((r) => normalizeRow(r, columns));
          filename = "isbt_grouped_selection_expanded_" + dateTag + ".csv";
        } else {
          rows = detailRows.map((r) => normalizeRow(r, columns));
          filename = "isbt_variant_selection_raw_" + dateTag + ".csv";
        }
      }

      saveCsv(rows, columns, filename);
    }

    function downloadModalSelection() {
      const dateTag = new Date().toISOString().slice(0, 10).replace(/-/g, "");
      if (state.modal.type === "variants") {
        const rows = checkedRowsByIndex(state.modal.rows, state.checked.modalVariants)
          .map((r) => normalizeRow(r, DATA.columns.variant_view));
        saveCsv(rows, DATA.columns.variant_view, "variants_popup_checked_rows_" + dateTag + ".csv");
      } else if (state.modal.type === "isbtVariants") {
        const columns = (DATA.isbt && DATA.isbt.columns && Array.isArray(DATA.isbt.columns.variant_view))
          ? DATA.isbt.columns.variant_view
          : [];
        const rows = checkedRowsByIndex(state.modal.rows, state.checked.modalVariants)
          .map((r) => normalizeRow(r, columns));
        saveCsv(rows, columns, "isbt_variants_popup_checked_rows_" + dateTag + ".csv");
      } else if (state.modal.type === "alleles") {
        const rows = checkedRowsByIndex(state.modal.rows, state.checked.modalAlleles)
          .map((r) => normalizeRow(r, DATA.columns.allele_view));
        saveCsv(rows, DATA.columns.allele_view, "alleles_popup_checked_rows_" + dateTag + ".csv");
      } else if (state.modal.type === "isbtAlleles") {
        const columns = (DATA.isbt && DATA.isbt.columns && Array.isArray(DATA.isbt.columns.grouped_view))
          ? DATA.isbt.columns.grouped_view
          : [];
        const rows = checkedRowsByIndex(state.modal.rows, state.checked.modalAlleles)
          .map((r) => normalizeRow(r, columns));
        saveCsv(rows, columns, "isbt_alleles_popup_checked_rows_" + dateTag + ".csv");
      }
    }

    function clearModalSelection() {
      if (state.modal.type === "variants" || state.modal.type === "isbtVariants") {
        state.checked.modalVariants.clear();
      } else if (state.modal.type === "alleles" || state.modal.type === "isbtAlleles") {
        state.checked.modalAlleles.clear();
      }
      renderModalTable();
      updateDownloadContext();
    }

    function formatUtcDateTime(value) {
      const raw = text(value).trim();
      if (!raw) return "Unknown";
      const dt = new Date(raw);
      if (Number.isNaN(dt.getTime())) return raw;
      return dt.toLocaleString(undefined, {
        year: "numeric",
        month: "short",
        day: "2-digit",
        hour: "2-digit",
        minute: "2-digit",
        hour12: false,
        timeZoneName: "short"
      });
    }

    function renderLatestUpdates() {
      const root = document.getElementById("latest-updates");
      root.innerHTML = "";
      const ul = document.createElement("ul");
      const isbtMeta = (DATA.isbt && DATA.isbt.metadata) ? DATA.isbt.metadata : {};
      const tableUpdates = (DATA.isbt && Array.isArray(DATA.isbt.table_updates)) ? DATA.isbt.table_updates : [];

      const lines = [];
      if (tableUpdates.length > 0) {
        lines.push("ISBT table revision history (latest first):");
        tableUpdates.slice(0, 12).forEach((item) => {
          const eventAt = formatUtcDateTime(item.eventAt || item.updatedAt || "");
          const group = text(item.Group || "");
          const operation = text(item.operation || "update");
          const allele = text(item.ISBT_Allele || item.Allele_id || "?");
          const variantId = text(item.Variant_id || "?");
          const dna = text(item["DNA Change"] || "-");
          const version = text(item.version || "");
          const actor = text(item.actor || "");

          let line = eventAt
            + " | " + group
            + " | " + operation
            + " | " + allele
            + " | " + variantId
            + " | " + dna;
          if (version) line += " | v" + version;
          if (actor) line += " | by " + actor;
          lines.push(line);
        });
      } else {
        lines.push("No ISBT table revision updates are currently available.");
      }

      lines.push("ISBT data pulled: " + text(isbtMeta.updated_at || "Not available"));
      lines.push("Cache contingency: cached ISBT data is reused on refresh failure and kept up to 30 days.");

      const pullErrors = Array.isArray(isbtMeta.pull_errors) ? isbtMeta.pull_errors : [];
      if (pullErrors.length > 0) {
        lines.push("ISBT refresh warnings: " + pullErrors.length.toLocaleString() + " fetch errors.");
      }
      if (isbtMeta.refresh_error) {
        lines.push("ISBT refresh error: " + text(isbtMeta.refresh_error));
      }
      if (isbtMeta.error) {
        lines.push("ISBT load error: " + text(isbtMeta.error));
      }
      if (isbtMeta.cache_write_error) {
        lines.push("ISBT cache write warning: " + text(isbtMeta.cache_write_error));
      }

      lines.forEach((line) => {
        const li = document.createElement("li");
        li.textContent = line;
        ul.appendChild(li);
      });
      root.appendChild(ul);
    }

    function renderFeedbackStatus(entry) {
      const box = document.getElementById("feedback-status");
      if (!entry) {
        box.className = "last-feedback muted";
        box.textContent = "No submissions yet.";
        return;
      }
      box.className = "last-feedback";
      box.innerHTML = "";
      const title = document.createElement("strong");
      title.textContent = "Last submission";
      box.appendChild(title);
      const rows = [
        entry.submitted_at,
        "Email: " + (entry.email || "(blank)"),
        "State: " + (entry.state || ""),
        "Comment length: " + text(entry.comment || "").length
      ];
      rows.forEach((line) => {
        const div = document.createElement("div");
        div.textContent = line;
        box.appendChild(div);
      });
    }

    async function submitFeedback() {
      const payload = {
        state: document.getElementById("feedback-state").value,
        email: document.getElementById("feedback-email").value,
        comment: document.getElementById("feedback-comment").value
      };
      const res = await fetch("/api/feedback", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(payload)
      });
      const data = await res.json();
      if (data.ok) {
        document.getElementById("feedback-email").value = "";
        document.getElementById("feedback-comment").value = "";
        document.getElementById("feedback-flash").textContent = "Feedback submitted.";
        renderFeedbackStatus(data.last);
        setTimeout(() => {
          document.getElementById("feedback-flash").textContent = "";
        }, 2200);
      }
    }

    function bindEvents() {
      const bind = (id, eventName, handler) => {
        const node = document.getElementById(id);
        if (node) {
          node.addEventListener(eventName, handler);
        }
      };

      bind("isbt-group-search", "input", (ev) => {
        state.isbtGroupSearch = ev.target.value || "";
        ensureSelectedIsbtGroup();
        state.checked.isbtDetail.clear();
        renderIsbtMode();
        updateStateField();
        updateDownloadContext();
      });
      bind("isbt-table-search", "input", (ev) => {
        state.isbtTableSearch = ev.target.value || "";
        state.checked.isbtDetail.clear();
        renderIsbtMode();
        updateStateField();
        updateDownloadContext();
      });
      bind("isbt-view-allele-btn", "click", () => {
        if (state.granularity !== "Allele") {
          if (state.modal.open) closeRelationModal();
          state.granularity = "Allele";
          state.checked.isbtDetail.clear();
          renderIsbtMode();
          updateStateField();
          updateDownloadContext();
        }
      });
      bind("isbt-view-variant-btn", "click", () => {
        if (state.granularity !== "Variant") {
          if (state.modal.open) closeRelationModal();
          state.granularity = "Variant";
          state.checked.isbtDetail.clear();
          renderIsbtMode();
          updateStateField();
          updateDownloadContext();
        }
      });
      bind("nav-help-btn", "click", openHelpModal);
      bind("toggle-panel-b", "click", () => {
        state.panelBCollapsed = !state.panelBCollapsed;
        applyLayoutState();
      });
      bind("toggle-isbt-selector", "click", () => {
        state.isbtSelectorCollapsed = !state.isbtSelectorCollapsed;
        applyLayoutState();
      });
      bind("download-main-btn", "click", downloadMainSelection);
      bind("clear-main-btn", "click", clearAllCheckedRows);
      bind("feedback-submit", "click", submitFeedback);

      bind("modal-close-btn", "click", closeRelationModal);
      bind("modal-clear-btn", "click", clearModalSelection);
      bind("modal-download-btn", "click", downloadModalSelection);
      bind("relation-modal", "click", (ev) => {
        if (ev.target.id === "relation-modal") closeRelationModal();
      });
      bind("viz-close-btn", "click", closeVizModal);
      bind("viz-zoom-in-btn", "click", () => {
        setVizScale(state.vizScale * 1.25);
      });
      bind("viz-zoom-out-btn", "click", () => {
        setVizScale(state.vizScale / 1.25);
      });
      bind("viz-zoom-reset-btn", "click", () => {
        setVizScale(1.0);
      });
      bind("viz-modal", "click", (ev) => {
        if (ev.target.id === "viz-modal") closeVizModal();
      });
      bind("help-close-btn", "click", closeHelpModal);
      bind("help-modal", "click", (ev) => {
        if (ev.target.id === "help-modal") closeHelpModal();
      });
    }

    async function init() {
      const response = await fetch("/api/init");
      DATA = await response.json();
      state.mode = "IT";
      state.activeTable = "isbtDetail";
      state.selectedIsbtGroup = (DATA.isbt && Array.isArray(DATA.isbt.groups) && DATA.isbt.groups.length > 0)
        ? DATA.isbt.groups[0]
        : "";
      bindEvents();
      renderLatestUpdates();
      renderFeedbackStatus(null);
      refreshMainPanels();
    }

    init();
  </script>
</body>
</html>
"""


class AppHandler(BaseHTTPRequestHandler):
    def _send_json(self, payload, status=200):
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_html(self, html_text, status=200):
        body = html_text.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_json_head(self, payload, status=200):
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()

    def _send_html_head(self, html_text, status=200):
        body = html_text.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()

    def do_GET(self):
        parsed = urlparse(self.path)
        if parsed.path in ("/", "/index.html"):
            self._send_html(INDEX_HTML)
            return
        if parsed.path == "/api/init":
            self._send_json(INIT_PAYLOAD)
            return
        if parsed.path == "/api/feedback/latest":
            latest = FEEDBACK_LOG[-1] if FEEDBACK_LOG else None
            self._send_json({"ok": True, "last": latest, "count": len(FEEDBACK_LOG)})
            return
        self._send_json({"ok": False, "error": "Not found"}, status=404)

    def do_HEAD(self):
        parsed = urlparse(self.path)
        if parsed.path in ("/", "/index.html"):
            self._send_html_head(INDEX_HTML)
            return
        if parsed.path == "/api/init":
            self._send_json_head(INIT_PAYLOAD)
            return
        if parsed.path == "/api/feedback/latest":
            latest = FEEDBACK_LOG[-1] if FEEDBACK_LOG else None
            self._send_json_head({"ok": True, "last": latest, "count": len(FEEDBACK_LOG)})
            return
        self._send_json_head({"ok": False, "error": "Not found"}, status=404)

    def do_POST(self):
        parsed = urlparse(self.path)
        length = int(self.headers.get("Content-Length", "0"))
        raw = self.rfile.read(length) if length > 0 else b"{}"
        try:
            payload = json.loads(raw.decode("utf-8"))
        except Exception:
            self._send_json({"ok": False, "error": "Invalid JSON"}, status=400)
            return

        if parsed.path == "/api/viz":
            source = clean_text(payload.get("source", ""))
            row = payload.get("row", {})
            result = build_row_viz_payload(source, row)
            status = 200 if result.get("ok") else 400
            self._send_json(result, status=status)
            return

        if parsed.path != "/api/feedback":
            self._send_json({"ok": False, "error": "Not found"}, status=404)
            return

        entry = {
            "submitted_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "state": clean_text(payload.get("state", "")),
            "email": clean_text(payload.get("email", "")),
            "comment": clean_text(payload.get("comment", "")),
        }
        FEEDBACK_LOG.append(entry)
        self._send_json({"ok": True, "last": entry, "count": len(FEEDBACK_LOG)})

    def log_message(self, format, *args):
        return


def main():
    host = "0.0.0.0"
    port = int(os.environ.get("PORT", "10000"))
    server = ThreadingHTTPServer((host, port), AppHandler)
    print(f"Serving app at http://{host}:{port}")
    server.serve_forever()


if __name__ == "__main__":
    main()

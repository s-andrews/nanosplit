import csv
import gzip
import io
import json
import os
import random
import string
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from flask import Flask, jsonify, redirect, render_template, request, send_file, url_for

app = Flask(__name__)

BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "data"
DATA_DIR.mkdir(exist_ok=True)

executor = ProcessPoolExecutor(max_workers=2)


def random_job_id(n: int = 10) -> str:
    chars = string.ascii_lowercase + string.digits
    return "".join(random.choice(chars) for _ in range(n))


def rc(seq: str) -> str:
    table = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(table)[::-1]


def open_fastq(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def iter_fastq(path: Path) -> Iterable[Tuple[str, str, str, str]]:
    with open_fastq(path) as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline().rstrip("\n")
            p = f.readline()
            q = f.readline().rstrip("\n")
            if not q:
                break
            yield h.rstrip("\n"), s, p.rstrip("\n"), q


def save_job(job_dir: Path, config: Dict):
    (job_dir / "job.json").write_text(json.dumps(config, indent=2))


def load_job(job_dir: Path) -> Dict:
    p = job_dir / "job.json"
    if not p.exists():
        return {}
    return json.loads(p.read_text())


def output_base_name(filename: str) -> str:
    path = Path(filename or "input.fastq")
    name = path.name
    for suffix in (".gz", ".fastq", ".fq"):
        if name.lower().endswith(suffix):
            name = name[:-len(suffix)]
    return safe_filename_part(name) or "input"


def safe_filename_part(value: str) -> str:
    cleaned = []
    for ch in value.strip():
        if ch.isalnum() or ch in ("-", "_", "."):
            cleaned.append(ch)
        elif ch.isspace():
            cleaned.append("_")
    return "".join(cleaned).strip("._") or "unnamed"


def get_output_base(cfg: Dict) -> str:
    return cfg.get("output_base") or output_base_name(cfg.get("fastq_name", "input.fastq"))


def filtered_fastq_path(job_dir: Path, cfg: Dict) -> Path:
    filename = cfg.get("filtered_name") or f"{get_output_base(cfg)}_filtered.fastq"
    path = job_dir / filename
    legacy_path = job_dir / "filtered.fastq"
    if path.exists() or not legacy_path.exists():
        return path
    return legacy_path


def barcode_label_map(cfg: Dict) -> Dict[str, str]:
    labels = {}
    for row in cfg.get("named_barcodes", []):
        if len(row) < 2:
            continue
        name, barcode = row[0], row[1]
        if not name or not barcode:
            continue
        label = f"{safe_filename_part(name)}_{safe_filename_part(barcode)}"
        labels[barcode] = label
        labels[rc(barcode)] = label
        labels[min(barcode, rc(barcode))] = label
    return labels


def output_label_for_key(key: str, cfg: Dict) -> str:
    labels = barcode_label_map(cfg)
    parts = key.split("|")
    return "__".join(labels.get(part, safe_filename_part(part)) for part in parts)


def add_start_barcode(grouped: Dict, lookup: Dict, barcode: str, count: int = 1):
    match = lookup.get(barcode)
    if match is None:
        key = barcode
        grouped[key] = {"forward": 0, "reverse": 0}
        lookup[barcode] = (key, "forward")
        reverse = rc(barcode)
        if reverse != barcode:
            lookup[reverse] = (key, "reverse")
        match = lookup[barcode]

    key, orientation = match
    grouped[key][orientation] += count


def unique_filename(filename: str, used: set) -> str:
    path = Path(filename)
    candidate = filename
    i = 2
    while candidate in used:
        candidate = f"{path.stem}_{i}{path.suffix}"
        i += 1
    used.add(candidate)
    return candidate


def update_status(job_dir: Path, step: str, status: str):
    st = {}
    sp = job_dir / "status.json"
    if sp.exists():
        st = json.loads(sp.read_text())
    st[step] = status
    sp.write_text(json.dumps(st, indent=2))


def analyze_lengths(job_dir_s: str):
    job_dir = Path(job_dir_s)
    update_status(job_dir, "lengths", "running")
    config = load_job(job_dir)
    fastq_path = job_dir / config["fastq_name"]
    lengths = [len(seq) for _, seq, _, _ in iter_fastq(fastq_path)]
    (job_dir / "lengths.json").write_text(json.dumps(lengths))
    update_status(job_dir, "lengths", "done")


def filter_and_barcodes(job_dir_s: str, min_len: int, max_len: int):
    job_dir = Path(job_dir_s)
    update_status(job_dir, "filter_barcodes", "running")
    cfg = load_job(job_dir)
    in_path = job_dir / cfg["fastq_name"]
    cfg["filtered_name"] = f"{get_output_base(cfg)}_filtered.fastq"
    save_job(job_dir, cfg)
    out_path = filtered_fastq_path(job_dir, cfg)
    bc_len = int(cfg["barcode_length"])
    grouped = {}
    barcode_lookup = {}
    total = 0
    with open(out_path, "wt") as out:
        for h, s, p, q in iter_fastq(in_path):
            l = len(s)
            if min_len <= l <= max_len:
                total += 1
                out.write(f"{h}\n{s}\n{p}\n{q}\n")
                add_start_barcode(grouped, barcode_lookup, s[:bc_len])
    barcodes = []
    for k, v in grouped.items():
        count = v["forward"] + v["reverse"]
        barcodes.append({
            "barcode": k,
            "forward": v["forward"],
            "reverse": v["reverse"],
            "count": count,
            "pct": (100.0 * count / total) if total else 0,
        })
    barcodes.sort(key=lambda x: x["count"], reverse=True)
    (job_dir / "barcode_step1.json").write_text(json.dumps({"total": total, "barcodes": barcodes}, indent=2))
    update_status(job_dir, "filter_barcodes", "done")


def end_match(seq: str, barcode: str) -> Optional[str]:
    candidates = [barcode, rc(barcode)]
    for offset in (0, 1, 2, 3):
        if len(seq) < len(barcode) + offset:
            continue
        frag = seq[offset:offset + len(barcode)]
        if frag in candidates:
            return frag
        frag2 = seq[len(seq) - len(barcode) - offset: len(seq) - offset if offset else len(seq)]
        if frag2 in candidates:
            return frag2
    return None


def quantitate_final(job_dir_s: str, selected: List[str]):
    job_dir = Path(job_dir_s)
    update_status(job_dir, "final_quant", "running")
    cfg = load_job(job_dir)
    dual = cfg["barcode_mode"] == "both"
    counts = Counter()
    path = filtered_fastq_path(job_dir, cfg)
    for _, s, _, _ in iter_fastq(path):
        left = None
        right = None
        for b in selected:
            if left is None and end_match(s, b):
                left = b
            if right is None and end_match(s[::-1], b):
                right = b
        if dual:
            if left and right:
                counts[f"{left}|{right}"] += 1
        else:
            if left:
                counts[left] += 1
    total = sum(counts.values())
    rows = [{"key": k, "count": v, "pct": (100.0 * v / total) if total else 0} for k, v in counts.items()]
    rows.sort(key=lambda x: x["count"], reverse=True)
    (job_dir / "final_quant.json").write_text(json.dumps({"total": total, "rows": rows}, indent=2))
    update_status(job_dir, "final_quant", "done")


def split_fastq(job_dir_s: str, selected_keys: List[str]):
    job_dir = Path(job_dir_s)
    update_status(job_dir, "split", "running")
    cfg = load_job(job_dir)
    dual = cfg["barcode_mode"] == "both"
    chosen = set(selected_keys)
    output_base = get_output_base(cfg)
    outs = {}
    split_files = []
    used_filenames = set()
    for k in chosen:
        filename = unique_filename(f"{output_base}_{output_label_for_key(k, cfg)}.fastq", used_filenames)
        split_files.append(filename)
        outs[k] = open(job_dir / filename, "wt")
    unmatched_name = unique_filename(f"{output_base}_unmatched.fastq", used_filenames)
    split_files.append(unmatched_name)
    unmatched = open(job_dir / unmatched_name, "wt")

    selected_barcodes = cfg.get("selected_step4", [])
    for h, s, p, q in iter_fastq(filtered_fastq_path(job_dir, cfg)):
        left = None
        right = None
        for b in selected_barcodes:
            if left is None and end_match(s, b):
                left = b
            if right is None and end_match(s[::-1], b):
                right = b
        key = f"{left}|{right}" if dual and left and right else left
        target = outs.get(key)
        if target:
            target.write(f"{h}\n{s}\n{p}\n{q}\n")
        else:
            unmatched.write(f"{h}\n{s}\n{p}\n{q}\n")

    for o in outs.values():
        o.close()
    unmatched.close()
    cfg["split_files"] = sorted(split_files)
    save_job(job_dir, cfg)
    update_status(job_dir, "split", "done")


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/create_job", methods=["POST"])
def create_job():
    fq = request.files["fastq_file"]
    barcode_len = int(request.form["barcode_length"])
    barcode_mode = request.form.get("barcode_mode", "one")
    named = request.form.get("named_barcodes", "")

    while True:
        job_id = random_job_id()
        job_dir = DATA_DIR / job_id
        if not job_dir.exists():
            break
    job_dir.mkdir(parents=True, exist_ok=False)

    filename = fq.filename or "input.fastq"
    saved_name = Path(filename).name
    fq.save(job_dir / saved_name)

    parsed_named = []
    for line in named.splitlines():
        if not line.strip():
            continue
        parts = [x.strip() for x in line.split(",")]
        parsed_named.append(parts)

    cfg = {
        "job_id": job_id,
        "fastq_name": saved_name,
        "output_base": output_base_name(saved_name),
        "barcode_length": barcode_len,
        "barcode_mode": barcode_mode,
        "named_barcodes": parsed_named,
    }
    save_job(job_dir, cfg)
    update_status(job_dir, "lengths", "queued")
    executor.submit(analyze_lengths, str(job_dir))
    return redirect(url_for("job_page", job_id=job_id))


@app.route("/job/<job_id>")
def job_page(job_id):
    job_dir = DATA_DIR / job_id
    if not job_dir.exists():
        return "Job not found", 404
    cfg = load_job(job_dir)
    return render_template("job.html", job=cfg)


@app.route("/api/job/<job_id>/status")
def job_status(job_id):
    job_dir = DATA_DIR / job_id
    sp = job_dir / "status.json"
    if not sp.exists():
        return jsonify({})
    return jsonify(json.loads(sp.read_text()))


@app.route("/api/job/<job_id>/lengths")
def job_lengths(job_id):
    p = DATA_DIR / job_id / "lengths.json"
    if not p.exists():
        return jsonify({"ready": False})
    return jsonify({"ready": True, "lengths": json.loads(p.read_text())})


@app.route("/api/job/<job_id>/submit_length_filter", methods=["POST"])
def submit_length_filter(job_id):
    job_dir = DATA_DIR / job_id
    min_len = int(request.json["min_len"])
    max_len = int(request.json["max_len"])
    cfg = load_job(job_dir)
    cfg["min_len"] = min_len
    cfg["max_len"] = max_len
    save_job(job_dir, cfg)
    update_status(job_dir, "filter_barcodes", "queued")
    executor.submit(filter_and_barcodes, str(job_dir), min_len, max_len)
    return jsonify({"ok": True})


@app.route("/api/job/<job_id>/barcode_step1")
def barcode_step1(job_id):
    p = DATA_DIR / job_id / "barcode_step1.json"
    if not p.exists():
        return jsonify({"ready": False})
    return jsonify({"ready": True, **json.loads(p.read_text())})


@app.route("/api/job/<job_id>/submit_step4", methods=["POST"])
def submit_step4(job_id):
    selected = request.json.get("selected", [])
    job_dir = DATA_DIR / job_id
    cfg = load_job(job_dir)
    cfg["selected_step4"] = selected
    save_job(job_dir, cfg)
    update_status(job_dir, "final_quant", "queued")
    executor.submit(quantitate_final, str(job_dir), selected)
    return jsonify({"ok": True})


@app.route("/api/job/<job_id>/final_quant")
def final_quant(job_id):
    p = DATA_DIR / job_id / "final_quant.json"
    if not p.exists():
        return jsonify({"ready": False})
    return jsonify({"ready": True, **json.loads(p.read_text())})


@app.route("/api/job/<job_id>/downloads")
def job_downloads(job_id):
    job_dir = DATA_DIR / job_id
    if not job_dir.exists():
        return jsonify({"ready": False}), 404

    files = []
    cfg = load_job(job_dir)
    split_files = cfg.get("split_files")
    filtered_path = filtered_fastq_path(job_dir, cfg)
    paths = [filtered_path] if filtered_path.exists() else []
    if split_files:
        paths.extend(job_dir / name for name in split_files)
    else:
        paths.extend(sorted(job_dir.glob("subset_*.fastq")))

    seen = set()
    for path in paths:
        if not path.exists() or path.name in seen:
            continue
        seen.add(path.name)
        files.append({
            "name": path.name,
            "url": url_for("download_file", job_id=job_id, filename=path.name),
        })
    return jsonify({"ready": True, "files": files})


@app.route("/api/job/<job_id>/submit_step5", methods=["POST"])
def submit_step5(job_id):
    selected_keys = request.json.get("selected", [])
    job_dir = DATA_DIR / job_id
    cfg = load_job(job_dir)
    cfg["selected_step5"] = selected_keys
    save_job(job_dir, cfg)
    update_status(job_dir, "split", "queued")
    executor.submit(split_fastq, str(job_dir), selected_keys)
    return jsonify({"ok": True})


@app.route("/job/<job_id>/download/<path:filename>")
def download_file(job_id, filename):
    p = DATA_DIR / job_id / filename
    if not p.exists():
        return "Not found", 404
    return send_file(p, as_attachment=True)


if __name__ == "__main__":
    app.run(debug=True)

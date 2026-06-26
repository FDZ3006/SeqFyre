"""
app.py - SeqFyre Web Application (Flask)

Antarmuka web untuk pipeline SeqFyre. Pengguna cukup membuka tautan,
mengunggah berkas FASTA/FASTQ/ZIP atau mencoba dataset demo, lalu
memperoleh tabel hasil, 4 grafik EDA, dan unduhan ZIP tanpa instalasi.

Semua keluaran (CSV, PNG, ZIP) dibuat in-memory sehingga tidak
menyentuh disk.
"""

import io
import os
import uuid
from collections import OrderedDict

from flask import Flask, jsonify, render_template, request, send_file
from werkzeug.utils import secure_filename

from seqfyre import Analyzer, Parser, ParserError, build_result_zip

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 16 * 1024 * 1024  # batas unggah 16 MB

# Path dataset demo bawaan
DEMO_DATASET = os.path.join(
    os.path.dirname(__file__), "data", "dataset_16S_rRNA_SeqFyre.fasta"
)

# Cache hasil bertoken, maksimum 32 entri
_RESULT_CACHE = OrderedDict()
_CACHE_LIMIT = 32


def _cache_put(analyzer):
    # Simpan analyzer ke cache, kembalikan token aksesnya
    token = uuid.uuid4().hex
    _RESULT_CACHE[token] = analyzer
    _RESULT_CACHE.move_to_end(token)
    while len(_RESULT_CACHE) > _CACHE_LIMIT:
        _RESULT_CACHE.popitem(last=False)
    return token


def _build_response(analyzer, top_n):
    # Susun payload JSON lengkap untuk frontend
    token = _cache_put(analyzer)
    top = analyzer.top_n(top_n)
    return {
        "token": token,
        "engine": Parser().engine,
        "summary": analyzer.summary(),
        "top": [r.to_row(rank=i) for i, r in enumerate(top, 1)],
        "table": analyzer.table_rows(),
        "plots": analyzer.plots_base64(),
    }


# Rute

@app.route("/")
def index():
    return render_template("index.html")


@app.route("/analyze", methods=["POST"])
def analyze():
    # Terima unggahan berkas, jalankan pipeline, balas JSON
    top_n = _safe_top_n(request.form.get("top_n"))

    if "file" not in request.files or request.files["file"].filename == "":
        return jsonify({"error": "Tidak ada berkas yang diunggah."}), 400

    upload = request.files["file"]
    filename = secure_filename(upload.filename)
    data = upload.read()

    try:
        records = Parser().parse_upload(filename, data)
        analyzer = Analyzer(records)
    except ParserError as exc:
        return jsonify({"error": str(exc)}), 400
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400
    except Exception as exc:
        return jsonify({"error": "Kesalahan tak terduga: " + str(exc)}), 500

    return jsonify(_build_response(analyzer, top_n))


@app.route("/demo", methods=["POST"])
def demo():
    # Jalankan pipeline pada dataset bawaan
    top_n = _safe_top_n(request.form.get("top_n"))
    try:
        records = Parser().parse_path(DEMO_DATASET)
        analyzer = Analyzer(records)
    except (ParserError, FileNotFoundError) as exc:
        return jsonify({"error": "Dataset demo gagal dibaca: " + str(exc)}), 500
    return jsonify(_build_response(analyzer, top_n))


@app.route("/download/<token>")
def download(token):
    # Unduh ZIP berisi CSV dan 4 grafik untuk hasil dengan token tertentu
    analyzer = _RESULT_CACHE.get(token)
    if analyzer is None:
        return jsonify({"error": "Hasil kedaluwarsa. Silakan analisis ulang."}), 404
    zip_bytes = build_result_zip(analyzer)
    return send_file(
        io.BytesIO(zip_bytes),
        mimetype="application/zip",
        as_attachment=True,
        download_name="SeqFyre_hasil.zip",
    )


@app.route("/health")
def health():
    # Endpoint untuk health check
    return jsonify({"status": "ok", "app": "SeqFyre"})


# Helper

def _safe_top_n(raw):
    # Validasi parameter top_n, hanya izinkan 3, 5, atau 10
    try:
        value = int(raw) if raw is not None else 3
    except (TypeError, ValueError):
        value = 3
    if value in (3, 5, 10):
        return value
    return 3


@app.errorhandler(413)
def too_large(_err):
    return jsonify({"error": "Berkas terlalu besar (maksimum 16 MB)."}), 413


if __name__ == "__main__":
    # Untuk pengembangan lokal
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port, debug=True)

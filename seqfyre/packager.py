"""
packager.py
-----------
Mengemas hasil analisis (CSV + 4 grafik PNG) menjadi satu arsip ZIP
yang dibangun sepenuhnya di memori menggunakan BytesIO.
"""

import io
import zipfile

from .analyzer import Analyzer


def build_result_zip(analyzer, prefix="SeqFyre_hasil"):
    # Bentuk arsip ZIP berisi CSV hasil dan keempat grafik EDA
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as zf:
        # 1) CSV hasil analisis
        zf.writestr(prefix + "/hasil_analisis.csv", analyzer.to_csv_bytes())
        # 2) Empat grafik EDA
        for filename, png in analyzer.plots_png_bytes().items():
            zf.writestr(prefix + "/grafik/" + filename, png)
        # 3) README ringkas di dalam arsip
        zf.writestr(prefix + "/README.txt", _readme_text(analyzer))
    buffer.seek(0)
    return buffer.read()


def _readme_text(analyzer):
    # Teks ringkas yang menyertai arsip hasil
    s = analyzer.summary()
    lines = [
        "Hasil Analisis SeqFyre",
        "======================",
        "",
        "Total sekuens     : " + str(s["total_sequences"]),
        "Rata-rata GC       : " + str(s["mean_gc"]) + "%",
        "GC min / maks      : " + str(s["min_gc"]) + "% / " + str(s["max_gc"]) + "%",
        "Termofil / Mesofil : " + str(s["thermophile_count"]) + " / " + str(s["mesophile_count"]),
        "",
        "Isi arsip:",
        "  - hasil_analisis.csv : tabel lengkap seluruh sekuens (terurut GC)",
        "  - grafik/            : 4 grafik EDA (histogram, bar, komposisi, scatter)",
    ]
    return "\n".join(lines)

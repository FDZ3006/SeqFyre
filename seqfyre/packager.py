"""
packager.py
===========
Mengemas hasil analisis (CSV + 4 grafik PNG) menjadi satu arsip ZIP yang
dibangun sepenuhnya di memori (BytesIO).

Pendekatan in-memory dipilih karena Render (dan PaaS sejenis) memakai
ephemeral filesystem — menulis ke disk tidak andal dan bisa hilang.
"""

from __future__ import annotations

import io
import zipfile

from .analyzer import Analyzer


def build_result_zip(analyzer: Analyzer, prefix: str = "SeqFyre_hasil") -> bytes:
    """Bentuk arsip ZIP berisi CSV hasil + keempat grafik EDA.

    Args:
        analyzer: objek Analyzer yang sudah berisi data.
        prefix:   awalan nama berkas di dalam ZIP.

    Returns:
        bytes arsip ZIP, siap dikirim sebagai unduhan oleh Flask.
    """
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as zf:
        # 1) CSV hasil analisis
        zf.writestr(f"{prefix}/hasil_analisis.csv", analyzer.to_csv_bytes())
        # 2) Empat grafik EDA
        for filename, png in analyzer.plots_png_bytes().items():
            zf.writestr(f"{prefix}/grafik/{filename}", png)
        # 3) README kecil di dalam arsip
        zf.writestr(f"{prefix}/README.txt", _readme_text(analyzer))
    buffer.seek(0)
    return buffer.read()


def _readme_text(analyzer: Analyzer) -> str:
    """Teks ringkas yang menyertai arsip hasil."""
    s = analyzer.summary()
    return (
        "Hasil Analisis SeqFyre\n"
        "======================\n\n"
        f"Total sekuens     : {s['total_sequences']}\n"
        f"Rata-rata GC       : {s['mean_gc']}%\n"
        f"GC min / maks      : {s['min_gc']}% / {s['max_gc']}%\n"
        f"Termofil / Mesofil : {s['thermophile_count']} / {s['mesophile_count']}\n\n"
        "Isi arsip:\n"
        "  - hasil_analisis.csv : tabel lengkap seluruh sekuens (terurut GC)\n"
        "  - grafik/            : 4 grafik EDA (histogram, bar, komposisi, scatter)\n"
    )

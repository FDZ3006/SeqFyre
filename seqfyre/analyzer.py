"""
analyzer.py
===========
Kelas Analyzer: otak pipeline. Menerima List[SequenceRecord] lalu:
    - mengurutkan berdasarkan GC content (memenuhi syarat WAJIB),
    - mengambil Top-N sekuens,
    - menghitung statistik agregat & frekuensi nukleotida total (Dictionary),
    - membuat CSV secara in-memory (ramah filesystem ephemeral Render),
    - menghasilkan 4 grafik EDA sebagai PNG base64 (dikirim ke frontend).

Semua keluaran biner (CSV, PNG) dibuat di memori (BytesIO/StringIO) sehingga
tidak menyentuh disk — penting untuk deployment di Render.
"""

from __future__ import annotations

import base64
import csv
import io

import matplotlib

# Backend non-interaktif: wajib di server tanpa display (headless).
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .models import SequenceRecord  # noqa: E402

# Palet warna konsisten dengan tema "fyre" (api/termal).
COLOR_THERMO = "#e8590c"   # oranye-api untuk termofil
COLOR_MESO = "#1c7ed6"     # biru untuk mesofil
COLOR_BASES = {
    "A": "#f08c00",
    "T": "#e8590c",
    "G": "#1c7ed6",
    "C": "#0c8599",
}


class Analyzer:
    """Mengelola koleksi SequenceRecord dan menurunkan semua hasil analisis."""

    def __init__(self, records: list[SequenceRecord], gc_threshold: float = 57.0):
        if not records:
            raise ValueError("Analyzer membutuhkan minimal satu sekuens.")
        # 'records' disimpan sebagai List — struktur data inti yang disyaratkan.
        self.records: list[SequenceRecord] = list(records)
        self.gc_threshold = gc_threshold

    # ------------------------------------------------------------------ #
    # Pengurutan & seleksi
    # ------------------------------------------------------------------ #
    def sort_by_gc(self, descending: bool = True) -> list[SequenceRecord]:
        """Urutkan sekuens berdasarkan GC content.

        Memakai Timsort bawaan Python (sorted) dengan key lambda pada properti
        gc_content. Kompleksitas O(M log M), M = jumlah sekuens.
        """
        return sorted(self.records, key=lambda r: r.gc_content, reverse=descending)

    def top_n(self, n: int = 3) -> list[SequenceRecord]:
        """Ambil N sekuens dengan GC tertinggi (array slicing pada List)."""
        n = max(1, min(n, len(self.records)))
        return self.sort_by_gc(descending=True)[:n]

    # ------------------------------------------------------------------ #
    # Statistik agregat
    # ------------------------------------------------------------------ #
    def total_nucleotide_frequency(self) -> dict[str, int]:
        """Akumulasi frekuensi nukleotida seluruh dataset (Dictionary)."""
        total: dict[str, int] = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0}
        for rec in self.records:
            for base, count in rec.nucleotide_frequency().items():
                total[base] = total.get(base, 0) + count
        return total

    def summary(self) -> dict:
        """Ringkasan statistik untuk panel dashboard."""
        gc_values = [r.gc_content for r in self.records]
        lengths = [r.length for r in self.records]
        thermo = [r for r in self.records if r.gc_content >= self.gc_threshold]
        meso = [r for r in self.records if r.gc_content < self.gc_threshold]
        return {
            "total_sequences": len(self.records),
            "mean_gc": round(sum(gc_values) / len(gc_values), 2),
            "min_gc": round(min(gc_values), 2),
            "max_gc": round(max(gc_values), 2),
            "mean_length": round(sum(lengths) / len(lengths), 1),
            "min_length": min(lengths),
            "max_length": max(lengths),
            "thermophile_count": len(thermo),
            "mesophile_count": len(meso),
            "total_frequency": self.total_nucleotide_frequency(),
        }

    def table_rows(self) -> list[dict]:
        """Daftar baris hasil (terurut, dengan rank) untuk tabel & CSV."""
        ordered = self.sort_by_gc(descending=True)
        return [rec.to_row(rank=i) for i, rec in enumerate(ordered, start=1)]

    # ------------------------------------------------------------------ #
    # Ekspor CSV (in-memory)
    # ------------------------------------------------------------------ #
    def to_csv_string(self) -> str:
        """Hasilkan teks CSV dari seluruh sekuens terurut."""
        rows = self.table_rows()
        buffer = io.StringIO()
        fieldnames = list(rows[0].keys())
        writer = csv.DictWriter(buffer, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
        return buffer.getvalue()

    def to_csv_bytes(self) -> bytes:
        """Versi bytes dari CSV (untuk dimasukkan ke ZIP in-memory)."""
        # utf-8-sig agar Excel membaca karakter dengan benar.
        return self.to_csv_string().encode("utf-8-sig")

    # ------------------------------------------------------------------ #
    # Visualisasi — 4 grafik EDA
    # ------------------------------------------------------------------ #
    @staticmethod
    def _fig_to_base64(fig) -> str:
        """Render Figure matplotlib menjadi data-URI PNG base64."""
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        encoded = base64.b64encode(buf.read()).decode("ascii")
        return f"data:image/png;base64,{encoded}"

    @staticmethod
    def _fig_to_bytes(fig) -> bytes:
        """Render Figure menjadi bytes PNG mentah (untuk ZIP)."""
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        return buf.read()

    def _build_figures(self) -> dict[str, "plt.Figure"]:
        """Bangun keempat figure EDA dan kembalikan sebagai dict nama->Figure."""
        ordered = self.sort_by_gc(descending=True)
        gc_values = [r.gc_content for r in self.records]
        figs: dict[str, plt.Figure] = {}

        # (1) HISTOGRAM distribusi GC -> menonjolkan pola bimodal.
        fig1, ax1 = plt.subplots(figsize=(6.5, 4))
        ax1.hist(gc_values, bins=10, color="#f76707", edgecolor="white")
        ax1.axvline(self.gc_threshold, color="#495057", linestyle="--",
                    linewidth=1.2, label=f"Ambang {self.gc_threshold:.0f}%")
        ax1.set_title("Distribusi GC Content (pola bimodal)")
        ax1.set_xlabel("GC Content (%)")
        ax1.set_ylabel("Jumlah sekuens")
        ax1.legend()
        figs["histogram_gc"] = fig1

        # (2) BAR CHART horizontal GC per sekuens (terurut).
        fig2, ax2 = plt.subplots(figsize=(6.5, max(4, 0.4 * len(ordered))))
        labels = [r.id for r in ordered]
        values = [r.gc_content for r in ordered]
        colors = [COLOR_THERMO if v >= self.gc_threshold else COLOR_MESO
                  for v in values]
        ax2.barh(labels, values, color=colors)
        ax2.invert_yaxis()  # rank 1 di atas
        ax2.set_title("GC Content per Sekuens (terurut)")
        ax2.set_xlabel("GC Content (%)")
        for i, v in enumerate(values):
            ax2.text(v + 0.3, i, f"{v:.1f}", va="center", fontsize=8)
        figs["barchart_gc"] = fig2

        # (3) KOMPOSISI nukleotida (stacked bar A/T/G/C per sekuens).
        fig3, ax3 = plt.subplots(figsize=(7, 4))
        ids = [r.id for r in ordered]
        comps = [r.composition_percent() for r in ordered]
        bottom = [0.0] * len(ordered)
        for base in ("A", "T", "G", "C"):
            heights = [c[base] for c in comps]
            ax3.bar(ids, heights, bottom=bottom, label=base,
                    color=COLOR_BASES[base])
            bottom = [b + h for b, h in zip(bottom, heights)]
        ax3.set_title("Komposisi Nukleotida per Sekuens")
        ax3.set_ylabel("Persentase (%)")
        ax3.set_xticks(range(len(ids)))
        ax3.set_xticklabels(ids, rotation=45, ha="right", fontsize=8)
        ax3.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.18))
        figs["composition"] = fig3

        # (4) SCATTER panjang vs GC -> deteksi outlier & homogenitas ~1500 bp.
        fig4, ax4 = plt.subplots(figsize=(6.5, 4))
        for r in self.records:
            c = COLOR_THERMO if r.gc_content >= self.gc_threshold else COLOR_MESO
            ax4.scatter(r.length, r.gc_content, color=c, s=60,
                        edgecolor="white", zorder=3)
        ax4.set_title("Panjang Sekuens vs GC Content")
        ax4.set_xlabel("Panjang (bp)")
        ax4.set_ylabel("GC Content (%)")
        # Legenda manual untuk dua kelompok ekologis.
        ax4.scatter([], [], color=COLOR_THERMO, label="Termofil")
        ax4.scatter([], [], color=COLOR_MESO, label="Mesofil")
        ax4.legend()
        figs["scatter_len_gc"] = fig4

        return figs

    def plots_base64(self) -> dict[str, str]:
        """Keempat grafik EDA sebagai data-URI base64 (untuk frontend web)."""
        return {name: self._fig_to_base64(fig)
                for name, fig in self._build_figures().items()}

    def plots_png_bytes(self) -> dict[str, bytes]:
        """Keempat grafik EDA sebagai bytes PNG (untuk dimasukkan ke ZIP)."""
        return {f"{name}.png": self._fig_to_bytes(fig)
                for name, fig in self._build_figures().items()}

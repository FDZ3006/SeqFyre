"""
analyzer.py
-----------
Kelas Analyzer: menerima List berisi SequenceRecord lalu:
    - mengurutkan berdasarkan GC content
    - mengambil Top-N sekuens
    - menghitung statistik agregat dan frekuensi nukleotida total
    - membuat CSV secara in-memory
    - menghasilkan 4 grafik EDA sebagai PNG base64
"""

import base64
import csv
import io

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .models import SequenceRecord

# Warna untuk dua kelompok ekologis
COLOR_THERMO = "#e8590c"
COLOR_MESO = "#1c7ed6"
COLOR_BASES = {
    "A": "#f08c00",
    "T": "#e8590c",
    "G": "#1c7ed6",
    "C": "#0c8599",
}


class Analyzer:

    def __init__(self, records, gc_threshold=57.0):
        if not records:
            raise ValueError("Analyzer membutuhkan minimal satu sekuens.")
        # Simpan koleksi sekuens dalam struktur List
        self.records = list(records)
        self.gc_threshold = gc_threshold

    # Pengurutan dan seleksi

    def sort_by_gc(self, descending=True):
        # Urutkan sekuens berdasarkan GC content menggunakan sorted()
        # Timsort bawaan Python, kompleksitas O(M log M)
        return sorted(self.records, key=lambda r: r.gc_content(), reverse=descending)

    def top_n(self, n=3):
        # Ambil N sekuens dengan GC tertinggi menggunakan slicing List
        if n < 1:
            n = 1
        if n > len(self.records):
            n = len(self.records)
        return self.sort_by_gc(descending=True)[:n]

    # Statistik agregat

    def total_nucleotide_frequency(self):
        # Akumulasi frekuensi nukleotida seluruh dataset menggunakan Dictionary
        total = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0}
        for rec in self.records:
            freq = rec.nucleotide_frequency()
            for base in freq:
                total[base] = total.get(base, 0) + freq[base]
        return total

    def summary(self):
        # Ringkasan statistik dataset
        gc_values = [r.gc_content() for r in self.records]
        lengths = [r.get_length() for r in self.records]
        thermo = [r for r in self.records if r.gc_content() >= self.gc_threshold]
        meso = [r for r in self.records if r.gc_content() < self.gc_threshold]
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

    def table_rows(self):
        # Daftar baris hasil terurut dengan rank untuk tabel dan CSV
        ordered = self.sort_by_gc(descending=True)
        rows = []
        for i, rec in enumerate(ordered, start=1):
            rows.append(rec.to_row(rank=i))
        return rows

    # Ekspor CSV

    def to_csv_string(self):
        # Hasilkan teks CSV dari seluruh sekuens terurut
        rows = self.table_rows()
        buffer = io.StringIO()
        fieldnames = list(rows[0].keys())
        writer = csv.DictWriter(buffer, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
        return buffer.getvalue()

    def to_csv_bytes(self):
        # Versi bytes dari CSV, pakai utf-8-sig agar Excel membaca dengan benar
        return self.to_csv_string().encode("utf-8-sig")

    # Visualisasi - 4 grafik EDA

    def _fig_to_base64(self, fig):
        # Render Figure matplotlib menjadi data-URI PNG base64
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        encoded = base64.b64encode(buf.read()).decode("ascii")
        return "data:image/png;base64," + encoded

    def _fig_to_bytes(self, fig):
        # Render Figure menjadi bytes PNG mentah untuk ZIP
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        return buf.read()

    def _build_figures(self):
        # Bangun keempat figure EDA dan kembalikan sebagai dictionary
        ordered = self.sort_by_gc(descending=True)
        gc_values = [r.gc_content() for r in self.records]
        figs = {}

        # 1) Histogram distribusi GC - menunjukkan pola bimodal
        fig1, ax1 = plt.subplots(figsize=(6.5, 4))
        ax1.hist(gc_values, bins=10, color="#f76707", edgecolor="white")
        ax1.axvline(
            self.gc_threshold,
            color="#495057",
            linestyle="--",
            linewidth=1.2,
            label="Ambang " + str(int(self.gc_threshold)) + "%"
        )
        ax1.set_title("Distribusi GC Content (pola bimodal)")
        ax1.set_xlabel("GC Content (%)")
        ax1.set_ylabel("Jumlah sekuens")
        ax1.legend()
        figs["histogram_gc"] = fig1

        # 2) Bar chart horizontal GC per sekuens terurut
        fig2, ax2 = plt.subplots(figsize=(6.5, max(4, 0.4 * len(ordered))))
        labels = [r.id for r in ordered]
        values = [r.gc_content() for r in ordered]
        colors = []
        for v in values:
            if v >= self.gc_threshold:
                colors.append(COLOR_THERMO)
            else:
                colors.append(COLOR_MESO)
        ax2.barh(labels, values, color=colors)
        ax2.invert_yaxis()
        ax2.set_title("GC Content per Sekuens (terurut)")
        ax2.set_xlabel("GC Content (%)")
        for i, v in enumerate(values):
            ax2.text(v + 0.3, i, str(round(v, 1)), va="center", fontsize=8)
        figs["barchart_gc"] = fig2

        # 3) Komposisi nukleotida stacked bar A/T/G/C per sekuens
        fig3, ax3 = plt.subplots(figsize=(7, 4))
        ids = [r.id for r in ordered]
        comps = [r.composition_percent() for r in ordered]
        bottom = [0.0] * len(ordered)
        for base in ("A", "T", "G", "C"):
            heights = [c[base] for c in comps]
            ax3.bar(ids, heights, bottom=bottom, label=base, color=COLOR_BASES[base])
            bottom = [b + h for b, h in zip(bottom, heights)]
        ax3.set_title("Komposisi Nukleotida per Sekuens")
        ax3.set_ylabel("Persentase (%)")
        ax3.set_xticks(range(len(ids)))
        ax3.set_xticklabels(ids, rotation=45, ha="right", fontsize=8)
        ax3.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.18))
        figs["composition"] = fig3

        # 4) Scatter panjang vs GC - deteksi outlier dan homogenitas panjang
        fig4, ax4 = plt.subplots(figsize=(6.5, 4))
        for r in self.records:
            if r.gc_content() >= self.gc_threshold:
                c = COLOR_THERMO
            else:
                c = COLOR_MESO
            ax4.scatter(r.get_length(), r.gc_content(), color=c, s=60, edgecolor="white", zorder=3)
        ax4.set_title("Panjang Sekuens vs GC Content")
        ax4.set_xlabel("Panjang (bp)")
        ax4.set_ylabel("GC Content (%)")
        ax4.scatter([], [], color=COLOR_THERMO, label="Termofil")
        ax4.scatter([], [], color=COLOR_MESO, label="Mesofil")
        ax4.legend()
        figs["scatter_len_gc"] = fig4

        return figs

    def plots_base64(self):
        # Keempat grafik EDA sebagai data-URI base64 untuk frontend web
        result = {}
        for name, fig in self._build_figures().items():
            result[name] = self._fig_to_base64(fig)
        return result

    def plots_png_bytes(self):
        # Keempat grafik EDA sebagai bytes PNG untuk dimasukkan ke ZIP
        result = {}
        for name, fig in self._build_figures().items():
            result[name + ".png"] = self._fig_to_bytes(fig)
        return result

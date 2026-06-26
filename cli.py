#!/usr/bin/env python3
"""
cli.py - SeqFyre Command-Line Pipeline

Skrip mandiri yang memenuhi seluruh persyaratan wajib mini project:
    - Membaca berkas FASTA atau FASTQ
    - Menyimpan data dalam struktur List
    - Menghitung frekuensi nukleotida menggunakan Dictionary
    - Mengurutkan sekuens berdasarkan GC content
    - Menampilkan 3 sekuens terbaik
    - Visualisasi grafik berdasarkan nilai GC
    - Menulis hasil ke berkas CSV

Cara pakai:
    python cli.py data/dataset_16S_rRNA_SeqFyre.fasta
    python cli.py data/dataset_16S_rRNA_SeqFyre.fasta --top 5 --outdir hasil

Keluaran:
    <outdir>/hasil_analisis.csv
    <outdir>/grafik/*.png
"""

import argparse
import os
import sys

from seqfyre import Analyzer, Parser, ParserError


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="SeqFyre - pipeline analisis GC content (CLI)."
    )
    p.add_argument("input", help="Berkas FASTA/FASTQ/ZIP yang dianalisis.")
    p.add_argument("--top", type=int, default=3,
                   help="Jumlah sekuens terbaik yang ditampilkan (default 3).")
    p.add_argument("--outdir", default="hasil",
                   help="Direktori keluaran (default 'hasil').")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    # 1. Baca berkas FASTA/FASTQ, simpan ke List
    parser = Parser()
    print("[*] Mesin parsing : " + parser.engine)
    try:
        records = parser.parse_path(args.input)
    except (ParserError, FileNotFoundError) as exc:
        print("[!] Gagal membaca berkas: " + str(exc))
        return 1
    print("[*] Sekuens terbaca: " + str(len(records)) + " (disimpan dalam List)")

    # 2. Hitung frekuensi nukleotida menggunakan Dictionary
    print("\n[*] Contoh frekuensi nukleotida (Dictionary) - sekuens pertama:")
    print("    " + records[0].id + ": " + str(records[0].nucleotide_frequency()))

    # 3. Urutkan berdasarkan GC content, tampilkan Top-N
    analyzer = Analyzer(records)
    top = analyzer.top_n(args.top)
    print("\n=== " + str(len(top)) + " SEKUENS TERBAIK (GC tertinggi) ===")
    print(f"{'RANK':<5}{'GC%':<8}{'ACCESSION':<14}{'KELAS':<10}ORGANISME")
    for i, rec in enumerate(top, 1):
        print(f"{i:<5}{rec.gc_content():<8.2f}{rec.id:<14}"
              f"{rec.classify_by_gc():<10}{rec.organism}")

    # 4. Siapkan direktori keluaran
    os.makedirs(os.path.join(args.outdir, "grafik"), exist_ok=True)

    # 5. Tulis hasil ke CSV
    csv_path = os.path.join(args.outdir, "hasil_analisis.csv")
    with open(csv_path, "w", encoding="utf-8-sig", newline="") as fh:
        fh.write(analyzer.to_csv_string())
    print("\n[*] CSV ditulis   : " + csv_path)

    # 6. Visualisasi grafik berdasarkan GC
    for filename, png in analyzer.plots_png_bytes().items():
        path = os.path.join(args.outdir, "grafik", filename)
        with open(path, "wb") as fh:
            fh.write(png)
    print("[*] 4 grafik EDA  : " + os.path.join(args.outdir, "grafik") + "/")

    # Ringkasan
    s = analyzer.summary()
    print("\n=== RINGKASAN ===")
    print("    Rata-rata GC       : " + str(s["mean_gc"]) + "%")
    print("    Termofil / Mesofil : " + str(s["thermophile_count"]) + " / " + str(s["mesophile_count"]))
    print("[v] Selesai.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
